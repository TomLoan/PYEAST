# Copyright CSIRO 2025. Thomas Loan 
# See LICENSE for full GpLv2 license. 

# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License or 
# (at your option) any later version. 

# This program is distributed in the hope that it will be useful;, 
# but WITHOUT ANY WARRENTY; without even the implied warranty of 
# MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details. 

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc., 
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 

# ===========================================================================

# src/pyeast/cli/main.py
""""
CLI program for PYEAST program. 
"""

# ===========================================================================

import os
import io
import datetime
from pathlib import Path
from typing import Optional, List
import click
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from prompt_toolkit import PromptSession
from prompt_toolkit.completion import WordCompleter
from PIL import Image
from Bio import SeqIO
from datetime import datetime

from ..core.tar import TARDesigner
from ..core.integration import IntegrationDesigner
from ..core.deletion import DeletionDesigner
from ..core.replace import ReplaceDesigner
from ..core.batch import BatchDesigner

from .. utils.visualisation import visualise_genbank, save_figure

console = Console()


def get_output_prefix() -> str:
    """Get output prefix from user - just the name, not path"""
    session = PromptSession()
    output_dir = Path("output")
    
    # Ensure output directory exists
    output_dir.mkdir(exist_ok=True)
    
    while True:
        try:
            console.print("\n[blue]Enter a name prefix for your output files[/blue]")
            console.print("[dim]Files will be saved in the output directory[/dim]")
            user_input = session.prompt("Name prefix: ").strip()
            
            # Basic validation
            if not user_input:
                console.print("[red]Please enter a name[/red]")
                continue
                
            # Remove any path separators for safety
            safe_name = "".join(c for c in user_input if c not in r'\/.')
            
            if safe_name != user_input:
                console.print(f"[yellow]Name adjusted to: {safe_name}[/yellow]")
            
            # Check if files would be overwritten
            potential_files = [
                output_dir / f"{safe_name}.gb",
                output_dir / f"{safe_name}_map.png",
                output_dir / f"{safe_name}_instructions.tsv",
                output_dir / f"{safe_name}_missing_primers.tsv",
                output_dir / f"{safe_name}_all_primers.tsv",
                output_dir / f"{safe_name}_screening_primers.tsv",
                output_dir / f"{safe_name}.fasta"
            ]
            
            existing = [f for f in potential_files if f.exists()]
            if existing:
                console.print("[yellow]Warning: The following files already exist:[/yellow]")
                for f in existing:
                    console.print(f"  - {f.name}")
                if not click.confirm("Overwrite these files?"):
                    continue
            
            return str(output_dir / safe_name)
                
        except KeyboardInterrupt:
            if click.confirm("\nCancel output name selection?"):
                raise click.Abort()

def handle_machine_instructions(designer: BatchDesigner, timestamp: str) -> None:
    """Ask user if they want to generate machine instructions and handle the response.
    
    Args:
        designer: BatchDesigner instance with completed assembly instructions
        timestamp: Current timestamp for file naming
    """
    if click.confirm("\nWould you like to generate machine instructions for liquid handling?"):
        # For now we only have epMotion, but structure allows easy expansion
        machines = ['epMotion']
        
        print("\nAvailable liquid handling machines:")
        for i, machine in enumerate(machines, 1):
            print(f"{i}. {machine}")
            
        if len(machines) == 1:
            print("\nCurrently only epMotion is supported.")
            if click.confirm("Generate epMotion instructions?"):
                try:
                    output_file = designer.generate_epmotion_instructions(timestamp)
                    print(f"\nMachine instructions saved to: {output_file}")
                except Exception as e:
                    print(f"\nError generating machine instructions: {str(e)}")


def get_component_dir() -> Path:
    """Get component directory with simple autocompletion"""
    # Default to src/pyeast/data/component libraries
    base_dir = Path("data/component libraries")
    console.print(base_dir)
    
    if not base_dir.exists():
        console.print("[red]Error: Default components directory not found[/red]")
        raise click.Abort()
        
    # Create completer from subdirectories
    subdirs = [d.name for d in base_dir.iterdir() if d.is_dir()]
    dir_completer = WordCompleter(subdirs, ignore_case=True)
    
    session = PromptSession()
    
    # Show available directories
    table = Table(title="Available Component Directories")
    table.add_column("Name", style="cyan")
    for subdir in subdirs:
        table.add_row(subdir)
    console.print(table)
    
    while True:
        try:
            user_input = session.prompt(
                "\nSelect components directory (TAB for suggestions): ",
                completer=dir_completer
            )
            
            selected_dir = base_dir / user_input
            if selected_dir.exists() and selected_dir.is_dir():
                return selected_dir
            else:
                console.print("[red]Invalid directory selection[/red]")
                
        except KeyboardInterrupt:
            if click.confirm("\nCancel directory selection?"):
                raise click.Abort()

def run_tar_interactive_mode(designer: TARDesigner):
    """Run TAR design in interactive mode"""
    try:
        # Get components directory
        components_dir = get_component_dir()
        
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            # Load and process sequences
            task_id = progress.add_task("Loading sequences...", total=None)
            sequences = designer.load_and_get_sequences(components_dir)
            progress.update(task_id, completed=True)
            
            if not sequences:
                console.print("[red]No sequences found in directory[/red]")
                return
            
        # Display sequences and get assembly order
        designer.print_sequence_grid(sequences)
        assembly_order = designer.get_assembly_order(sequences)
        
        if not assembly_order:
            return
        
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            # Set assembly order and design primers
            task_id = progress.add_task("Designing primers...", total=None)
            designer.set_assembly_order(assembly_order)
            designer.design_tar_primers()
            progress.update(task_id, completed=True)
            
            # Check primer locations
            task_id = progress.add_task("Checking primer locations...", total=None)
            designer.check_primer_locations(Path("data/primers"))
            progress.update(task_id, completed=True)
            
            # Find templates
            task_id = progress.add_task("Finding templates...", total=None)
            designer.find_templates(Path("data/templates"))
            progress.update(task_id, completed=True)
            
            # Rationalize selections
            task_id = progress.add_task("Rationalizing selections...", total=None)
            designer.rationalize_selections()
            progress.update(task_id, completed=True)
            
            # Generate instructions
            task_id = progress.add_task("Generating instructions...", total=None)
            instructions = designer.write_instructions()
            progress.update(task_id, completed=True)
        
        # Display instructions and confirm
        designer.display_instructions(instructions)
        
        if not click.confirm("\nProceed with assembly?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return
        
        output_prefix = get_output_prefix()
        designer.console.print(output_prefix)
       
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            # Generate assembly
            task_id = progress.add_task("Generating assembly...", total=None)
            assembly = designer.create_assembly()
            assembly.name = output_prefix.split("\\")[1]
            progress.update(task_id, completed=True)
            
            # Get output prefix after design is confirmed
        
        
        # Save outputs
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)
            
            try:
                # Save GenBank file
                SeqIO.write(assembly, f"{output_prefix}.gb", "genbank")
                
                # Save all primers
                with open(f"{output_prefix}_all_primers.tsv", 'w') as f:
                    f.write("Name\tSequence\n")
                    for name, primer in designer.primers.items():
                        f.write(f"{name}\t{primer}\n")
                
                # Save missing primers if any
                if designer.missing_primers:
                    with open(f"{output_prefix}_missing_primers.tsv", 'w') as f:
                        f.write("Name\tSequence\n")
                        for name, info in designer.missing_primers.items():
                            f.write(f"{name}\t{info[0]['sequence']}\n")
                
                # Save instructions
                with open(f"{output_prefix}_instructions.tsv", 'w') as f:
                    f.write("Part\tF_Plate\tF_Name\tF_Well\tR_Plate\tR_Name\tR_Well\tTemplate\tSize\n")
                    for row in instructions:
                        f.write("\t".join(map(str, row)) + "\n")
                
                # Generate and save map
                img_data, fig = visualise_genbank(f"{output_prefix}.gb")
                save_figure(fig, f"{output_prefix}_map.png")
                
                progress.update(task_id, completed=True)
                
                # Show summary of saved files
                console.print("\n[bold green]Files saved:[/bold green]")
                console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
                console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
                console.print(f"[green]Assembly instructions: {output_prefix}_instructions.tsv[/green]")
                console.print(f"[green]All primers: {output_prefix}_all_primers.tsv[/green]")
                if designer.missing_primers:
                    console.print(f"[green]Missing primers: {output_prefix}_missing_primers.tsv[/green]")
                
                # Show map
                img = Image.open(io.BytesIO(img_data.getvalue()))
                img.show()
                
            except Exception as e:
                console.print(f"[red]Error saving files: {str(e)}[/red]")
                raise
            
            console.print("\n[bold green]✓[/bold green] Plasmid design complete!")
            
    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise

def run_integration_interactive_mode(designer: IntegrationDesigner):
    """Run integration design in interactive mode"""
    try:
        # Get components directory
        components_dir = get_component_dir()
        
        # Load sequences with progress bar
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Loading sequences...", total=None)
            designer.load_sequences(components_dir)
            progress.update(task_id, completed=True)
            
        if not designer.components or not designer.int_sites:
            console.print("[red]No sequences found[/red]")
            return

        # Get user selections (interactive, no progress bar needed)
        designer.get_assembly_selections()
        if not designer.assembly_sequences:
            return

        # Design process with progress bar
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            # Design primers
            task_id = progress.add_task("Designing primers...", total=None)
            designer.primers = designer.design_integration_primers()
            progress.update(task_id, completed=True)

            # Find existing primers
            task_id = progress.add_task("Checking primer locations...", total=None)
            designer.check_primer_locations(Path("data/primers"))
            progress.update(task_id, completed=True)

            # Find templates
            task_id = progress.add_task("Finding templates...", total=None)
            designer.find_templates(Path("data/templates"))
            progress.update(task_id, completed=True)
            
            # Rationalize selections
            task_id = progress.add_task("Rationalizing selections...", total=None)
            designer.rationalize_selections()            
            progress.update(task_id, completed=True)            
            
            # Generate instructions
            task_id = progress.add_task("Generating instructions...", total=None)
            instructions = designer.write_instructions()
            progress.update(task_id, completed=True)
        
        # Display instructions and confirm
        designer.display_instructions(instructions)
            

        if not click.confirm("\nProceed with assembly?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return
        
        # Get output prefix and save results
        output_prefix = get_output_prefix()

        # Create assembly
        task_id = progress.add_task("Creating assembly...", total=None)
        assembly = designer.create_linear_assembly()
        assembly.name = output_prefix.split("\\")[1]
        progress.update(task_id, completed=True)
        
        try:
                # Save GenBank file
                SeqIO.write(assembly, f"{output_prefix}.gb", "genbank")
                
                # Save all primers
                with open(f"{output_prefix}_all_primers.tsv", 'w') as f:
                    f.write("Name\tSequence\n")
                    for name, primer in designer.primers.items():
                        f.write(f"{name}\t{primer}\n")
                
                # Save missing primers if any
                if designer.missing_primers:
                    with open(f"{output_prefix}_missing_primers.tsv", 'w') as f:
                        f.write("Name\tSequence\n")
                        for name, info in designer.missing_primers.items():
                            f.write(f"{name}\t{info[0]['sequence']}\n")
                
                # Save instructions
                with open(f"{output_prefix}_instructions.tsv", 'w') as f:
                    f.write("Part\tF_Plate\tF_Name\tF_Well\tR_Plate\tR_Name\tR_Well\tTemplate\tSize\n")
                    for row in instructions:
                        f.write("\t".join(map(str, row)) + "\n")
                
                # Generate and save map
                img_data, fig = visualise_genbank(f"{output_prefix}.gb")
                save_figure(fig, f"{output_prefix}_map.png")
                
                progress.update(task_id, completed=True)
                
                # Show summary of saved files
                console.print("\n[bold green]Files saved:[/bold green]")
                console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
                console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
                console.print(f"[green]Assembly instructions: {output_prefix}_instructions.tsv[/green]")
                console.print(f"[green]All primers: {output_prefix}_all_primers.tsv[/green]")
                if designer.missing_primers:
                    console.print(f"[green]Missing primers: {output_prefix}_missing_primers.tsv[/green]")
                
                # Show map
                img = Image.open(io.BytesIO(img_data.getvalue()))
                img.show()
                
        except Exception as e:
            console.print(f"[red]Error saving files: {str(e)}[/red]")
            raise
                

        console.print("\n[bold green]✓[/bold green] Integration design complete!")
        
    except KeyboardInterrupt:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()
    
def run_deletion_interactive_mode(designer: DeletionDesigner):
    """Run deletion design in interactive mode."""
    try:
        console = Console()
        
        # Get target sequence from user
        while True:
            console.print("\n[blue]Enter the DNA sequence you want to delete:[/blue]")
            console.print("[dim]Use only A, T, G, and C[/dim]")
            sequence = designer.session.prompt("Sequence: ").upper().strip()
            
            if not set(sequence).issubset({'A', 'T', 'G', 'C'}):
                console.print("[red]Invalid DNA sequence. Please use only A, T, G, and C.[/red]")
                continue
            if len(sequence) < designer.downstream_homology_len:
                console.print("[red]Target sequence is shorter than the homology lenghts. Adjust this parameter with the --downstream_homology_len option") 
                continue  
            
            #all chacks passed
            designer.target_sequence = sequence
            break
        # Find target in genome
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Searching genome for target sequence...", total=None)
            designer.genome_location = designer.find_target_sequence(
                designer.genome_file,
                designer.target_sequence
            )
            progress.update(task_id, completed=True)
            
        if not designer.genome_location:
            console.print("[red]Target sequence not found in genome.[/red]")
            return
            
        chrom_id, start, end, orientation = designer.genome_location
        console.print(f"\n[green]Target sequence found:[/green]")
        console.print(f"Chromosome: {chrom_id}")
        console.print(f"Position: {start}-{end}")
        console.print(f"Orientation: {orientation}")
        
        if not click.confirm("\nProceed with deletion design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return
            
        # Create deletion cassette
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing deletion cassette...", total=None)
            designer.deletion_cassette = designer.make_deletion_cassette(
                designer.genome_file,
                designer.ura3_file
            )
            progress.update(task_id, completed=True)
            
            # Design screening primers
            task_id = progress.add_task("Designing screening strategy...", total=None)
            forward_primer, reverse_primer, product_sizes = designer.design_screening_strategy(
                designer.genome_file
            )
            progress.update(task_id, completed=True)
        
        # Show screening strategy
        console.print("\n[bold cyan]Screening Strategy:[/bold cyan]")
        table = Table()
        table.add_column("Stage", style="cyan")
        table.add_column("Product Size", style="green")
        
        table.add_row("Parent strain", str(product_sizes['parent']))
        table.add_row("After transformation", str(product_sizes['deletion_intermediate']))
        table.add_row("Final deletion", str(product_sizes['final_deletion']))
        
        console.print(table)
        console.print("\n[bold cyan]Screening Primers:[/bold cyan]")
        console.print(f"Forward: {forward_primer}")
        console.print(f"Reverse: {reverse_primer}")
        
        if not click.confirm("\nSave deletion design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return
            
        # Get output prefix
        output_prefix = get_output_prefix()

        #rename the deletion cassette to the user input 
        designer.deletion_cassette.name = output_prefix.split("\\")[1]
        
        # Save results
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)

            if not designer.deletion_cassette:
                raise ValueError("No deletion cassette designed. Run make_deletion_cassette first.")
            
            if not designer.screening_primers:
                raise ValueError("No screening primers designed. Run design_screening_strategy first.")
                
            # Save cassette as GenBank and FASTA
            SeqIO.write(designer.deletion_cassette, f"{output_prefix}.gb", "genbank")
            SeqIO.write(designer.deletion_cassette, f"{output_prefix}.fasta", "fasta")
            
            # Save primers
            forward_primer, reverse_primer = designer.screening_primers
            with open(f"{output_prefix}_screening_primers.tsv", 'w') as f:
                f.write(f"{output_prefix}_ScreenF\t{forward_primer}\n")
                f.write(f"{output_prefix}_ScreenR\t{reverse_primer}")
                
            # Generate and save map
            img_data, fig = visualise_genbank(f"{output_prefix}.gb")
            save_figure(fig, f"{output_prefix}_map.png")

            # Show map
            img = Image.open(io.BytesIO(img_data.getvalue()))
            img.show()
            
            # Print summary
            console.print("\n[bold green]Files saved:[/bold green]")
            console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
            console.print(f"[green]FASTA file: {output_prefix}.fasta[/green]")
            console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
            console.print(f"[green]Screening primers: {output_prefix}_screening_primers.tsv[/green]")

            
        console.print("\n[bold green]✓[/bold green] Deletion design complete!")
        
    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()

def run_replace_interactive_mode(designer: ReplaceDesigner):
    """Run replace design in interactive mode."""
    try:
        console = Console()
        
        # Get target sequence from user
        while True:
            console.print("\n[blue]Enter the DNA sequence you want to replace:[/blue]")
            console.print("[dim]Use only A, T, G, and C[/dim]")
            sequence = designer.session.prompt("Sequence: ").upper().strip()
            
            #check for valid DNA sequence
            if not set(sequence).issubset({'A', 'T', 'G', 'C'}):
                console.print("[red]Invalid DNA sequence. Please use only A, T, G, and C.[/red]")
                continue

            #min target length limited by the size of homology lengths
            min_target_len = min(designer.downstream_homology_len, designer.upstream_homology_len)              

            if len(sequence) < min_target_len:
                console.print("[red]Target sequence is shorter than the homology lenghts. Adjust these parameters with the --downstream_homology_len and --upstream_homology_len options") 
                continue 
            #all checks passed
            designer.target_sequence = sequence
            break
                            
        # Find target in genome
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Searching genome for target sequence...", total=None)
            designer.genome_location = designer.find_target_sequence(
                designer.genome_file,
                designer.target_sequence
            )
            progress.update(task_id, completed=True)
            
        if not designer.genome_location:
            console.print("[red]Target sequence not found in genome.[/red]")
            return
            
        chrom_id, start, end, orientation = designer.genome_location
        console.print(f"\n[green]Target sequence found:[/green]")
        console.print(f"Chromosome: {chrom_id}")
        console.print(f"Position: {start}-{end}")
        console.print(f"Orientation: {orientation}")
        
        # Get components directory
        components_dir = get_component_dir()
        
        # Get replacement sequence selection
        designer.replacement_sequence = designer.get_replacement_selection(components_dir)
        if not designer.replacement_sequence:
            console.print("[yellow]No replacement sequence selected. Design cancelled.[/yellow]")
            return
            
        # Get marker position
        while True:
            console.print("\n[blue]Do you want the URA3 marker upstream or downstream of the replacement?[/blue]")
            position = designer.session.prompt("Position (up/down): ", complete_in_thread= True, completer = WordCompleter(['up', 'down'])).lower().strip()
            if position in ['up', 'down']:
                designer.marker_position = 'upstream' if position == 'up' else 'downstream'
                break
            else:
                console.print("[red]Invalid position. Please enter 'up' or 'down'.[/red]")
            
        # Create replacement cassette
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing replacement cassette...", total=None)
            designer.replacement_cassette = designer.make_replacement_cassette(
                designer.genome_file,
                designer.ura3_file
            )
            progress.update(task_id, completed=True)
            
            # Design screening primers
            task_id = progress.add_task("Designing screening strategy...", total=None)
            forward_primer, reverse_primer, product_sizes = designer.design_screening_strategy(
                designer.genome_file
            )
            progress.update(task_id, completed=True)
        
        # Show screening strategy
        console.print("\n[bold cyan]Screening Strategy:[/bold cyan]")
        table = Table()
        table.add_column("Stage", style="cyan")
        table.add_column("Product Size", style="green")
        
        table.add_row("Parent strain", str(product_sizes['parent']))
        table.add_row("After transformation", str(product_sizes['replacement_intermediate']))
        table.add_row("Final replacement", str(product_sizes['final_replacement']))
        
        console.print(table)
        console.print("\n[bold cyan]Screening Primers:[/bold cyan]")
        console.print(f"Forward: {forward_primer}")
        console.print(f"Reverse: {reverse_primer}")
        
        if not click.confirm("\nSave replacement design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return
            
        # Get output prefix
        output_prefix = get_output_prefix()

        designer.replacement_cassette.name = output_prefix.split("\\")[1]
        
        # Save results
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)
            
            # Save cassette as GenBank and FASTA
            SeqIO.write(designer.replacement_cassette, f"{output_prefix}.gb", "genbank")
            SeqIO.write(designer.replacement_cassette, f"{output_prefix}.fasta", "fasta")
            
            # Save primers
            forward_primer, reverse_primer = designer.screening_primers
            with open(f"{output_prefix}_screening_primers.tsv", 'w') as f:
                f.write(f"{output_prefix}_ScreenF\t{forward_primer}\n")
                f.write(f"{output_prefix}_ScreenR\t{reverse_primer}")
                
            # Generate and save map
            img_data, fig = visualise_genbank(f"{output_prefix}.gb")
            save_figure(fig, f"{output_prefix}_map.png")
            
            # Show map
            img = Image.open(io.BytesIO(img_data.getvalue()))
            img.show()
            
            # Print summary
            console.print("\n[bold green]Files saved:[/bold green]")
            console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
            console.print(f"[green]FASTA file: {output_prefix}.fasta[/green]")
            console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
            console.print(f"[green]Screening primers: {output_prefix}_screening_primers.tsv[/green]")
            
            progress.update(task_id, completed=True)
            
        console.print("\n[bold green]✓[/bold green] Replacement design complete!")
        
    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()

def run_batch_interactive_mode(designer: BatchDesigner):
    """Run batch assembly design in interactive mode"""
    try:
        # Load available constructs
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Loading available constructs...", total=None)
            try:
                designer.load_constructs()
            except ValueError as e:
                console.print(f"[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)
       
        # Get construct selections
        try:
            designer.get_selections()
        except KeyboardInterrupt:
            console.print("\n[yellow]Operation cancelled[/yellow]")
            return
        
        # Validate selections
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Validating constructs...", total=None)
            try:
                designer.validate_constructs()
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)
            
        # Process constructs
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Processing selected constructs...", total=None)
            try:
                designer.process_selected_constructs()
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        # Find primers and templates
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Finding primers and templates...", total=None)
            try:
                designer.find_primers_and_templates()
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        # Organize PCR batches
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Organizing PCR batches...", total=None)
            try:
                designer.organize_pcr_batches()
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        if not click.confirm("\nProceed with generating instructions?"):
            console.print("[yellow]Operation cancelled[/yellow]")
            return

        # Get output prefix
        output_prefix = get_output_prefix()

        # Generate instructions
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            # Generate human-readable instructions
            task_id = progress.add_task("Generating instructions...", total=None)
            try:
                designer.generate_human_instructions(output_prefix)
                designer.generate_assembly_groups(output_prefix)
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        # Optionally generate machine instructions
        if click.confirm("\nGenerate machine instructions for epMotion?"):
            task_id = progress.add_task("Generating machine instructions...", total=None)
            try:
                timestamp = datetime.now().strftime("%H-%M-%d-%b-%Y").upper()
                designer.generate_epmotion_instructions(timestamp)
                designer.generate_machine_assembly_instructions('epmotion', timestamp)
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        console.print("\n[bold green]✓[/bold green] Batch design complete!")
        
    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()
        
    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()
    
@click.group()
def cli():
    """PYeast: Python tools for yeast genetic engineering
    
    Created by Tom Loan"""
    pass

@cli.command()
@click.option('--homology_length', 
            type = int, 
            default = 25, 
            help='Length of homology added to each primer (default: 25). Note the resulting homologous junction will be 2x this number')
def tar(homology_length):
    """Design TAR cloning experiments in Saccharomyces cerevisiae
    \b\n
    Mechamism:

    \b
    5'==Component1==3' + 5'==Component2==3' + ...	PCR Products
		||
	**transformation**
		\\/
	Assembled plasmid                               Plasmid DNA
    """
    try:
        designer = TARDesigner(homology_length = homology_length)
        run_tar_interactive_mode(designer)
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()
    
@cli.command()
@click.option('--homology_length', 
              type=int, 
              default=25,
              help='Length of homology for each junction (default: 25) Note the resulting homologous junction will be 2x this number')
def integrate(homology_length):
    """Design genomic integrations in Saccharomyces cerevisiae
    \b\n
    Mechamism:

    \b
    upstream_seq + Component1 + Component2 + ... + downstream_seq	PCR products
	      X                                          X
    upstream_seq==================================downsteam_seq	        gDNA
		||
	**transformation**
		\\/
    upstream+component1==component1==...==downstream_seq                gDNA    
    """
    
    try:
        designer = IntegrationDesigner(homology_length=homology_length)
        run_integration_interactive_mode(designer)      
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()

@cli.command()
@click.option('--upstream_homology_len', 
              type=int, 
              default=300,
              help='Length of upstream homology for recombination (default: 300)')
@click.option('--downstream_homology_len',
              type=int,
              default=200,
              help='Length of downstream homology for recombination (default: 200)')
@click.option('--repeat_length',
              type=int,
              default=160,
              help='Length of repeat sequence for marker removal (default: 160)')
@click.option('--genome_file',
              type=click.Path(exists=True, path_type=Path),
              default=Path("data/templates/BY4741_Toronto_2012.fsa"),
              help='Path to genome file (default: BY4741_Toronto_2012.fsa)')
@click.option('--ura3_file',
              type=click.Path(exists=True, path_type=Path),
              default=Path("data/component libraries/Saccharomyces cerevisiae/URA3.fasta"),
              help='Path to URA3 marker file. Default /data/component libraries/Saccharomyces cerevisiae/URA3.fasta')
def delete(upstream_homology_len, downstream_homology_len, repeat_length, genome_file, ura3_file):
    """Design scarless deletions in Saccharomyces cerevisiae.  
    \b\n
    Mechanism:\n
    \b
    5' ====up===[Target=====down]==repeat==3'   gDNA
           X                  X
    5' ====up==repeat==ura3==down==3'           Cassette
                ||
        **transformation**
                \\/
    5' =====up=repeat=Ura3====||
                X             ||                gDNA
         3' ===repeat====down=||
                ||
        **FOA counter selection**
                \\/
    5' =====up====repeat=== 3'                  gDNA
    """
    try:
        designer = DeletionDesigner(
            upstream_homology_len=upstream_homology_len,
            downstream_homology_len=downstream_homology_len,
            repeat_length=repeat_length,
            genome_file=genome_file,
            ura3_file=ura3_file
        )
        run_deletion_interactive_mode(designer)
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()
    
@cli.command()
@click.option('--upstream_homology_len', 
              type=int, 
              default=200,
              help='Length of upstream homology for recombination (default: 200)')
@click.option('--downstream_homology_len',
              type=int,
              default=200,
              help='Length of downstream homology for recombination (default: 200)')
@click.option('--repeat_length',
              type=int,
              default=160,
              help='Length of repeat sequence for marker removal (default: 160)')
@click.option('--genome_file',
              type=click.Path(exists=True, path_type=Path),
              default=Path("data/templates/BY4741_Toronto_2012.fsa"),
              help='Path to genome file (default: BY4741_Toronto_2012.fsa)')
@click.option('--ura3_file',
              type=click.Path(exists=True, path_type=Path),
              default=Path("data/component libraries/Saccharomyces cerevisiae/URA3.fasta"),
              help='Path to URA3 marker file. Default /data/component libraries/Saccharomyces cerevisiae/URA3.fasta')
def replace(upstream_homology_len, downstream_homology_len, repeat_length, genome_file, ura3_file):
    """Design pop-in/pop-out replacements in Saccharomyces cerevisiae
    \b\n
    Note URA can be positions up or downstream of ura3
    Target must be longer than 
    Mechanism:
    
    \b
    5'====up==[target=====================down]==repeat==3'     gDNA
          X                                 X
    5'====up==replacement==repeat==ura3===down===3'             cassette
                      ||
                 **transformation**
                      \\/
    5'===up==replacement==repeat==ura3==||
                            X           ||                      gDNA              
                      3'==repeat==down==||
                      ||
                **FOA counter selection**
                      \\/
    5'==up==replacement==repeate===3'                           gDNA   
    """
    try:
        designer = ReplaceDesigner(
            upstream_homology_len=upstream_homology_len,
            downstream_homology_len=downstream_homology_len,
            repeat_length=repeat_length,
            genome_file=genome_file,
            ura3_file=ura3_file
        )
        run_replace_interactive_mode(designer)
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()

@cli.command()
@click.option('--reuse_limit', 
              type=int, 
              default=5,
              help='Maximum times a PCR product can be reused (default: 5)')
def batch(reuse_limit):
    """Design batched assemblies for Saccharomyces cerevisiae transformations
    \b\n
    This command helps organize multiple DNA assemblies into efficient batches
    for parallel processing. 
    \b\n
    It generates:\b\n
    - Rationalized PCR instructions minimizing redundant reactions\b\n
    - Organized batches keeping construct PCRs together\b\n
    - Machine-readable instructions for liquid handling robots\b\n
    
    The input constructs must already exist as GenBank files in the output directory
    from previous tar/integrate/replace operations. For liquid handling instructions 
    templates should be aliquoted into the TemPlate plate, and the name added to the 
    coresponding well in data/templates/TemPlates.xlsx 
    
    The GenBank files must contain:
    - misc_feature annotations for components
    - primer_bind annotations for assembly primers
    """
    try:
        designer = BatchDesigner(reuse_limit=reuse_limit)
        run_batch_interactive_mode(designer)
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()
    
if __name__ == '__main__':
    cli()
    