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

import io
from datetime import datetime
from pathlib import Path

import click
from Bio import SeqIO
from PIL import Image
from prompt_toolkit import PromptSession
from prompt_toolkit.completion import WordCompleter
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table
from tabulate import tabulate

from pyeast.core.batch import BatchDesigner
from pyeast.core.deletion import DeletionDesigner
from pyeast.core.gg import ggDesigner
from pyeast.core.integration import IntegrationDesigner
from pyeast.core.replace import ReplaceDesigner
from pyeast.core.tar import TARDesigner
from pyeast.utils.path_utils import ensure_output_dir_exists, get_component_libraries_path, get_primers_path, get_templates_path
from pyeast.utils.visualisation import save_figure, visualise_genbank

console = Console()

DATA_REPO_URL = "https://github.com/TomLoan/PYEAST_data.git"
_DATA_REPO_CLONE_DIR = Path.home() / ".pyeast" / "data-repo"


def get_output_prefix() -> Path:
    """Get output prefix from user - creates a subfolder and returns path to files within it.

    This function:
    1. Asks user for a name
    2. Checks if a subfolder already exists with that name
    3. If it exists, warns user and asks for confirmation to overwrite
    4. Creates a subfolder in output/ with that name
    5. Returns a Path object: output/name/name

    Example:
        If user inputs "my_construct", this returns:
        Path("output/my_construct/my_construct")

        Then files are saved as:
        - output/my_construct/my_construct.gb
        - output/my_construct/my_construct_primers.tsv
        etc.

    Returns:
        Path object pointing to files within the created subfolder.
        Use .name property to get just the filename without path.
    """
    session = PromptSession()
    output_dir = ensure_output_dir_exists()

    # Output directory is ensured by ensure_output_dir_exists()

    while True:
        try:
            console.print("\n[blue]Enter a name for your output files[/blue]")
            console.print("[dim]A subfolder will be created: output/your_name/[/dim]")
            user_input = session.prompt("Name: ").strip()

            # Basic validation
            if not user_input:
                console.print("[red]Please enter a name[/red]")
                continue

            # Replace spaces with underscores first
            user_input = user_input.replace(" ", "_")

            # Remove any path separators and other invalid filename characters
            # Keep alphanumeric, hyphens, underscores
            safe_name = "".join(c for c in user_input if c not in r'\/.:*?"<>|')

            if not safe_name:
                console.print("[red]Invalid name after removing special characters[/red]")
                continue

            # Check if subfolder already exists
            output_subfolder = output_dir / safe_name

            if output_subfolder.exists():
                # List files in the existing directory
                existing_files = list(output_subfolder.glob("*"))

                console.print(f"\n[yellow]Warning: Directory '{safe_name}' already exists[/yellow]")
                if existing_files:
                    console.print(f"[yellow]It contains {len(existing_files)} file(s):[/yellow]")
                    # Show first few files
                    for file in existing_files[:5]:
                        console.print(f"[dim]  - {file.name}[/dim]")
                    if len(existing_files) > 5:
                        console.print(f"[dim]  ... and {len(existing_files) - 5} more[/dim]")

                # Ask user if they want to overwrite
                if not click.confirm("\nOverwrite existing experiment?"):
                    console.print("[yellow]Please choose a different name[/yellow]")
                    continue

                console.print("[yellow]Proceeding with overwrite...[/yellow]")

            # Create subfolder (or use existing one)
            output_subfolder.mkdir(exist_ok=True)

            console.print(f"[dim]Files will be saved to: {output_subfolder}/[/dim]")

            # Return path pointing to files within the subfolder
            # e.g., Path("output/my_construct/my_construct")
            return output_subfolder / safe_name

        except KeyboardInterrupt:
            raise click.Abort()



def handle_machine_instructions(designer: BatchDesigner, output_prefix: str) -> None:
    """Ask user if they want to generate machine instructions and handle the response.

    Args:
        designer: BatchDesigner instance with completed assembly instructions
        output_prefix: Path prefix for output files (same as human instructions)
    """
    session = PromptSession()

    if click.confirm("\nWould you like to generate machine instructions for liquid handling?"):
        # Select from available machines
        machines = ['epMotion', 'Janus', 'Hamilton']

        print("\nAvailable liquid handling machines:")
        for i, machine in enumerate(machines, 1):
            print(f"{i}. {machine}")

        completer = WordCompleter(machines, ignore_case=True)

        selection = session.prompt("Which machine would you like to use? \n", completer = completer).strip().lower()

        if selection == "epmotion":
            if click.confirm("Generate epMotion instructions?"):
                try:
                    timestamp = datetime.now().strftime("%H-%M-%d-%b-%Y").upper()
                    # Write instructions for the PCR set up
                    designer.generate_epmotion_instructions(output_prefix, timestamp)
                    # Write instructions for the assembly transformations
                    designer.generate_machine_assembly_instructions(output_prefix, 'epmotion', timestamp)
                except Exception as e:
                    console.print(f"[red]Error generating machine instructions: {str(e)}[/red]")

        if selection == "janus"or selection =='hamilton':
            if click.confirm(f"Generate worklist for {selection}?"):
                try:
                    timestamp = datetime.now().strftime("%H-%M-%d-%b-%Y").upper()
                    # Write instructions for the PCR set up
                    designer.generate_janus_instructions(output_prefix, timestamp)
                    # Write instructions for the assembly transformations
                    designer.generate_machine_assembly_instructions(output_prefix, 'janus', timestamp)
                except Exception as e:
                    console.print(f"[red]Error generating machine instructions: {str(e)}[/red]")

    else:
        pass
# Could edit this function so the user can select more than one library if desired?
# Not sure What I'd do about name conflicts in this situation though, more trouble than it's worth?
def get_component_dir() -> Path:
    """Get component directory with simple autocompletion, including private libraries"""
    # Check both public and private base directories
    base_dir = get_component_libraries_path()
    private_base_dir = get_component_libraries_path(private=True)

    if not base_dir.exists():
        console.print("[red]Error: Default components directory not found[/red]")
        raise click.Abort()

    # Collect subdirectories from both locations
    subdirs = set()

    # Add public directories
    for d in base_dir.iterdir():
        if d.is_dir():
            subdirs.add(d.name)

    # Add private directories (if private base exists)
    if private_base_dir.exists():
        for d in private_base_dir.iterdir():
            if d.is_dir():
                subdirs.add(d.name)

    # Convert to sorted list for display
    subdirs = sorted(list(subdirs))

    if not subdirs:
        console.print("[red]No component libraries found[/red]")
        raise click.Abort()

    # Create completer
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

            # Always return the public path - load_sequences will handle both
            selected_dir = base_dir / user_input

            # Verify at least one location exists
            private_dir = private_base_dir / user_input
            if selected_dir.exists() or private_dir.exists():
                return selected_dir

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
            designer.check_primer_locations(get_primers_path())
            progress.update(task_id, completed=True)

            # Find templates
            task_id = progress.add_task("Finding templates...", total=None)
            designer.find_templates(get_templates_path())
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
            assembly.name = output_prefix.name
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
            designer.check_primer_locations(get_primers_path())
            progress.update(task_id, completed=True)

            # Find templates
            task_id = progress.add_task("Finding templates...", total=None)
            designer.find_templates(get_templates_path())
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
        assembly.name = output_prefix.name
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
        console.print("\n[green]Target sequence found:[/green]")
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
        designer.deletion_cassette.name = output_prefix.name

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
                f.write(f"{output_prefix.name}_ScreenF\t{forward_primer}\n")
                f.write(f"{output_prefix.name}_ScreenR\t{reverse_primer}")

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
        console.print("\n[green]Target sequence found:[/green]")
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

        designer.replacement_cassette.name = output_prefix.name

        # Save results
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)

            # Save cassette as GenBank and FASTA
            SeqIO.write(designer.replacement_cassette, f"{output_prefix}.gb", "genbank")
            SeqIO.write(designer.replacement_cassette, f"{output_prefix}.fasta", "fasta")

            # Save primers
            forward_primer, reverse_primer = designer.screening_primers
            with open(f"{output_prefix}_screening_primers.tsv", 'w') as f:
                f.write(f"{output_prefix.name}_ScreenF\t{forward_primer}\n")
                f.write(f"{output_prefix.name}_ScreenR\t{reverse_primer}")

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
                # Save input record
                designer.save_input_record(output_prefix)
            except ValueError as e:
                progress.update(task_id, completed=True)
                console.print(f"\n[red]{str(e)}[/red]")
                return
            progress.update(task_id, completed=True)

        # Optionally generate machine instructions
        progress.update(task_id, completed=True)

        handle_machine_instructions(designer, output_prefix)

        progress.update(task_id, completed=True)

        console.print("\n[bold green]✓[/bold green] Batch design complete!")

    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()

def run_gg_interactive_mode(designer: ggDesigner):
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
        # Main workflow loop - returns suer to selection on failure
        while True:

            try:
                # Display sequences and get assembly order
                designer.print_sequence_grid(sequences)
                assembly_order = designer.get_assembly_order(sequences)
                #console.print(assembly_order)

                if not assembly_order:
                    console.print("[yellow]No assembly selected[/yellow]")
                    return

                # Get an output prefix
                output_path = get_output_prefix()
                #console.print(output_path)
                prefix = output_path.name

                # Convert to SeqRecord objects and map to plasmids
                with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
                    task_id = progress.add_task("Assembling selected sequences...", total = None)
                    designer.get_plasmid_names()
                    assembly_sim = designer.gg_assembly(prefix)
                    progress.update(task_id, completed = True)
                if assembly_sim.errors:
                    if click.confirm("Return to part selection?"):
                        #reset state and try again
                        designer.assemblies_names = None
                        designer.assembly_sequences = []
                        designer.plasmid_names = []
                        designer.part_to_plasmid_mapping = {}
                        designer.assembly = None
                        designer.assembly_sim = None
                        continue
                    else:
                        console.print("[yellow]Golden Gate design cancelled[/yellow]")
                        return
                # For Successful assembly ask about saving outputs
                console.print("[green]Assembly Successful![/green]")
                if click.confirm("Save Outputs?"):
                    console.print(output_path, prefix)
                    designer.gg_save_output(str(output_path.parent))
                    designer.gg_instructions(str(output_path.parent), prefix)

                    #Save the input to file for future reference
                    if len(assembly_order) == 1:
                        assemblies_file = output_path.parent / "input.txt"
                    elif len(assembly_order) > 1:
                        assemblies_file = output_path.parent / "inputs.txt"

                    with open(assemblies_file, 'w') as f:
                        f.write(tabulate(assembly_order))
                #exit loop on successful assembly
                break

            except KeyboardInterrupt:
                console.print("[yellow]Operation canceled[/yellow]")
                return
            except Exception as e:
                console.print(f"[red]Unexpected error: {str(e)}[/red]")
                if click.confirm("Return to part selection and try again?"):
                    continue
                else:
                    raise

    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")

def _write_config(data_dir: Path, output_dir) -> None:
    """Write ~/.pyeast/config.yaml with the given data and optional output paths."""
    import yaml
    config_file = Path.home() / ".pyeast" / "config.yaml"
    config_file.parent.mkdir(parents=True, exist_ok=True)
    config_data = {'data_dir': str(data_dir)}
    if output_dir:
        config_data['output_dir'] = str(output_dir)
    with open(config_file, 'w') as f:
        yaml.dump(config_data, f)


def _clone_data_repo(subprocess) -> None:
    """Clone the PYEAST data repository to ~/.pyeast/data-repo/ and write config."""
    clone_dir = _DATA_REPO_CLONE_DIR

    if clone_dir.exists():
        console.print(f"[green]Data repository already cloned at {clone_dir}[/green]")
        console.print(f"[dim]To update: git -C {clone_dir} pull[/dim]")
        _write_config(clone_dir, None)
        console.print(f"[green]Config updated to use data at {clone_dir}[/green]")
        return

    console.print(f"Cloning PYEAST data repository from {DATA_REPO_URL} ...")
    result = subprocess.run(
        ["git", "clone", DATA_REPO_URL, str(clone_dir)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        console.print(f"[red]✗ Clone failed:[/red]\n{result.stderr}")
        raise click.Abort()

    _write_config(clone_dir, None)
    console.print(f"[green]Data repository cloned to {clone_dir}[/green]")
    console.print(f"[green]Configured PYEAST to use data at {clone_dir}[/green]")
    console.print(f"[dim]Config saved to: {Path.home() / '.pyeast' / 'config.yaml'}[/dim]")


@click.group()
def cli():
    """PYeast: Python tools for yeast genetic engineering

    Created by Tom Loan"""
    pass

@cli.command()
@click.option('--data-dir', type=click.Path(), help='Path to existing data directory to use')
@click.option('--output-dir', type=click.Path(), help='Path to output directory for results (optional)')
def init(data_dir, output_dir):
    """Initialize PYEAST data directory configuration.

    Without --data-dir: clones the PYEAST data repository to ~/.pyeast/data-repo/
    automatically, or shows current configuration if already set up.

    With --data-dir: registers an existing data directory.

    Examples:
        pyeast init
        pyeast init --data-dir /path/to/PYEAST_data
        pyeast init --data-dir /path/to/data --output-dir /path/to/output
    """
    import subprocess

    from pyeast.config import get_config

    config = get_config()

    if data_dir:
        # User provided an existing data location
        target = Path(data_dir)
        if not target.exists():
            console.print(f"[red]✗ Directory not found: {target}[/red]")
            raise click.Abort()

        if output_dir:
            output_target = Path(output_dir)
            if not output_target.exists():
                console.print(f"[yellow]Output directory does not exist: {output_target}[/yellow]")
                console.print("[dim]It will be created when needed.[/dim]")

        _write_config(target.resolve(), Path(output_dir).resolve() if output_dir else None)
        console.print(f"[green]Configured PYEAST to use data at {target}[/green]")
        if output_dir:
            console.print(f"[green]Output directory set to {output_dir}[/green]")
        console.print(f"[dim]Config saved to: {Path.home() / '.pyeast' / 'config.yaml'}[/dim]")
        return

    # No --data-dir provided
    if config.data_dir.exists():
        # Already configured — show current state and offer to reconfigure
        console.print(f"[green]Data directory:[/green] {config.data_dir}")
        console.print(f"[dim]Output directory: {config.output_dir}[/dim]")

        if not click.confirm("\nWould you like to reconfigure?", default=False):
            return

        console.print("\n[cyan]How would you like to configure the data directory?[/cyan]")
        console.print("  1. Enter a path to an existing data directory")
        console.print("  2. Clone the public PYEAST data repository to ~/.pyeast/data-repo/")
        choice = click.prompt("Choice", type=click.Choice(["1", "2"]))

        if choice == "1":
            new_path = click.prompt("Data directory path", type=click.Path())
            target = Path(new_path)
            if not target.exists():
                console.print(f"[red]✗ Directory not found: {target}[/red]")
                raise click.Abort()
            _write_config(target.resolve(), None)
            console.print(f"[green]Configured PYEAST to use data at {target}[/green]")
            return
        # choice == "2" falls through to clone logic below

    # Clone the public data repo
    _clone_data_repo(subprocess)

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
              default=100,
              help='Length of downstream homology for recombination (default: 100)')
@click.option('--repeat_length',
              type=int,
              default=80,
              help='Length of repeat sequence for marker removal (default: 80)')
@click.option('--genome_file',
              type=click.Path(exists=True, path_type=Path),
              default=None,
              help='Path to genome file (default: BY4741_Toronto_2012.fsa from data directory)')
@click.option('--ura3_file',
              type=click.Path(exists=True, path_type=Path),
              default=None,
              help='Path to URA3 marker file (default: URA3.fasta from data/component_libraries/Saccharomyces_cerevisiae/)')
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
        # Set defaults if not provided
        if genome_file is None:
            default_genome = get_templates_path() / "BY4741_Toronto_2012.fsa"
            if not default_genome.exists():
                console.print(f"[red]Default genome file not found: {default_genome}[/red]")
                console.print("[yellow]Hint: Run 'pyeast init' to configure data directory[/yellow]")
                raise click.Abort()
            genome_file = default_genome
        if ura3_file is None:
            default_ura3 = get_component_libraries_path() / "Saccharomyces_cerevisiae" / "URA3.fasta"
            if not default_ura3.exists():
                console.print(f"[red]Default URA3 file not found: {default_ura3}[/red]")
                console.print("[yellow]Hint: Run 'pyeast init' to configure data directory[/yellow]")
                raise click.Abort()
            ura3_file = default_ura3

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
              default=80,
              help='Length of repeat sequence for marker removal (default: 80)')
@click.option('--genome_file',
              type=click.Path(exists=True, path_type=Path),
              default=None,
              help='Path to genome file (default: BY4741_Toronto_2012.fsa from data directory)')
@click.option('--ura3_file',
              type=click.Path(exists=True, path_type=Path),
              default=None,
              help='Path to URA3 marker file (default: URA3.fasta from data/component_libraries/Saccharomyces_cerevisiae/)')
def replace(upstream_homology_len, downstream_homology_len, repeat_length, genome_file, ura3_file):
    """Design pop-in/pop-out replacements in Saccharomyces cerevisiae
    \b\n
    Note URA can be positions up or downstream of ura3
    Target must be longer than both homology arms
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
        # Set defaults if not provided
        if genome_file is None:
            default_genome = get_templates_path() / "BY4741_Toronto_2012.fsa"
            if not default_genome.exists():
                console.print(f"[red]Default genome file not found: {default_genome}[/red]")
                console.print("[yellow]Hint: Run 'pyeast init' to configure data directory[/yellow]")
                raise click.Abort()
            genome_file = default_genome
        if ura3_file is None:
            default_ura3 = get_component_libraries_path() / "Saccharomyces_cerevisiae" / "URA3.fasta"
            if not default_ura3.exists():
                console.print(f"[red]Default URA3 file not found: {default_ura3}[/red]")
                console.print("[yellow]Hint: Run 'pyeast init' to configure data directory[/yellow]")
                raise click.Abort()
            ura3_file = default_ura3

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

@cli.command()
@click.option('--library',
              type = bool,
              default = False,
              help = 'Assemble selections in a single reaction to create a library of constructs (default: False)')
def gg(library):
    """===Experimental===
    Design golden gate cloning experiments in Saccharomyces cerevisiae
    \b\n
    Still being developed, use with caution. supports multiplex and library type assemblies
    Use / to separate component names you want to multiplex , or input /allX to select all components of type X
    The designer can handle parts input out of order, although this can make it hard to see what you're doing, and will
    idenify the correct enzyme for the parts you've selected automatically.
    """
    try:
        designer = ggDesigner(is_library = library)
        run_gg_interactive_mode(designer)
    except click.Abort:
        raise
    except Exception as e:
        console.print(f"[red]Error: {str(e)}[/red]")
        raise click.Abort()

if __name__ == '__main__':
    cli()
