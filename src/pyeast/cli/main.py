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

def print_construct_grid(available_constructs: dict) -> None:
    """Display available constructs for batch processing in a formatted Rich table."""
    table = Table(title="Available Constructs for Batch Processing")
    table.add_column("Name", style="cyan")
    table.add_column("Type", style="blue")
    table.add_column("Source", style="yellow")
    table.add_column("Topology", style="green")
    table.add_column("Length", justify="right", style="blue")
    table.add_column("Components", justify="right", style="magenta")
    table.add_column("Parts", style="white")

    for name, record in available_constructs.items():
        is_circular = record.annotations.get('topology', '').lower() == 'circular'
        source = record.annotations.get('source_folder', 'output/')
        source_type = record.annotations.get('source_type', 'unknown')
        components = [
            feature.qualifiers.get("label", ["Unlabeled"])[0]
            for feature in record.features
            if feature.type == "misc_feature"
        ]
        component_str = ", ".join(components)
        if len(component_str) > 120:
            component_str = component_str[:117] + "..."
        table.add_row(
            name,
            source_type,
            source,
            "Circular" if is_circular else "Linear",
            f"{len(record):,} bp",
            str(len(components)),
            component_str or "No components found",
        )
    console.print(table)


def get_batch_selections(available_constructs: dict) -> list[str] | None:
    """Interactively select constructs for batch assembly.

    Returns a list of construct names, or None if the user cancels.
    """
    session = PromptSession()
    completer = WordCompleter(list(available_constructs.keys()), ignore_case=True)

    while True:
        try:
            user_input = session.prompt(
                "\nEnter names of constructs to assemble (space-separated): ",
                completer=completer,
            ).strip()

            if not user_input:
                console.print("[yellow]No constructs selected[/yellow]")
                continue

            selected = user_input.split()
            invalid = [name for name in selected if name not in available_constructs]

            if invalid:
                console.print(f"[red]Invalid construct(s): {', '.join(invalid)}[/red]")
                continue

            console.print("\n[green]Selected constructs:[/green]")
            for name in selected:
                record = available_constructs[name]
                is_circular = record.annotations.get('topology', '').lower() == 'circular'
                console.print(f"  {name} ({'Circular' if is_circular else 'Linear'})")

            if click.confirm("\nProceed with these constructs?"):
                return selected

        except KeyboardInterrupt:
            if click.confirm("\nCancel construct selection?"):
                return None
            continue


def get_gg_assembly_order(sequences: dict) -> tuple[list[list[str]], bool] | None:
    """Interactively get GG assembly order with multiplex support.

    Returns (assemblies, is_multiplex) or None if the user cancels.
    assemblies is a list of assemblies; each assembly is a list of sequence names.
    """
    session = PromptSession()
    sequence_completer = WordCompleter(list(sequences.keys()), ignore_case=True)
    is_multiplex = False

    while True:
        try:
            console.print("\n[blue]Enter sequences to assemble (space-separated)[/blue]")
            console.print("[dim]Use TAB for autocompletion[/dim]")
            console.print("[dim]Remember that golden gate parts need to be assembled in order[/dim]")
            console.print("[blue]Separate components with / to multiplex, use /all[red]X[/red] to select all parts of type [red]x[/red][/blue]")

            user_input = session.prompt("Sequences: ", completer=sequence_completer)
            selected = user_input.split()

            if not selected:
                console.print("[yellow]No sequences selected[/yellow]")
                continue

            invalid_selections = []
            for selection in selected:
                if selection.startswith("/all"):
                    part_type = selection[4:].strip()
                    if not any(name.split("_")[0] == part_type for name in sequences):
                        invalid_selections.append(f"{selection} (no parts found with type '{part_type}')")
                elif "/" in selection:
                    for part_name in selection.split("/"):
                        if part_name not in sequences:
                            invalid_selections.append(f"{part_name} (from multiplex selection '{selection}')")
                else:
                    if selection not in sequences:
                        invalid_selections.append(selection)

            if invalid_selections:
                console.print(f"[red]Invalid selection(s): {', '.join(invalid_selections)}[/red]")
                continue

            assemblies = [[]]
            is_multiplex = False
            for selection in selected:
                if selection.startswith("/all"):
                    is_multiplex = True
                    part_type = selection[4:].strip()
                    multiplex_parts = [n for n in sequences if n.split("_")[0] == part_type]
                    console.print(f"{len(multiplex_parts)} type {part_type} parts for multiplexing")
                    new_assemblies = []
                    for assembly in assemblies:
                        for part in multiplex_parts:
                            new_assemblies.append(assembly + [part])
                    assemblies = new_assemblies
                elif "/" in selection:
                    is_multiplex = True
                    part_names = selection.split("/")
                    console.print(f"{len(part_names)} parts for multiplexing at this position")
                    new_assemblies = []
                    for assembly in assemblies:
                        for part in part_names:
                            new_assemblies.append(assembly + [part])
                    assemblies = new_assemblies
                else:
                    for assembly in assemblies:
                        assembly.append(selection)

            if is_multiplex:
                console.print(f"{len(assemblies)} constructs for assembly")
            else:
                console.print("One construct for assembly")

            return assemblies, is_multiplex

        except KeyboardInterrupt:
            if click.confirm("\nDo you want to exit?"):
                return None
            continue


def get_gg_liquid_handler(instruments: list) -> str:
    """Interactively select a liquid handler for GG instruction generation.

    Returns the selected instrument name as a lowercase string.
    """
    session = PromptSession()
    console.print("\n[blue]Available liquid handlers:[/blue]")
    for i, instrument in enumerate(instruments, 1):
        console.print(f"  {i}: {instrument}")
    completer = WordCompleter(instruments, ignore_case=True)
    selection = session.prompt(
        "Select a liquid handler for instruction formatting: ",
        completer=completer,
    ).strip().lower()
    return selection


def print_sequence_grid(sequences: dict, title: str = "Available Sequences"):
    """Display sequences in a formatted Rich table."""
    table = Table(title=title)
    table.add_column("Name", style="cyan")
    table.add_column("Length", justify="right", style="green")
    table.add_column("Description", style="white")

    for name, seq in sequences.items():
        table.add_row(
            name,
            f"{len(seq)} bp",
            seq.description[:150] + "..." if len(seq.description) > 149 else seq.description
        )

    console.print(table)


def print_integration_sites(sites: dict, title: str = "Available Integration Sites"):
    """Display integration sites in a formatted Rich table."""
    table = Table(title=title)
    table.add_column("Name", style="cyan")
    table.add_column("Up", justify="right", style="green")
    table.add_column("Down", justify="right", style="green")
    table.add_column("Description", style="white")

    for name, (up_seq, down_seq) in sites.items():
        table.add_row(
            name,
            f"{len(up_seq)} bp",
            f"{len(down_seq)} bp",
            up_seq.description[:150] + "..." if len(up_seq.description) > 149 else up_seq.description
        )

    console.print(table)


def display_instructions(instructions: list[list[str]]):
    """Display assembly instructions in a formatted Rich table."""
    table = Table(title="Assembly Instructions")
    table.add_column("Part Name", style="bold cyan")
    table.add_column("Fwd Primer", style="green")
    table.add_column("Fwd Plate", style="bold cyan")
    table.add_column("Well", style="bold blue")
    table.add_column("Rev Primer", style="green")
    table.add_column("Rev Plate", style="bold cyan")
    table.add_column("Well", style="bold blue")
    table.add_column("Template", style="magenta")
    table.add_column("Size", justify="right")

    for line in instructions:
        table.add_row(
            str(line[0]),
            str(line[2]),
            str(line[1] if line[1] != "N/A" else "[bold red]Not found[/bold red]"),
            str(line[3] if line[3] != "N/A" else "[bold red]Not found[/bold red]"),
            str(line[5]),
            str(line[4] if line[4] != "N/A" else "[bold red]Not found[/bold red]"),
            str(line[6] if line[6] != "N/A" else "[bold red]Not found[/bold red]"),
            str(line[7] if line[7] != "Not found" else "[bold red]Not found[/bold red]"),
            str(line[8])
        )

    console.print(table)


def get_tar_assembly_order(sequences: dict) -> list[str] | None:
    """Interactively get the TAR assembly order from the user with autocomplete."""
    session = PromptSession()
    sequence_completer = WordCompleter(sequences.keys(), ignore_case=True)

    while True:
        try:
            console.print("\n[blue]Enter sequences to assemble (space-separated)[/blue]")
            console.print("[dim]Use TAB for autocompletion[/dim]")

            user_input = session.prompt("Sequences: ", completer=sequence_completer)
            selected = user_input.split()

            if not selected:
                console.print("[yellow]No sequences selected[/yellow]")
                continue

            invalid = [name for name in selected if name not in sequences]
            if invalid:
                console.print(f"[red]Invalid sequence(s): {', '.join(invalid)}[/red]")
                continue

            console.print("\n[green]Selected sequences:[/green]")
            for i, name in enumerate(selected, 1):
                console.print(f"{i}. {name}")

            if click.confirm("\nProceed with these sequences?"):
                return selected

        except KeyboardInterrupt:
            if click.confirm("\nDo you want to exit?"):
                return None
            continue


def get_integration_selections(
    components: dict,
    int_sites: dict,
) -> tuple[list[str], str] | None:
    """Interactively get component order and integration site from the user."""
    session = PromptSession()
    component_completer = WordCompleter(list(components.keys()), ignore_case=True)
    int_site_completer = WordCompleter(list(int_sites.keys()), ignore_case=True)

    while True:
        try:
            print_sequence_grid(components, "Available Components")
            print_integration_sites(int_sites)

            user_input = session.prompt(
                "\nEnter the names of the components you want to assemble, in order (space-separated): ",
                completer=component_completer,
            )

            component_names = []
            valid = True
            for name in user_input.split():
                matches = [k for k in components if k.lower() == name.lower()]
                if matches:
                    component_names.append(matches[0])
                else:
                    console.print(f"[red]Invalid component name: {name}[/red]")
                    valid = False
                    break

            if not valid or not component_names:
                if not component_names:
                    console.print("[yellow]No components selected[/yellow]")
                continue

            while True:
                int_site = session.prompt(
                    "\nEnter the name of the integration site: ",
                    completer=int_site_completer,
                ).strip()
                matches = [s for s in int_sites if s.lower() == int_site.lower()]
                if matches:
                    int_site_name = matches[0]
                    break
                console.print(f"[red]Invalid integration site name: {int_site}[/red]")

            console.print("\n[green]Selected assembly:[/green]")
            console.print("Components:")
            for i, name in enumerate(component_names, 1):
                console.print(f"{i}. {name}")
            console.print(f"Integration site: {int_site_name}")

            if click.confirm("Is this correct?"):
                return component_names, int_site_name

            console.print("Let's try again.")

        except ValueError as e:
            console.print(f"[red]Error: {e}[/red]")
        except KeyboardInterrupt:
            if click.confirm("\nCancel selection?"):
                raise click.Abort()


def run_tar_interactive_mode(designer: TARDesigner):
    """Run TAR design in interactive mode"""
    try:
        components_dir = get_component_dir()

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Loading sequences...", total=None)
            sequences = designer.load_and_get_sequences(components_dir)
            progress.update(task_id, completed=True)

        if not sequences:
            console.print("[red]No sequences found in directory[/red]")
            return

        print_sequence_grid(sequences)
        assembly_order = get_tar_assembly_order(sequences)

        if not assembly_order:
            return

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing...", total=None)
            result = designer.design(
                library_path=components_dir,
                assembly_order=assembly_order,
                primer_folder=get_primers_path(),
                template_folder=get_templates_path(),
                name="assembly",
            )
            progress.update(task_id, completed=True)

        display_instructions(result.instructions)

        if not click.confirm("\nProceed with assembly?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return

        output_prefix = get_output_prefix()
        result.assembly.name = output_prefix.name
        result.name = output_prefix.name

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)

            try:
                result.save(output_prefix)
                img_data, fig = visualise_genbank(f"{output_prefix}.gb")
                save_figure(fig, f"{output_prefix}_map.png")
                progress.update(task_id, completed=True)

                console.print("\n[bold green]Files saved:[/bold green]")
                console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
                console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
                console.print(f"[green]Assembly instructions: {output_prefix}_instructions.tsv[/green]")
                console.print(f"[green]All primers: {output_prefix}_all_primers.tsv[/green]")
                if result.missing_primers:
                    console.print(f"[green]Missing primers: {output_prefix}_missing_primers.tsv[/green]")

                img = Image.open(io.BytesIO(img_data.getvalue()))
                img.show()

            except Exception as e:
                console.print(f"[red]Error saving files: {str(e)}[/red]")
                raise

        console.print("\n[bold green]checkmark[/bold green] Plasmid design complete!")

    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise

def run_integration_interactive_mode(designer: IntegrationDesigner):
    """Run integration design in interactive mode"""
    try:
        components_dir = get_component_dir()

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Loading sequences...", total=None)
            designer.load_sequences(components_dir)
            progress.update(task_id, completed=True)

        if not designer.components or not designer.int_sites:
            console.print("[red]No sequences found[/red]")
            return

        selections = get_integration_selections(designer.components, designer.int_sites)
        if not selections:
            return
        assembly_order, integration_site_name = selections

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing...", total=None)
            result = designer.design(
                components_dir=components_dir,
                assembly_order=assembly_order,
                integration_site_name=integration_site_name,
                primer_folder=get_primers_path(),
                template_folder=get_templates_path(),
                name="assembly",
            )
            progress.update(task_id, completed=True)

        display_instructions(result.instructions)

        if not click.confirm("\nProceed with assembly?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return

        output_prefix = get_output_prefix()
        result.assembly.name = output_prefix.name
        result.name = output_prefix.name

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)

            try:
                result.save(output_prefix)
                img_data, fig = visualise_genbank(f"{output_prefix}.gb")
                save_figure(fig, f"{output_prefix}_map.png")
                progress.update(task_id, completed=True)

                console.print("\n[bold green]Files saved:[/bold green]")
                console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
                console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
                console.print(f"[green]Assembly instructions: {output_prefix}_instructions.tsv[/green]")
                console.print(f"[green]All primers: {output_prefix}_all_primers.tsv[/green]")
                if result.missing_primers:
                    console.print(f"[green]Missing primers: {output_prefix}_missing_primers.tsv[/green]")

                img = Image.open(io.BytesIO(img_data.getvalue()))
                img.show()

            except Exception as e:
                console.print(f"[red]Error saving files: {str(e)}[/red]")
                raise

        console.print("\n[bold green]checkmark[/bold green] Integration design complete!")

    except KeyboardInterrupt:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()

def get_replacement_selection(components_dir: Path):
    """Interactively select a replacement sequence from a component library.

    Args:
        components_dir: Directory containing sequence FASTA files.

    Returns:
        Selected SeqRecord, or None if the user cancels.
    """
    from pyeast.utils.sequence_utils import load_sequences as _load_sequences
    sequences = _load_sequences(components_dir)
    if not sequences:
        console.print("[red]No sequences found in directory[/red]")
        return None

    print_sequence_grid(sequences, "Available Replacement Sequences")
    session = PromptSession()
    sequence_completer = WordCompleter(list(sequences.keys()), ignore_case=True)

    while True:
        try:
            console.print("\n[blue]Select a sequence to insert:[/blue]")
            user_input = session.prompt(
                "Sequence name: ",
                completer=sequence_completer,
            ).strip()

            matches = [seq for name, seq in sequences.items() if name.lower() == user_input.lower()]
            if matches:
                return matches[0]
            console.print("[red]Invalid sequence name[/red]")

        except KeyboardInterrupt:
            return None


def run_deletion_interactive_mode(designer: DeletionDesigner):
    """Run deletion design in interactive mode."""
    try:
        session = PromptSession()

        # Get target sequence from user
        while True:
            console.print("\n[blue]Enter the DNA sequence you want to delete:[/blue]")
            console.print("[dim]Use only A, T, G, and C[/dim]")
            sequence = session.prompt("Sequence: ").upper().strip()

            if not set(sequence).issubset({'A', 'T', 'G', 'C'}):
                console.print("[red]Invalid DNA sequence. Please use only A, T, G, and C.[/red]")
                continue
            if len(sequence) < designer.downstream_homology_len:
                console.print("[red]Target sequence is shorter than the homology lengths. Adjust with --downstream_homology_len[/red]")
                continue
            break

        # Find target in genome (pre-run for display/confirm before committing)
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Searching genome for target sequence...", total=None)
            location = designer.find_target_sequence(designer.genome_file, sequence)
            progress.update(task_id, completed=True)

        if not location:
            console.print("[red]Target sequence not found in genome.[/red]")
            return

        chrom_id, start, end, orientation = location
        console.print("\n[green]Target sequence found:[/green]")
        console.print(f"Chromosome: {chrom_id}")
        console.print(f"Position: {start}-{end}")
        console.print(f"Orientation: {orientation}")

        if not click.confirm("\nProceed with deletion design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return

        # Run full design (pass genome_location to skip second genome search)
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing deletion cassette...", total=None)
            result = designer.design(
                target_sequence=sequence,
                name="deletion_cassette",
                genome_location=location,
            )
            progress.update(task_id, completed=True)

        # Show screening strategy
        table = Table()
        table.add_column("Stage", style="cyan")
        table.add_column("Product Size", style="green")
        table.add_row("Parent strain", str(result.product_sizes['parent']))
        table.add_row("After transformation", str(result.product_sizes['deletion_intermediate']))
        table.add_row("Final deletion", str(result.product_sizes['final_deletion']))
        console.print("\n[bold cyan]Screening Strategy:[/bold cyan]")
        console.print(table)
        console.print("\n[bold cyan]Screening Primers:[/bold cyan]")
        console.print(f"Forward: {result.forward_primer}")
        console.print(f"Reverse: {result.reverse_primer}")

        if not click.confirm("\nSave deletion design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return

        output_prefix = get_output_prefix()
        result.cassette.name = output_prefix.name
        result.name = output_prefix.name

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)
            result.save(output_prefix)
            img_data, fig = visualise_genbank(f"{output_prefix}.gb")
            save_figure(fig, f"{output_prefix}_map.png")
            img = Image.open(io.BytesIO(img_data.getvalue()))
            img.show()
            progress.update(task_id, completed=True)

            console.print("\n[bold green]Files saved:[/bold green]")
            console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
            console.print(f"[green]FASTA file: {output_prefix}.fasta[/green]")
            console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
            console.print(f"[green]Screening primers: {output_prefix}_screening_primers.tsv[/green]")

        console.print("\n[bold green]checkmark[/bold green] Deletion design complete!")

    except click.Abort:
        console.print("\n[yellow]Operation cancelled[/yellow]")
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {str(e)}")
        raise click.Abort()

def run_replace_interactive_mode(designer: ReplaceDesigner):
    """Run replace design in interactive mode."""
    try:
        session = PromptSession()

        # Get target sequence from user
        while True:
            console.print("\n[blue]Enter the DNA sequence you want to replace:[/blue]")
            console.print("[dim]Use only A, T, G, and C[/dim]")
            sequence = session.prompt("Sequence: ").upper().strip()

            if not set(sequence).issubset({'A', 'T', 'G', 'C'}):
                console.print("[red]Invalid DNA sequence. Please use only A, T, G, and C.[/red]")
                continue

            min_target_len = min(designer.downstream_homology_len, designer.upstream_homology_len)
            if len(sequence) < min_target_len:
                console.print("[red]Target sequence is shorter than the homology lengths. Adjust with --downstream_homology_len and --upstream_homology_len[/red]")
                continue
            break

        # Find target in genome (pre-run for display/confirm before committing)
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Searching genome for target sequence...", total=None)
            location = designer.find_target_sequence(designer.genome_file, sequence)
            progress.update(task_id, completed=True)

        if not location:
            console.print("[red]Target sequence not found in genome.[/red]")
            return

        chrom_id, start, end, orientation = location
        console.print("\n[green]Target sequence found:[/green]")
        console.print(f"Chromosome: {chrom_id}")
        console.print(f"Position: {start}-{end}")
        console.print(f"Orientation: {orientation}")

        # Get replacement sequence
        components_dir = get_component_dir()
        replacement_sequence = get_replacement_selection(components_dir)
        if not replacement_sequence:
            console.print("[yellow]No replacement sequence selected. Design cancelled.[/yellow]")
            return

        # Get marker position
        while True:
            console.print("\n[blue]Do you want the URA3 marker upstream or downstream of the replacement?[/blue]")
            position = session.prompt(
                "Position (up/down): ",
                completer=WordCompleter(['up', 'down']),
            ).lower().strip()
            if position in ('up', 'down'):
                marker_position = 'upstream' if position == 'up' else 'downstream'
                break
            console.print("[red]Invalid position. Please enter 'up' or 'down'.[/red]")

        # Run full design (pass genome_location to skip second genome search)
        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Designing replacement cassette...", total=None)
            result = designer.design(
                target_sequence=sequence,
                replacement_sequence=replacement_sequence,
                marker_position=marker_position,
                name="replacement_cassette",
                genome_location=location,
            )
            progress.update(task_id, completed=True)

        # Show screening strategy
        table = Table()
        table.add_column("Stage", style="cyan")
        table.add_column("Product Size", style="green")
        table.add_row("Parent strain", str(result.product_sizes['parent']))
        table.add_row("After transformation", str(result.product_sizes['replacement_intermediate']))
        table.add_row("Final replacement", str(result.product_sizes['final_replacement']))
        console.print("\n[bold cyan]Screening Strategy:[/bold cyan]")
        console.print(table)
        console.print("\n[bold cyan]Screening Primers:[/bold cyan]")
        console.print(f"Forward: {result.forward_primer}")
        console.print(f"Reverse: {result.reverse_primer}")

        if not click.confirm("\nSave replacement design?"):
            console.print("[yellow]Design cancelled[/yellow]")
            return

        output_prefix = get_output_prefix()
        result.cassette.name = output_prefix.name
        result.name = output_prefix.name

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
            task_id = progress.add_task("Saving outputs...", total=None)
            result.save(output_prefix)
            img_data, fig = visualise_genbank(f"{output_prefix}.gb")
            save_figure(fig, f"{output_prefix}_map.png")
            img = Image.open(io.BytesIO(img_data.getvalue()))
            img.show()
            progress.update(task_id, completed=True)

            console.print("\n[bold green]Files saved:[/bold green]")
            console.print(f"[green]GenBank file: {output_prefix}.gb[/green]")
            console.print(f"[green]FASTA file: {output_prefix}.fasta[/green]")
            console.print(f"[green]Sequence map: {output_prefix}_map.png[/green]")
            console.print(f"[green]Screening primers: {output_prefix}_screening_primers.tsv[/green]")

        console.print("\n[bold green]checkmark[/bold green] Replacement design complete!")

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
        print_construct_grid(designer.available_constructs)
        selected = get_batch_selections(designer.available_constructs)
        if selected is None:
            console.print("\n[yellow]Operation cancelled[/yellow]")
            return
        designer.selected_constructs = {
            name: designer.available_constructs[name] for name in selected
        }

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
                print_sequence_grid(sequences)
                result = get_gg_assembly_order(sequences)

                if result is None:
                    console.print("[yellow]No assembly selected[/yellow]")
                    return
                assembly_order, is_multiplex = result
                designer.assemblies_names = assembly_order
                designer.multiplex = is_multiplex

                # Get an output prefix
                output_path = get_output_prefix()
                prefix = output_path.name

                # Map parts to plasmids and run assembly simulation
                with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}")) as progress:
                    task_id = progress.add_task("Assembling selected sequences...", total=None)
                    designer.get_plasmid_names()
                    assembly_sim = designer.gg_assembly(prefix)
                    progress.update(task_id, completed=True)

                if assembly_sim.errors:
                    if click.confirm("Return to part selection?"):
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

                console.print("[green]Assembly Successful![/green]")
                if click.confirm("Save Outputs?"):
                    liquid_handler = get_gg_liquid_handler(designer.instruments)
                    designer.gg_save_output(str(output_path.parent))
                    designer.gg_instructions(str(output_path.parent), prefix, liquid_handler=liquid_handler)

                    if len(assembly_order) == 1:
                        assemblies_file = output_path.parent / "input.txt"
                    else:
                        assemblies_file = output_path.parent / "inputs.txt"

                    with open(assemblies_file, 'w') as f:
                        f.write(tabulate(assembly_order))

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
        console.print(f"[red]Clone failed:[/red]\n{result.stderr}")
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
            console.print(f"[red]Directory not found: {target}[/red]")
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
        console.print(f"[green]Output directory:[/green] {config.output_dir}")

        if not click.confirm("\nWould you like to reconfigure?", default=False):
            return

        console.print("\n[cyan]What would you like to reconfigure?[/cyan]")
        console.print("  1. Data directory (enter existing path)")
        console.print("  2. Data directory (clone PYEAST data repository)")
        console.print("  3. Output directory")
        choice = click.prompt("Choice", type=click.Choice(["1", "2", "3"]))

        if choice == "1":
            new_path = click.prompt("Data directory path", type=click.Path())
            target = Path(new_path)
            if not target.exists():
                console.print(f"[red]Directory not found: {target}[/red]")
                raise click.Abort()
            _write_config(target.resolve(), None)
            console.print(f"[green]Configured PYEAST to use data at {target}[/green]")
            console.print(f"[dim]Config saved to: {Path.home() / '.pyeast' / 'config.yaml'}[/dim]")
            return

        if choice == "3":
            new_output = click.prompt("Output directory path", type=click.Path())
            output_path = Path(new_output)
            _write_config(config.data_dir, output_path.resolve())
            console.print(f"[green]Output directory set to {output_path}[/green]")
            console.print(f"[dim]Config saved to: {Path.home() / '.pyeast' / 'config.yaml'}[/dim]")
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
