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

# src/pyeast/utils/visualisation.py
"""
DNA visualisation tools for PYEAST. 
"""

# ===========================================================================




import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, CircularGraphicRecord, GraphicRecord
from Bio import SeqIO
from io import BytesIO
import itertools

def visualise_genbank(gb_file_path):
    """
    Visualise a GenBank file as a plasmid map with improved label consistency and readability.

    Args:
        gb_file_path (str): Path to the GenBank file.

    Returns:
        tuple: (BytesIO object containing the image data, matplotlib figure object)
    """
    # Parse the GenBank file
    record = SeqIO.read(gb_file_path, "genbank")

    # Create a list to store the features
    features = []

    # Define your organization's color palette
    color_palette = [
        "#1E22AA", "#004B87", "#6D2077", "#007377", "#007A53", "#DF1995",
        "#E87722", "#FFB81C", "#9FAEE5", "#71CC98", "#78BE20", "#2DCCD3"
    ]
    color_cycle = itertools.cycle(color_palette)

    for feature in record.features:
        if feature.type not in ["source", "primer_bind"]:
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            
            label = (feature.qualifiers.get("label", []) or 
                     feature.qualifiers.get("gene", []) or 
                     feature.qualifiers.get("product", []) or 
                     [f"{feature.type}_{start}-{end}"])[0]
            
            color = next(color_cycle)
            
            features.append(GraphicFeature(
                start=start, end=end, strand=strand, 
                color=color, label=label
            ))

    # Create the graphic record
    if record.annotations['topology'] == 'circular':
        graphic_record = CircularGraphicRecord(sequence_length=len(record.seq), features=features, top_position=1, labels_spacing=1500)
    else: 
        graphic_record = GraphicRecord(sequence_length=len(record.seq), features=features)

    # Plot the graphic record
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    graphic_record.plot(ax=ax, with_ruler=True, strand_in_label_threshold=7)

    # Set title
    plt.title(f"Sequence Map: {record.name}", fontsize=16, y=1.05)

    # Save to BytesIO object
    img_data = BytesIO()
    plt.savefig(img_data, format='png', dpi=300, bbox_inches="tight")
    img_data.seek(0)  # Move to the start of the BytesIO buffer

    return img_data, fig

# Helper function to save the figure to a file
def save_figure(fig, output_path):
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)