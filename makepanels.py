import svgutils.transform as sg
from pathlib import Path

# Set your base directory
base_dir = Path("C:/Users/hayat/Downloads/R_files/graphs/publication_figures/Figure_6")

# List your files in order
files = [
    base_dir / "all_plasmids_cog_category_plot_filtered_Figure_6a.svg",
    base_dir / "widespread_plasmids_cog_category_plot_filtered_Figure_6b.svg",
    base_dir / "all_plasmids_top_pfam_plot_categorized_Figure_6c.svg",
    base_dir / "widespread_plasmids_top_pfam_plot_categorized_Figure_6d.svg",
    base_dir / "enriched_PFAMs_OR_with_CI_Figure_6e.svg"

]

# Panel dimensions
panel_width = 800
panel_height = 500


# Grid dimensions (3 columns × 2 rows)
ncols = 3
nrows = 3
fig_width = ncols * panel_width
fig_height = nrows * panel_height

# Create new figure
fig = sg.SVGFigure(f"{fig_width}px", f"{fig_height}px")

elements = []
labels = []

for i, f in enumerate(files):
    # Load panel
    panel = sg.fromfile(str(f)).getroot()

    # Wrap inside a GroupElement (so moveto works reliably)
    group = sg.GroupElement([panel])

    # Determine row/col
    row = i // ncols
    col = i % ncols

    # Move panel into grid position
    group.moveto(col * panel_width, row * panel_height)
    elements.append(group)

    # Add panel label
    label_char = chr(ord('a') + i)
    label = sg.TextElement(col * panel_width + 10, 
                           row * panel_height + 20,
                           f"({label_char})", 
                           size=16, weight="bold")
    labels.append(label)

# Append panels + labels
fig.append(elements + labels)

# Save into the same folder (or change path if you want elsewhere)
out_file = base_dir / "figure6.svg"
fig.save(str(out_file))
