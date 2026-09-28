"""
PDB Position Coloring Script for ChimeraX (Diversity Pi Version)
Processes "AA_diversity.xlsx" and colors "VP1_Sabin_cristal_cleaned.pdb" 
based on Diversity_Pi values using a gradient.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex, Normalize
import os
import argparse

def load_diversity_data(excel_file):
    """Load and process the Diversity Excel file."""
    try:
        # Read Excel file assuming headers exist
        df = pd.read_excel(excel_file)
        
        # Check required columns
        required_cols = ['Position', 'Diversity_Pi']
        if not all(col in df.columns for col in required_cols):
            # Try to map if names are slightly different
            df.columns = [c.strip() for c in df.columns]
            if not all(col in df.columns for col in required_cols):
                print(f"Error: Excel file must contain columns: {required_cols}")
                print(f"Found columns: {df.columns.tolist()}")
                return None
        
        # Clean data
        df = df.dropna(subset=['Position', 'Diversity_Pi'])
        df['Position'] = pd.to_numeric(df['Position'], errors='coerce')
        df['Diversity_Pi'] = pd.to_numeric(df['Diversity_Pi'], errors='coerce')
        df = df.dropna(subset=['Position', 'Diversity_Pi'])
        df['Position'] = df['Position'].astype(int)
        
        # Sort by position
        df = df.sort_values('Position').reset_index(drop=True)
        
        # Generate colors based on Diversity_Pi
        # Use Viridis reversed colormap
        cmap = plt.get_cmap('viridis_r')
        
        # Normalize data for color mapping
        # We can map min-max of the data to 0-1, or fixed range
        min_pi = df['Diversity_Pi'].min()
        max_pi = df['Diversity_Pi'].max()
        
        print(f"Diversity_Pi Range: {min_pi:.6f} to {max_pi:.6f}")
        
        norm = Normalize(vmin=min_pi, vmax=max_pi)
        
        # Calculate hex colors
        def get_hex_color(val):
            # Get rgba from cmap
            rgba = cmap(norm(val))
            # Convert to hex
            return rgb2hex(rgba)
            
        df['Colour'] = df['Diversity_Pi'].apply(get_hex_color)
        
        print(f"Loaded {len(df)} positions from {excel_file}")
        return df
        
    except Exception as e:
        print(f"Error loading Excel file {excel_file}: {e}")
        return None

def hex_to_rgb(hex_color):
    """Convert hex color to RGB values (0-1 range)."""
    try:
        from matplotlib.colors import hex2color
        # Remove # if present
        hex_color = hex_color.lstrip('#')
        # Convert to RGB
        rgb = hex2color('#' + hex_color)
        return rgb
    except Exception as e:
        print(f"Warning: Could not convert color {hex_color} to RGB: {e}")
        return (0.5, 0.5, 0.5)  # Default gray

def hex_to_rgb_255(hex_color):
    """Convert hex color to RGB values (0-255 range) for ChimeraX."""
    rgb_01 = hex_to_rgb(hex_color)
    return tuple(int(c * 255) for c in rgb_01)

def create_position_coloring_commands(df_variant, model_id="1"):
    """Create ChimeraX commands to color individual positions."""
    coloring_commands = []
    
    for _, row in df_variant.iterrows():
        position = int(row['Position'])
        hex_color = row['Colour']
        
        # Convert hex to RGB (0-255 range for ChimeraX)
        rgb_255 = hex_to_rgb_255(hex_color)
        r, g, b = rgb_255
        
        # Create position selection for ChimeraX
        selection = f"#{model_id}:{position}"
        
        # Create ChimeraX coloring command
        coloring_commands.append(f"color {selection} rgb({r},{g},{b})")
    
    return coloring_commands

def read_pdb_file(pdb_file):
    """Read PDB file and return lines."""
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        return lines
    except Exception as e:
        print(f"Error reading PDB file: {e}")
        return None

def write_colored_pdb(pdb_lines, df_variant, variant_name, output_dir="colored_pdbs"):
    """Write a new PDB file with B-factor column modified to represent Diversity_Pi."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_file = os.path.join(output_dir, f"{variant_name}_colored.pdb")
    
    # Create mapping from position to Diversity_Pi
    position_to_value = dict(zip(df_variant['Position'], df_variant['Diversity_Pi']))
    
    with open(output_file, 'w') as f:
        for line in pdb_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    res_num = int(line[22:26].strip())
                    # Get the value for this position (or 0 if not filtered)
                    value = position_to_value.get(res_num, 0.0)
                    
                    # Replace B-factor (columns 61-66) with the value
                    # Start index 60 (0-based) is column 61 (1-based)
                    # B-factor is 6 chars wide (61-66)
                    new_line = line[:60] + f"{value:6.2f}" + line[66:]
                    f.write(new_line)
                except (ValueError, IndexError):
                    f.write(line)
            else:
                f.write(line)
    
    print(f"Colored PDB saved: {output_file}")
    return output_file

def create_chimerax_commands_file(variant_name, pdb_file, df_variant, coloring_commands, output_dir="colored_pdbs"):
    """Create a ChimeraX commands file for this variant."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    commands_file = os.path.join(output_dir, f"{variant_name}_chimerax_commands.cxc")
    
    with open(commands_file, 'w') as f:
        f.write(f"# ChimeraX commands for {variant_name} position coloring\n")
        f.write(f"# Generated automatically from Diversity_Pi data\n")
        f.write(f"# Total positions: {len(df_variant)}\n\n")
        
        # Load PDB file
        f.write("# Load PDB file\n")
        f.write(f"open {os.path.abspath(pdb_file)}\n\n")
        
        # Setup display
        f.write("# Setup display\n")
        f.write("hide #1 atoms\n")
        f.write("show #1 cartoons\n")
        f.write("set bgColor white\n")
        f.write("lighting soft\n\n")
        
        # Coloring commands
        f.write("# Color individual positions according to Excel data\n")
        for i, cmd in enumerate(coloring_commands):
            f.write(f"{cmd}\n")
        f.write("\n")
        
        # View setup
        f.write("# Setup view\n")
        f.write("view #1\n")
        f.write("center #1\n\n")
        
        # High quality rendering setup
        f.write("# High quality rendering settings\n")
        f.write("graphics silhouettes true\n")
        f.write("graphics silhouetteWidth 2\n")
        f.write("set maxFrameRate 60\n\n")
        
        # Save commands
        f.write("# Save high-quality image and session\n")
        f.write(f"save {variant_name}_position_colored.png supersample 3\n")
        f.write(f"save {variant_name}_session.cxs\n")
    
    print(f"ChimeraX commands file saved: {commands_file}")
    return commands_file

def create_color_reference(df_variant, variant_name, output_dir="colored_pdbs"):
    """Create a reference image showing the color scheme."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Create visualization
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10))
    
    # Top plot: Bar plot of Diversity_Pi
    positions = df_variant['Position'].values
    diversity = df_variant['Diversity_Pi'].values
    colors = [hex_to_rgb(color) for color in df_variant['Colour'].values]
    
    ax1.bar(positions, diversity, color=colors, edgecolor='black', linewidth=0.5)
    
    ax1.set_xlabel('Residue Position')
    ax1.set_ylabel('Diversity Pi')
    ax1.set_title(f'{variant_name} - Diversity Pi Color Mapping')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Set x-axis limits with some padding
    if len(positions) > 0:
        pos_range = max(positions) - min(positions)
        padding = max(5, pos_range * 0.05)
        ax1.set_xlim(min(positions) - padding, max(positions) + padding)
    
    # Bottom plot: Color legend
    ax2.axis('off')
    
    # Create color swatches for top 20 most diverse positions or sampled across range
    # Let's show a gradient bar essentially
    
    import matplotlib as mpl
    # Create a gradient image
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    
    ax2.imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis_r'), extent=[0, 10, 0, 1])
    ax2.set_yticks([])
    ax2.set_xlabel('Diversity Pi (Normalized by Dataset Min/Max)')
    
    # Add text for min/max labels
    min_pi = diversity.min()
    max_pi = diversity.max()
    ax2.text(0, -0.2, f"Min: {min_pi:.6f}", ha='center')
    ax2.text(10, -0.2, f"Max: {max_pi:.6f}", ha='center')
    
    plt.tight_layout()
    
    color_ref_file = os.path.join(output_dir, f"{variant_name}_chimerax_color_reference.png")
    plt.savefig(color_ref_file, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"ChimeraX color reference saved: {color_ref_file}")

def main():
    parser = argparse.ArgumentParser(description='Color PDB file based on Diversity Pi from Excel - ChimeraX version')
    parser.add_argument('--pdb_file', default='VP1_Sabin_cristal_cleaned.pdb', help='PDB file to color (default: VP1_Sabin_cristal_cleaned.pdb)')
    parser.add_argument('--input_file', default='AA_diversity.xlsx', help='Excel file with Diversity data (default: AA_diversity.xlsx)')
    parser.add_argument('--output_dir', default='colored_pdbs_chimerax', help='Output directory for colored PDB files')
    parser.add_argument('--model_id', default='1', help='ChimeraX model ID (default: 1)')
    
    args = parser.parse_args()
    
    # Load PDB file
    print(f"Loading PDB file: {args.pdb_file}")
    if os.path.exists(args.pdb_file):
        pdb_lines = read_pdb_file(args.pdb_file)
    else:
        print(f"Warning: PDB file {args.pdb_file} not found. Will generate commands without reading PDB.")
        pdb_lines = None
    
    # Find Excel file
    print(f"Loading Excel file: {args.input_file}")
    if not os.path.exists(args.input_file):
        print(f"Error: Excel file {args.input_file} not found.")
        return
        
    variant_name = "VP1_Sabin" # Default name
    
    # Load data
    df_variant = load_diversity_data(args.input_file)
    if df_variant is None or len(df_variant) == 0:
        print(f"No valid data found in {args.input_file}")
        return
    
    print(f"Processing {len(df_variant)} positions...")
    
    # Create coloring commands for ChimeraX
    coloring_commands = create_position_coloring_commands(df_variant, args.model_id)
    
    # Write colored PDB file (only if we have the PDB lines)
    colored_pdb_file = args.pdb_file
    if pdb_lines:
        colored_pdb_file = write_colored_pdb(pdb_lines, df_variant, variant_name, args.output_dir)
    
    # Create ChimeraX commands file
    create_chimerax_commands_file(variant_name, colored_pdb_file, df_variant, 
                                    coloring_commands, args.output_dir)
    
    # Create color reference
    create_color_reference(df_variant, variant_name, args.output_dir)
    
    print(f"\n{'='*80}")
    print("SUMMARY - CHIMERAX DIVERSITY COLORING:")
    print(f"{'='*80}")
    print(f"Output directory: '{args.output_dir}'")
    print(f"\nFiles generated:")
    if pdb_lines:
        print(f"- {variant_name}_colored.pdb (Diversity in B-factor)")
    print(f"- {variant_name}_chimerax_commands.cxc (ChimeraX script)")
    print(f"- {variant_name}_chimerax_color_reference.png (Legend)")
    print("\nTo use in ChimeraX:")
    print(f"open {variant_name}_chimerax_commands.cxc")

if __name__ == "__main__":
    main()
