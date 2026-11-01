#!/usr/bin/env python3
"""
Voronoi Diagram Analysis for Lipid Bilayers
============================================
Generates representative Voronoi diagrams for top and bottom leaflets
of binary lipid mixtures from MD simulations.

For Jupyter Notebook / IPython use
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
import MDAnalysis as mda

# =============================================================================
# COLOR SCHEME
# =============================================================================

PC_COLORS = {
    'DOPC': '#3498db',  # Blue
    'POPC': '#2ecc71',  # Green  
    'SOPC': '#9b59b6',  # Purple
}

SECONDARY_COLORS = {
    'DOPA': '#e74c3c',  # Red
    'POPA': '#e67e22',  # Orange
    'SOPA': '#f39c12',  # Yellow-Orange
    'DOPG': '#e74c3c',  # Red
    'POPG': '#e67e22',  # Orange
    'SOPG': '#f39c12',  # Yellow-Orange
}

# =============================================================================
# FUNCTION: Extract Leaflet Positions
# =============================================================================

def extract_leaflet_positions(universe, pc_lipid, secondary_lipid, frame=-1):
    """
    Extract phosphorus positions for each lipid type, separated by leaflet.
    """
    universe.trajectory[frame]
    
    pc_P = universe.select_atoms(f"resname {pc_lipid} and name P")
    sec_P = universe.select_atoms(f"resname {secondary_lipid} and name P")
    all_P = universe.select_atoms("name P")
    
    membrane_center_z = all_P.positions[:, 2].mean()
    
    pc_upper = pc_P.positions[pc_P.positions[:, 2] > membrane_center_z][:, :2]
    pc_lower = pc_P.positions[pc_P.positions[:, 2] < membrane_center_z][:, :2]
    sec_upper = sec_P.positions[sec_P.positions[:, 2] > membrane_center_z][:, :2]
    sec_lower = sec_P.positions[sec_P.positions[:, 2] < membrane_center_z][:, :2]
    
    box = universe.dimensions[:2]
    
    return {
        'upper': {'pc': pc_upper, 'secondary': sec_upper},
        'lower': {'pc': pc_lower, 'secondary': sec_lower},
        'box': box,
        'frame': universe.trajectory.frame,
        'time': universe.trajectory.time,
    }

# =============================================================================
# FUNCTION: Plot Voronoi Leaflet
# =============================================================================

def plot_voronoi_leaflet(ax, pc_positions, sec_positions, box, pc_lipid, secondary_lipid, title=""):
    """
    Create a Voronoi diagram for one leaflet.
    """
    box_x, box_y = box
    
    all_positions = np.vstack([pc_positions, sec_positions])
    n_pc = len(pc_positions)
    
    vor = Voronoi(all_positions)
    
    pc_color = PC_COLORS.get(pc_lipid, '#3498db')
    sec_color = SECONDARY_COLORS.get(secondary_lipid, '#e74c3c')
    
    for i, region_idx in enumerate(vor.point_region):
        region = vor.regions[region_idx]
        if not region or -1 in region:
            continue
        
        polygon = np.array([vor.vertices[j] for j in region])
        
        if np.any(polygon < -10) or np.any(polygon[:, 0] > box_x + 10) or np.any(polygon[:, 1] > box_y + 10):
            continue
        
        color = pc_color if i < n_pc else sec_color
        ax.fill(*zip(*polygon), color=color, alpha=0.7, edgecolor='white', linewidth=0.5)
    
    ax.scatter(pc_positions[:, 0], pc_positions[:, 1], c='white', s=20, 
               edgecolors='black', linewidths=0.5, zorder=5)
    ax.scatter(sec_positions[:, 0], sec_positions[:, 1], c='black', s=20,
               edgecolors='white', linewidths=0.5, zorder=5)
    
    ax.set_xlim(0, box_x)
    ax.set_ylim(0, box_y)
    ax.set_aspect('equal')
    ax.set_xlabel('X (Å)', fontsize=11)
    ax.set_ylabel('Y (Å)', fontsize=11)
    ax.set_title(title, fontweight='bold', fontsize=12)
    
    return ax

# =============================================================================
# FUNCTION: Main Voronoi Analysis
# =============================================================================

def analyze_voronoi(tpr_file, xtc_file, pc_lipid, secondary_lipid, output_prefix, frame=-1):
    """
    Main function to generate Voronoi diagrams for both leaflets.
    """
    print(f"Loading trajectory: {xtc_file}")
    u = mda.Universe(tpr_file, xtc_file)
    
    print(f"Extracting positions for frame {frame}...")
    data = extract_leaflet_positions(u, pc_lipid, secondary_lipid, frame)
    
    print(f"  Upper leaflet: {len(data['upper']['pc'])} {pc_lipid}, {len(data['upper']['secondary'])} {secondary_lipid}")
    print(f"  Lower leaflet: {len(data['lower']['pc'])} {pc_lipid}, {len(data['lower']['secondary'])} {secondary_lipid}")
    print(f"  Box: {data['box'][0]:.1f} x {data['box'][1]:.1f} Å")
    print(f"  Time: {data['time']} ps")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    plot_voronoi_leaflet(axes[0], 
                         data['upper']['pc'], 
                         data['upper']['secondary'],
                         data['box'],
                         pc_lipid, secondary_lipid,
                         title=f'Upper Leaflet\n({len(data["upper"]["pc"])} {pc_lipid} + {len(data["upper"]["secondary"])} {secondary_lipid})')
    
    plot_voronoi_leaflet(axes[1],
                         data['lower']['pc'],
                         data['lower']['secondary'],
                         data['box'],
                         pc_lipid, secondary_lipid,
                         title=f'Lower Leaflet\n({len(data["lower"]["pc"])} {pc_lipid} + {len(data["lower"]["secondary"])} {secondary_lipid})')
    
    pc_color = PC_COLORS.get(pc_lipid, '#3498db')
    sec_color = SECONDARY_COLORS.get(secondary_lipid, '#e74c3c')
    
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=pc_color, edgecolor='white', label=pc_lipid, alpha=0.7),
        Patch(facecolor=sec_color, edgecolor='white', label=secondary_lipid, alpha=0.7),
    ]
    fig.legend(handles=legend_elements, loc='upper center', ncol=2, fontsize=12,
               bbox_to_anchor=(0.5, 0.98))
    
    plt.suptitle(f'Voronoi Diagram: {pc_lipid}/{secondary_lipid} Binary Mixture\n(Frame {data["frame"]}, t = {data["time"]/1000:.0f} ns)',
                 fontweight='bold', fontsize=14, y=1.05)
    
    plt.tight_layout()
    
    plt.savefig(f'{output_prefix}_voronoi.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f'{output_prefix}_voronoi.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    print(f"\nSaved: {output_prefix}_voronoi.png/pdf")
    
    plt.show()
    
    return data


# =============================================================================
# RUN ANALYSIS - MODIFY THESE PARAMETERS FOR YOUR SYSTEM
# =============================================================================

# Set your file paths and lipid names here:
tpr_file = "Prod_0.tpr"
xtc_file = "Prod_merged.xtc"
pc_lipid = "DOPC"
secondary_lipid = "DOPA"
output_prefix = "DOPC_DOPA"

# Run the analysis
data = analyze_voronoi(tpr_file, xtc_file, pc_lipid, secondary_lipid, output_prefix)
