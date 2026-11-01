#!/usr/bin/env python3
"""
Order Parameter (S_CD) - sn-1 Chain Comparison
===============================================
This script generates a single plot comparing sn-1 chain order parameters
for all 27 lipid membrane systems.

Systems are organized by PC lipid type with different line styles for groups.

Author: Generated for lipid membrane MD simulation analysis
Date: February 2026
"""

import matplotlib.pyplot as plt
import numpy as np

# =============================================================================
# PLOT SETTINGS
# =============================================================================

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# =============================================================================
# COLOR SCHEME
# =============================================================================

# Colors by PC lipid type
PC_COLORS = {
    'DOPC': '#3498db',  # Blue
    'POPC': '#2ecc71',  # Green
    'SOPC': '#9b59b6',  # Purple
}

# Line styles by group
GROUP_LINESTYLES = {
    'Group1': '-',    # Solid for PA (80:20)
    'Group2': '--',   # Dashed for PG (65:35)
    'Group3': ':',    # Dotted for PG (95:5)
}

# Markers by anionic lipid type
ANIONIC_MARKERS = {
    'DO': 'o',  # Circle for DOPA/DOPG
    'PO': 's',  # Square for POPA/POPG
    'SO': '^',  # Triangle for SOPA/SOPG
}

# =============================================================================
# DATA: Carbon positions
# =============================================================================

carbons_18 = np.arange(2, 17)  # C2-C16 for 18-carbon chains
carbons_16 = np.arange(2, 15)  # C2-C14 for 16-carbon chains

# =============================================================================
# DATA: PC sn-1 Order Parameters for All 27 Systems
# =============================================================================

# Organized by group and system
sn1_data = {
    # =========================================================================
    # GROUP 1: PC/PA (80:20)
    # =========================================================================
    'Group1': {
        'DOPC_DOPA': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'DOPA',
        },
        'DOPC_POPA': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'POPA',
        },
        'DOPC_SOPA': {
            'data': [0.144, 0.167, 0.165, 0.163, 0.148, 0.107, 0.073, 0.071, 0.067, 0.091, 0.088, 0.090, 0.075, 0.067, 0.050],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'SOPA',
        },
        'POPC_DOPA': {
            'data': [0.151, 0.185, 0.183, 0.194, 0.184, 0.181, 0.165, 0.157, 0.138, 0.124, 0.104, 0.087, 0.065],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'DOPA',
        },
        'POPC_POPA': {
            'data': [0.151, 0.184, 0.182, 0.193, 0.183, 0.180, 0.164, 0.156, 0.137, 0.123, 0.103, 0.086, 0.064],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'POPA',
        },
        'POPC_SOPA': {
            'data': [0.152, 0.185, 0.183, 0.194, 0.184, 0.181, 0.165, 0.157, 0.138, 0.124, 0.104, 0.087, 0.065],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'SOPA',
        },
        'SOPC_DOPA': {
            'data': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'DOPA',
        },
        'SOPC_POPA': {
            'data': [0.151, 0.185, 0.184, 0.196, 0.188, 0.188, 0.174, 0.165, 0.148, 0.136, 0.117, 0.104, 0.087, 0.073, 0.054],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'POPA',
        },
        'SOPC_SOPA': {
            'data': [0.152, 0.186, 0.186, 0.198, 0.190, 0.188, 0.175, 0.166, 0.148, 0.135, 0.117, 0.103, 0.086, 0.072, 0.054],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'SOPA',
        },
    },
    
    # =========================================================================
    # GROUP 2: PC/PG (65:35)
    # =========================================================================
    'Group2': {
        'DOPC_DOPG': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'DOPG',
        },
        'DOPC_POPG': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'POPG',
        },
        'DOPC_SOPG': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'SOPG',
        },
        'POPC_DOPG': {
            'data': [0.152, 0.186, 0.184, 0.195, 0.185, 0.182, 0.166, 0.158, 0.139, 0.125, 0.105, 0.088, 0.066],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'DOPG',
        },
        'POPC_POPG': {
            'data': [0.153, 0.187, 0.185, 0.197, 0.187, 0.184, 0.168, 0.163, 0.143, 0.128, 0.108, 0.090, 0.067],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'POPG',
        },
        'POPC_SOPG': {
            'data': [0.154, 0.188, 0.186, 0.198, 0.188, 0.186, 0.170, 0.166, 0.147, 0.132, 0.111, 0.092, 0.068],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'SOPG',
        },
        'SOPC_DOPG': {
            'data': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'DOPG',
        },
        'SOPC_POPG': {
            'data': [0.152, 0.186, 0.186, 0.199, 0.190, 0.189, 0.176, 0.170, 0.152, 0.140, 0.121, 0.107, 0.090, 0.076, 0.056],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'POPG',
        },
        'SOPC_SOPG': {
            'data': [0.152, 0.186, 0.186, 0.197, 0.189, 0.188, 0.175, 0.168, 0.149, 0.137, 0.118, 0.109, 0.088, 0.074, 0.058],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'SOPG',
        },
    },
    
    # =========================================================================
    # GROUP 3: PC/PG (95:5)
    # =========================================================================
    'Group3': {
        'DOPC_DOPG': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'DOPG',
        },
        'DOPC_POPG': {
            'data': [0.144, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.071, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'POPG',
        },
        'DOPC_SOPG': {
            'data': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'chain': '18:1', 'PC': 'DOPC', 'anionic': 'SOPG',
        },
        'POPC_DOPG': {
            'data': [0.151, 0.184, 0.182, 0.192, 0.182, 0.179, 0.163, 0.153, 0.134, 0.120, 0.101, 0.085, 0.063],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'DOPG',
        },
        'POPC_POPG': {
            'data': [0.151, 0.184, 0.182, 0.193, 0.183, 0.180, 0.164, 0.154, 0.135, 0.121, 0.101, 0.085, 0.063],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'POPG',
        },
        'POPC_SOPG': {
            'data': [0.152, 0.185, 0.184, 0.194, 0.184, 0.181, 0.165, 0.156, 0.138, 0.124, 0.104, 0.085, 0.061],
            'chain': '16:0', 'PC': 'POPC', 'anionic': 'SOPG',
        },
        'SOPC_DOPG': {
            'data': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'DOPG',
        },
        'SOPC_POPG': {
            'data': [0.150, 0.184, 0.184, 0.196, 0.189, 0.189, 0.175, 0.167, 0.150, 0.138, 0.120, 0.107, 0.089, 0.074, 0.055],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'POPG',
        },
        'SOPC_SOPG': {
            'data': [0.151, 0.186, 0.185, 0.197, 0.189, 0.189, 0.176, 0.167, 0.150, 0.139, 0.120, 0.107, 0.089, 0.075, 0.056],
            'chain': '18:0', 'PC': 'SOPC', 'anionic': 'SOPG',
        },
    },
}

GROUP_LABELS = {
    'Group1': 'PA (80:20)',
    'Group2': 'PG (65:35)',
    'Group3': 'PG (95:5)',
}


# =============================================================================
# PLOTTING FUNCTION: Single Plot with All 27 Systems
# =============================================================================

def plot_all_27_sn1_single():
    """
    Create a single plot comparing sn-1 order parameters for all 27 systems.
    Color = PC type, Line style = Group, Marker = Anionic type
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Plot each system
    for group in ['Group1', 'Group2', 'Group3']:
        for sys_name, sys_data in sn1_data[group].items():
            pc = sys_data['PC']
            anionic = sys_data['anionic']
            chain = sys_data['chain']
            data = sys_data['data']
            
            # Determine carbons based on chain type
            if chain == '16:0':
                carbons = carbons_16
            else:
                carbons = carbons_18
            
            # Get anionic type (DO, PO, SO)
            anionic_type = anionic[:2]
            
            # Create label
            label = f'{pc}/{anionic} ({GROUP_LABELS[group]})'
            
            # Plot
            ax.plot(carbons, data,
                    linestyle=GROUP_LINESTYLES[group],
                    marker=ANIONIC_MARKERS[anionic_type],
                    color=PC_COLORS[pc],
                    linewidth=1.5,
                    markersize=5,
                    markerfacecolor=PC_COLORS[pc],
                    markeredgecolor='white',
                    markeredgewidth=0.5,
                    alpha=0.8,
                    label=label)
    
    # Highlight double bond region
    ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow', label='Double bond region (18:1)')
    
    # Formatting
    ax.set_xlabel('Carbon Position', fontsize=14, fontweight='bold')
    ax.set_ylabel('$|S_{CD}|$', fontsize=14, fontweight='bold')
    ax.set_title('PC Lipid sn-1 Chain Order Parameters\nAll 27 Systems Compared', 
                 fontsize=16, fontweight='bold', pad=15)
    ax.set_xlim(1, 17)
    ax.set_ylim(0.04, 0.22)
    ax.set_xticks(range(2, 17))
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Create custom legend
    from matplotlib.lines import Line2D
    
    legend_elements = [
        # PC colors
        Line2D([0], [0], color=PC_COLORS['DOPC'], linewidth=3, label='DOPC (18:1)'),
        Line2D([0], [0], color=PC_COLORS['POPC'], linewidth=3, label='POPC (16:0)'),
        Line2D([0], [0], color=PC_COLORS['SOPC'], linewidth=3, label='SOPC (18:0)'),
        # Spacer
        Line2D([0], [0], color='white', linewidth=0, label=''),
        # Group line styles
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2, label='Group 1: PA (80:20)'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=2, label='Group 2: PG (65:35)'),
        Line2D([0], [0], color='gray', linestyle=':', linewidth=2, label='Group 3: PG (95:5)'),
        # Spacer
        Line2D([0], [0], color='white', linewidth=0, label=''),
        # Anionic markers
        Line2D([0], [0], color='gray', marker='o', linestyle='None', markersize=8, label='DO- anionic'),
        Line2D([0], [0], color='gray', marker='s', linestyle='None', markersize=8, label='PO- anionic'),
        Line2D([0], [0], color='gray', marker='^', linestyle='None', markersize=8, label='SO- anionic'),
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, 
              framealpha=0.95, edgecolor='black', ncol=1)
    
    # Add annotations
    ax.annotate('DOPC: 18:1 (unsaturated)\nDouble bond dip at C9-C10', 
                xy=(9.5, 0.068), xytext=(12, 0.055),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='#3498db', lw=1.5),
                color='#3498db')
    
    ax.annotate('POPC: 16:0 (saturated)\nMonotonic decrease', 
                xy=(8, 0.175), xytext=(11, 0.19),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='#2ecc71', lw=1.5),
                color='#2ecc71')
    
    ax.annotate('SOPC: 18:0 (saturated)\nHighest order', 
                xy=(6, 0.189), xytext=(2, 0.21),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='#9b59b6', lw=1.5),
                color='#9b59b6')
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn1_Single_Plot.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn1_Single_Plot.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn1_Single_Plot.png/pdf")
    plt.show()


def plot_all_27_sn1_by_pc():
    """
    Create three subplots, one for each PC lipid type, showing all 9 systems each.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    pc_lipids = ['DOPC', 'POPC', 'SOPC']
    chain_types = {'DOPC': '18:1', 'POPC': '16:0', 'SOPC': '18:0'}
    
    for ax, pc in zip(axes, pc_lipids):
        chain = chain_types[pc]
        carbons = carbons_16 if chain == '16:0' else carbons_18
        
        # Plot all 9 systems for this PC
        for group in ['Group1', 'Group2', 'Group3']:
            for sys_name, sys_data in sn1_data[group].items():
                if sys_data['PC'] == pc:
                    anionic = sys_data['anionic']
                    data = sys_data['data']
                    anionic_type = anionic[:2]
                    
                    label = f'{anionic} ({GROUP_LABELS[group]})'
                    
                    ax.plot(carbons, data,
                            linestyle=GROUP_LINESTYLES[group],
                            marker=ANIONIC_MARKERS[anionic_type],
                            color=PC_COLORS[pc],
                            linewidth=2,
                            markersize=6,
                            markerfacecolor=PC_COLORS[pc],
                            markeredgecolor='white',
                            markeredgewidth=0.5,
                            alpha=0.85,
                            label=label)
        
        # Highlight double bond region (only for 18:1)
        if chain == '18:1':
            ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        
        ax.set_xlabel('Carbon Position', fontsize=12, fontweight='bold')
        ax.set_title(f'{pc} sn-1 ({chain})\n9 Systems', fontsize=13, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0.04, 0.22)
        ax.set_xticks(range(2, 17, 2))
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, ncol=1, framealpha=0.9)
    
    axes[0].set_ylabel('$|S_{CD}|$', fontsize=13, fontweight='bold')
    
    fig.suptitle('PC Lipid sn-1 Chain Order Parameters by PC Type\nAll 27 Systems (9 per panel)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn1_By_PC.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn1_By_PC.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn1_By_PC.png/pdf")
    plt.show()


def plot_all_27_sn1_by_group():
    """
    Create three subplots, one for each group, showing all 9 systems each.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    groups = ['Group1', 'Group2', 'Group3']
    group_titles = ['Group 1: PC/PA (80:20)', 'Group 2: PC/PG (65:35)', 'Group 3: PC/PG (95:5)']
    
    for ax, group, title in zip(axes, groups, group_titles):
        for sys_name, sys_data in sn1_data[group].items():
            pc = sys_data['PC']
            anionic = sys_data['anionic']
            chain = sys_data['chain']
            data = sys_data['data']
            
            carbons = carbons_16 if chain == '16:0' else carbons_18
            anionic_type = anionic[:2]
            
            label = f'{pc}/{anionic}'
            
            ax.plot(carbons, data,
                    linestyle='-',
                    marker=ANIONIC_MARKERS[anionic_type],
                    color=PC_COLORS[pc],
                    linewidth=2,
                    markersize=6,
                    markerfacecolor=PC_COLORS[pc],
                    markeredgecolor='white',
                    markeredgewidth=0.5,
                    alpha=0.85,
                    label=label)
        
        # Highlight double bond region
        ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        
        ax.set_xlabel('Carbon Position', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=13, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0.04, 0.22)
        ax.set_xticks(range(2, 17, 2))
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, ncol=1, framealpha=0.9)
    
    axes[0].set_ylabel('$|S_{CD}|$', fontsize=13, fontweight='bold')
    
    fig.suptitle('PC Lipid sn-1 Chain Order Parameters by Group\nAll 27 Systems (9 per panel)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn1_By_Group.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn1_By_Group.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn1_By_Group.png/pdf")
    plt.show()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("sn-1 Order Parameter Comparison - All 27 Systems")
    print("=" * 70)
    print("\nGenerating plots...\n")
    
    # Plot 1: All 27 systems on a single plot
    print("-" * 50)
    print("1. All 27 Systems on Single Plot")
    print("-" * 50)
    plot_all_27_sn1_single()
    
    # Plot 2: Separated by PC type
    print("\n" + "-" * 50)
    print("2. Separated by PC Type (3 panels)")
    print("-" * 50)
    plot_all_27_sn1_by_pc()
    
    # Plot 3: Separated by group
    print("\n" + "-" * 50)
    print("3. Separated by Group (3 panels)")
    print("-" * 50)
    plot_all_27_sn1_by_group()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
    print("\nOutput files:")
    print("  1. Order_All_27_Systems_sn1_Single_Plot.png/pdf - All 27 on one plot")
    print("  2. Order_All_27_Systems_sn1_By_PC.png/pdf       - Separated by PC type")
    print("  3. Order_All_27_Systems_sn1_By_Group.png/pdf    - Separated by group")
    print("\nLegend Key:")
    print("  Colors:      DOPC=Blue (18:1), POPC=Green (16:0), SOPC=Purple (18:0)")
    print("  Line styles: Solid=PA(90:10), Dashed=PG(60:40), Dotted=PG(90:10)")
    print("  Markers:     Circle=DO-, Square=PO-, Triangle=SO-")
    print("=" * 70)
