#!/usr/bin/env python3
"""
Order Parameter (S_CD) - sn-2 Chain Comparison
===============================================
This script generates a single plot comparing sn-2 chain order parameters
for all 27 lipid membrane systems.

All sn-2 chains are 18:1 (oleoyl), so they all have the characteristic
double bond dip at C9-C10.

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
# COLOR SCHEME (Consistent with K_A plots)
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

carbons_18 = np.arange(2, 17)  # C2-C16 for 18-carbon chains (all sn-2 are 18:1)

# =============================================================================
# DATA: PC sn-2 Order Parameters for All 27 Systems
# All sn-2 chains are 18:1 (oleoyl)
# =============================================================================

sn2_data = {
    # =========================================================================
    # GROUP 1: PC/PA (80:20)
    # =========================================================================
    'Group1': {
        'DOPC_DOPA': {
            'data': [0.195, 0.186, 0.193, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'PC': 'DOPC', 'anionic': 'DOPA',
        },
        'DOPC_POPA': {
            'data': [0.195, 0.186, 0.193, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'PC': 'DOPC', 'anionic': 'POPA',
        },
        'DOPC_SOPA': {
            'data': [0.194, 0.186, 0.192, 0.169, 0.163, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'PC': 'DOPC', 'anionic': 'SOPA',
        },
        'POPC_DOPA': {
            'data': [0.197, 0.188, 0.193, 0.169, 0.162, 0.097, 0.083, 0.080, 0.075, 0.101, 0.102, 0.104, 0.091, 0.081, 0.061],
            'PC': 'POPC', 'anionic': 'DOPA',
        },
        'POPC_POPA': {
            'data': [0.197, 0.187, 0.192, 0.168, 0.161, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.060],
            'PC': 'POPC', 'anionic': 'POPA',
        },
        'POPC_SOPA': {
            'data': [0.197, 0.187, 0.192, 0.169, 0.162, 0.097, 0.083, 0.080, 0.075, 0.101, 0.102, 0.104, 0.091, 0.081, 0.061],
            'PC': 'POPC', 'anionic': 'SOPA',
        },
        'SOPC_DOPA': {
            'data': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'PC': 'SOPC', 'anionic': 'DOPA',
        },
        'SOPC_POPA': {
            'data': [0.200, 0.190, 0.195, 0.172, 0.165, 0.099, 0.086, 0.082, 0.079, 0.107, 0.109, 0.112, 0.100, 0.090, 0.069],
            'PC': 'SOPC', 'anionic': 'POPA',
        },
        'SOPC_SOPA': {
            'data': [0.200, 0.191, 0.196, 0.173, 0.166, 0.100, 0.086, 0.083, 0.078, 0.106, 0.108, 0.111, 0.098, 0.088, 0.068],
            'PC': 'SOPC', 'anionic': 'SOPA',
        },
    },
    
    # =========================================================================
    # GROUP 2: PC/PG (65:35)
    # =========================================================================
    'Group2': {
        'DOPC_DOPG': {
            'data': [0.194, 0.185, 0.192, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'PC': 'DOPC', 'anionic': 'DOPG',
        },
        'DOPC_POPG': {
            'data': [0.194, 0.186, 0.192, 0.170, 0.163, 0.097, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'PC': 'DOPC', 'anionic': 'POPG',
        },
        'DOPC_SOPG': {
            'data': [0.194, 0.185, 0.191, 0.169, 0.162, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'PC': 'DOPC', 'anionic': 'SOPG',
        },
        'POPC_DOPG': {
            'data': [0.197, 0.188, 0.193, 0.169, 0.162, 0.097, 0.084, 0.080, 0.076, 0.102, 0.103, 0.105, 0.092, 0.082, 0.062],
            'PC': 'POPC', 'anionic': 'DOPG',
        },
        'POPC_POPG': {
            'data': [0.198, 0.189, 0.194, 0.170, 0.163, 0.098, 0.084, 0.080, 0.076, 0.102, 0.103, 0.106, 0.093, 0.083, 0.063],
            'PC': 'POPC', 'anionic': 'POPG',
        },
        'POPC_SOPG': {
            'data': [0.198, 0.189, 0.194, 0.170, 0.164, 0.098, 0.084, 0.080, 0.076, 0.102, 0.104, 0.106, 0.093, 0.083, 0.063],
            'PC': 'POPC', 'anionic': 'SOPG',
        },
        'SOPC_DOPG': {
            'data': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'PC': 'SOPC', 'anionic': 'DOPG',
        },
        'SOPC_POPG': {
            'data': [0.200, 0.191, 0.196, 0.174, 0.167, 0.100, 0.087, 0.083, 0.079, 0.107, 0.109, 0.112, 0.100, 0.090, 0.069],
            'PC': 'SOPC', 'anionic': 'POPG',
        },
        'SOPC_SOPG': {
            'data': [0.200, 0.190, 0.196, 0.172, 0.166, 0.100, 0.086, 0.083, 0.079, 0.107, 0.108, 0.114, 0.100, 0.090, 0.070],
            'PC': 'SOPC', 'anionic': 'SOPG',
        },
    },
    
    # =========================================================================
    # GROUP 3: PC/PG (95:5)
    # =========================================================================
    'Group3': {
        'DOPC_DOPG': {
            'data': [0.194, 0.185, 0.192, 0.170, 0.163, 0.097, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'PC': 'DOPC', 'anionic': 'DOPG',
        },
        'DOPC_POPG': {
            'data': [0.194, 0.186, 0.192, 0.170, 0.163, 0.097, 0.081, 0.081, 0.076, 0.102, 0.102, 0.104, 0.090, 0.080, 0.061],
            'PC': 'DOPC', 'anionic': 'POPG',
        },
        'DOPC_SOPG': {
            'data': [0.194, 0.185, 0.191, 0.169, 0.162, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'PC': 'DOPC', 'anionic': 'SOPG',
        },
        'POPC_DOPG': {
            'data': [0.197, 0.187, 0.191, 0.167, 0.160, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.061],
            'PC': 'POPC', 'anionic': 'DOPG',
        },
        'POPC_POPG': {
            'data': [0.197, 0.186, 0.191, 0.167, 0.160, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.061],
            'PC': 'POPC', 'anionic': 'POPG',
        },
        'POPC_SOPG': {
            'data': [0.197, 0.187, 0.192, 0.168, 0.161, 0.096, 0.083, 0.079, 0.075, 0.101, 0.102, 0.103, 0.091, 0.081, 0.061],
            'PC': 'POPC', 'anionic': 'SOPG',
        },
        'SOPC_DOPG': {
            'data': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'PC': 'SOPC', 'anionic': 'DOPG',
        },
        'SOPC_POPG': {
            'data': [0.201, 0.191, 0.196, 0.172, 0.165, 0.099, 0.087, 0.083, 0.080, 0.108, 0.111, 0.113, 0.102, 0.091, 0.070],
            'PC': 'SOPC', 'anionic': 'POPG',
        },
        'SOPC_SOPG': {
            'data': [0.200, 0.190, 0.196, 0.172, 0.165, 0.100, 0.088, 0.084, 0.081, 0.109, 0.112, 0.115, 0.103, 0.092, 0.071],
            'PC': 'SOPC', 'anionic': 'SOPG',
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

def plot_all_27_sn2_single():
    """
    Create a single plot comparing sn-2 order parameters for all 27 systems.
    Color = PC type, Line style = Group, Marker = Anionic type
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Plot each system
    for group in ['Group1', 'Group2', 'Group3']:
        for sys_name, sys_data in sn2_data[group].items():
            pc = sys_data['PC']
            anionic = sys_data['anionic']
            data = sys_data['data']
            
            # Get anionic type (DO, PO, SO)
            anionic_type = anionic[:2]
            
            # Create label
            label = f'{pc}/{anionic} ({GROUP_LABELS[group]})'
            
            # Plot
            ax.plot(carbons_18, data,
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
    ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow', label='Double bond region (C9-C10)')
    
    # Formatting
    ax.set_xlabel('Carbon Position', fontsize=14, fontweight='bold')
    ax.set_ylabel('$|S_{CD}|$', fontsize=14, fontweight='bold')
    ax.set_title('PC Lipid sn-2 Chain Order Parameters (All 18:1 Oleoyl)\nAll 27 Systems Compared', 
                 fontsize=16, fontweight='bold', pad=15)
    ax.set_xlim(1, 17)
    ax.set_ylim(0.04, 0.22)
    ax.set_xticks(range(2, 17))
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Create custom legend
    from matplotlib.lines import Line2D
    
    legend_elements = [
        # PC colors
        Line2D([0], [0], color=PC_COLORS['DOPC'], linewidth=3, label='DOPC sn-2'),
        Line2D([0], [0], color=PC_COLORS['POPC'], linewidth=3, label='POPC sn-2'),
        Line2D([0], [0], color=PC_COLORS['SOPC'], linewidth=3, label='SOPC sn-2'),
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
    ax.annotate('Double bond dip\n(C9-C10)', 
                xy=(9.5, 0.076), xytext=(12.5, 0.055),
                fontsize=10, ha='center',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    
    ax.annotate('SOPC: Highest order\n(saturated sn-1 neighbor)', 
                xy=(5, 0.173), xytext=(2, 0.20),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='#9b59b6', lw=1.5),
                color='#9b59b6')
    
    ax.annotate('DOPC: Lowest order\n(unsaturated sn-1 neighbor)', 
                xy=(5, 0.163), xytext=(7, 0.145),
                fontsize=9, ha='left',
                arrowprops=dict(arrowstyle='->', color='#3498db', lw=1.5),
                color='#3498db')
    
    ax.annotate('All sn-2 chains are 18:1 (oleoyl)', 
                xy=(14, 0.045), fontsize=10, ha='center',
                style='italic', color='gray')
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn2_Single_Plot.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn2_Single_Plot.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn2_Single_Plot.png/pdf")
    plt.show()


def plot_all_27_sn2_by_pc():
    """
    Create three subplots, one for each PC lipid type, showing all 9 systems each.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    pc_lipids = ['DOPC', 'POPC', 'SOPC']
    
    for ax, pc in zip(axes, pc_lipids):
        # Plot all 9 systems for this PC
        for group in ['Group1', 'Group2', 'Group3']:
            for sys_name, sys_data in sn2_data[group].items():
                if sys_data['PC'] == pc:
                    anionic = sys_data['anionic']
                    data = sys_data['data']
                    anionic_type = anionic[:2]
                    
                    label = f'{anionic} ({GROUP_LABELS[group]})'
                    
                    ax.plot(carbons_18, data,
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
        
        # Highlight double bond region
        ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        
        ax.set_xlabel('Carbon Position', fontsize=12, fontweight='bold')
        ax.set_title(f'{pc} sn-2 (18:1)\n9 Systems', fontsize=13, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0.04, 0.22)
        ax.set_xticks(range(2, 17, 2))
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, ncol=1, framealpha=0.9)
    
    axes[0].set_ylabel('$|S_{CD}|$', fontsize=13, fontweight='bold')
    
    fig.suptitle('PC Lipid sn-2 Chain Order Parameters by PC Type (All 18:1)\nAll 27 Systems (9 per panel)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn2_By_PC.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn2_By_PC.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn2_By_PC.png/pdf")
    plt.show()


def plot_all_27_sn2_by_group():
    """
    Create three subplots, one for each group, showing all 9 systems each.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    
    groups = ['Group1', 'Group2', 'Group3']
    group_titles = ['Group 1: PC/PA (80:20)', 'Group 2: PC/PG (65:35)', 'Group 3: PC/PG (95:5)']
    
    for ax, group, title in zip(axes, groups, group_titles):
        for sys_name, sys_data in sn2_data[group].items():
            pc = sys_data['PC']
            anionic = sys_data['anionic']
            data = sys_data['data']
            anionic_type = anionic[:2]
            
            label = f'{pc}/{anionic}'
            
            ax.plot(carbons_18, data,
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
    
    fig.suptitle('PC Lipid sn-2 Chain Order Parameters by Group (All 18:1)\nAll 27 Systems (9 per panel)', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('Order_All_27_Systems_sn2_By_Group.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_sn2_By_Group.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_sn2_By_Group.png/pdf")
    plt.show()


def plot_pc_comparison_averaged():
    """
    Create a plot comparing average sn-2 order parameters for each PC type.
    Shows mean ± std across all 9 systems for each PC.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Calculate averages for each PC type
    for pc in ['DOPC', 'POPC', 'SOPC']:
        all_data = []
        for group in ['Group1', 'Group2', 'Group3']:
            for sys_name, sys_data in sn2_data[group].items():
                if sys_data['PC'] == pc:
                    all_data.append(sys_data['data'])
        
        all_data = np.array(all_data)
        mean = np.mean(all_data, axis=0)
        std = np.std(all_data, axis=0)
        
        ax.plot(carbons_18, mean, 'o-', color=PC_COLORS[pc], linewidth=2.5,
                markersize=8, label=f'{pc} (n=9)')
        ax.fill_between(carbons_18, mean - std, mean + std, 
                        color=PC_COLORS[pc], alpha=0.2)
    
    # Highlight double bond region
    ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow', label='Double bond region')
    
    ax.set_xlabel('Carbon Position', fontsize=14, fontweight='bold')
    ax.set_ylabel('$|S_{CD}|$', fontsize=14, fontweight='bold')
    ax.set_title('PC Lipid sn-2 Chain Order Parameters\nMean ± SD Across All Systems (n=9 each)', 
                 fontsize=15, fontweight='bold')
    ax.set_xlim(1, 17)
    ax.set_ylim(0.04, 0.22)
    ax.set_xticks(range(2, 17))
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=12, framealpha=0.95)
    
    # Add annotation
    ax.annotate('SOPC > POPC > DOPC\n(Effect of saturated sn-1 neighbor)', 
                xy=(4, 0.19), fontsize=10, ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig('Order_sn2_PC_Comparison_Averaged.png', dpi=300, facecolor='white')
    plt.savefig('Order_sn2_PC_Comparison_Averaged.pdf', facecolor='white')
    print("Saved: Order_sn2_PC_Comparison_Averaged.png/pdf")
    plt.show()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("sn-2 Order Parameter Comparison - All 27 Systems")
    print("=" * 70)
    print("\nGenerating plots...\n")
    
    # Plot 1: All 27 systems on single plot
    print("-" * 50)
    print("1. All 27 Systems on Single Plot")
    print("-" * 50)
    plot_all_27_sn2_single()
    
    # Plot 2: Separated by PC type
    print("\n" + "-" * 50)
    print("2. Separated by PC Type (3 panels)")
    print("-" * 50)
    plot_all_27_sn2_by_pc()
    
    # Plot 3: Separated by group
    print("\n" + "-" * 50)
    print("3. Separated by Group (3 panels)")
    print("-" * 50)
    plot_all_27_sn2_by_group()
    
    # Plot 4: PC comparison with averages
    print("\n" + "-" * 50)
    print("4. PC Type Comparison (Averaged)")
    print("-" * 50)
    plot_pc_comparison_averaged()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
    print("\nOutput files:")
    print("  1. Order_All_27_Systems_sn2_Single_Plot.png/pdf - All 27 on one plot")
    print("  2. Order_All_27_Systems_sn2_By_PC.png/pdf       - Separated by PC type")
    print("  3. Order_All_27_Systems_sn2_By_Group.png/pdf    - Separated by group")
    print("  4. Order_sn2_PC_Comparison_Averaged.png/pdf     - Averaged by PC type")
    print("\nLegend Key:")
    print("  Colors:      DOPC=Blue, POPC=Green, SOPC=Purple")
    print("  Line styles: Solid=PA(80:20), Dashed=PG(65:35), Dotted=PG(95:5)")
    print("  Markers:     Circle=DO-, Square=PO-, Triangle=SO-")
    print("\nKey Finding: All sn-2 chains are 18:1 (oleoyl)")
    print("  - All show characteristic double bond dip at C9-C10")
    print("  - SOPC > POPC > DOPC order (saturated sn-1 neighbor effect)")
    print("=" * 70)
