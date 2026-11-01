#!/usr/bin/env python3
"""
Decomposed Total Interaction Energies in Mixed Lipid Bilayers
==============================================================
Generates 4 publication-quality figures comparing:
(a) Secondary-Secondary interactions (PA-PA or PG-PG)
(b) PC-PC interactions
(c) PC-Secondary interactions (PC-PA or PC-PG)
(d) Total lipid-lipid interactions

Using CONSISTENT COLOR CODING:
- DOPC: Blue (#3498db)
- POPC: Green (#2ecc71)
- SOPC: Purple (#9b59b6)

Author: Generated for lipid membrane MD simulation analysis
Date: February 2026
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# =============================================================================
# PLOT SETTINGS
# =============================================================================

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# =============================================================================
# CONSISTENT COLOR SCHEME
# =============================================================================

PC_COLORS = {
    'DOPC': '#3498db',  # Blue
    'POPC': '#2ecc71',  # Green
    'SOPC': '#9b59b6',  # Purple
}

# =============================================================================
# DATA: All 27 Systems - Decomposed Interaction Energies (kJ/mol)
# =============================================================================

# Group 1: PC/PA (80:20)
group1_data = {
    'name': 'Group 1: PC/PA (80:20)',
    'short_name': '80:20 PC/PA',
    'systems': ['DOPC/DOPA', 'DOPC/POPA', 'DOPC/SOPA',
                'POPC/DOPA', 'POPC/POPA', 'POPC/SOPA',
                'SOPC/DOPA', 'SOPC/POPA', 'SOPC/SOPA'],
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
    'secondary_type': ['DOPA', 'POPA', 'SOPA', 'DOPA', 'POPA', 'SOPA', 'DOPA', 'POPA', 'SOPA'],
    'PC_PC': [9088, 8888, 9065, 13082, 12863, 13180, 11859, 11673, 12003],
    'PC_Secondary': [-17639, -16906, -17421, -17205, -16623, -17129, -17339, -16673, -17263],
    'Secondary_Secondary': [-38791, -38509, -38620, -38406, -38673, -38709, -39136, -39183, -38821],
    'Total': [-47342, -46527, -46976, -42529, -42433, -42658, -44616, -44183, -44081],
}

# Group 2: PC/PG (65:35)
group2_data = {
    'name': 'Group 2: PC/PG (65:35)',
    'short_name': '65:35 PC/PG',
    'systems': ['DOPC/DOPG', 'DOPC/POPG', 'DOPC/SOPG',
                'POPC/DOPG', 'POPC/POPG', 'POPC/SOPG',
                'SOPC/DOPG', 'SOPC/POPG', 'SOPC/SOPG'],
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
    'secondary_type': ['DOPG', 'POPG', 'SOPG', 'DOPG', 'POPG', 'SOPG', 'DOPG', 'POPG', 'SOPG'],
    'PC_PC': [17266, 16890, 17440, 20024, 19839, 20082, 19553, 18846, 19798],
    'PC_Secondary': [-25981, -25007, -25898, -25410, -24762, -24991, -26098, -24432, -25987],
    'Secondary_Secondary': [-110579, -109618, -109475, -110612, -109393, -109515, -110392, -109649, -109210],
    'Total': [-119294, -117735, -117933, -115997, -114315, -114424, -116937, -115235, -115400],
}

# Group 3: PC/PG (95:5)
group3_data = {
    'name': 'Group 3: PC/PG (95:5)',
    'short_name': '95:5 PC/PG',
    'systems': ['DOPC/DOPG', 'DOPC/POPG', 'DOPC/SOPG',
                'POPC/DOPG', 'POPC/POPG', 'POPC/SOPG',
                'SOPC/DOPG', 'SOPC/POPG', 'SOPC/SOPG'],
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
    'secondary_type': ['DOPG', 'POPG', 'SOPG', 'DOPG', 'POPG', 'SOPG', 'DOPG', 'POPG', 'SOPG'],
    'PC_PC': [8753, 8874, 8808, 13400, 13454, 13575, 12233, 12012, 12120],
    'PC_Secondary': [-5121, -5140, -5257, -4975, -5004, -5130, -5170, -5074, -5219],
    'Secondary_Secondary': [-14858, -14671, -14633, -14875, -14691, -14666, -14810, -14679, -14647],
    'Total': [-11226, -10938, -11082, -6451, -6241, -6222, -7746, -7740, -7746],
}

ALL_GROUPS = [group1_data, group2_data, group3_data]


# =============================================================================
# FUNCTION: Create Single Interaction Plot
# =============================================================================

def create_interaction_plot(ax, data_key, title, scale_factor=10000):
    """
    Create a grouped bar plot for a specific interaction type.
    Systems grouped by PC type within each group.
    """
    bar_width = 0.08
    
    # Define x positions for each group
    group_positions = [0, 1.2, 2.4]  # Spacing between groups
    
    for group_idx, (group, group_pos) in enumerate(zip(ALL_GROUPS, group_positions)):
        # For each PC type
        for pc_idx, pc in enumerate(['DOPC', 'POPC', 'SOPC']):
            # Get indices for this PC type
            indices = [i for i, p in enumerate(group['PC_type']) if p == pc]
            
            # Calculate x positions for this PC type's bars
            pc_offset = pc_idx * 0.28  # Offset for each PC type
            
            for i, idx in enumerate(indices):
                x_pos = group_pos + pc_offset + i * bar_width
                value = group[data_key][idx] / scale_factor
                
                bar = ax.bar(x_pos, value, bar_width * 0.9, 
                            color=PC_COLORS[pc], edgecolor='black', linewidth=0.7)
    
    # Customize axes
    ax.set_ylabel('Interaction Energy (×10$^4$ kJ/mol)', fontsize=12)
    ax.set_title(title, fontweight='bold', fontsize=14)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Set x-ticks at group centers
    group_centers = [pos + 0.35 for pos in group_positions]
    ax.set_xticks(group_centers)
    ax.set_xticklabels([g['short_name'] for g in ALL_GROUPS], fontsize=11)
    
    # Add legend
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=PC_COLORS[pc], 
                                      edgecolor='black', label=pc)
                       for pc in ['DOPC', 'POPC', 'SOPC']]
    ax.legend(handles=legend_elements, loc='best', title='PC Lipid', 
              framealpha=0.9, fontsize=10)
    
    return ax


# =============================================================================
# FUNCTION: Create 4-Panel Figure
# =============================================================================

def create_four_panel_figure():
    """
    Create a 2x2 figure with all four interaction types.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # (A) Secondary-Secondary interactions
    create_interaction_plot(
        axes[0, 0], 
        'Secondary_Secondary',
        '(A) Secondary-Secondary Lipid Total Interaction Energies'
    )
    
    # (B) PC-PC interactions
    create_interaction_plot(
        axes[0, 1], 
        'PC_PC',
        '(B) PC-PC Lipd Total Interaction Energies'
    )
    
    # (C) PC-Secondary interactions
    create_interaction_plot(
        axes[1, 0], 
        'PC_Secondary',
        '(C) PC-Secondary Lipid Total Interaction Energies'
    )
    
    # (D) Total interactions
    create_interaction_plot(
        axes[1, 1], 
        'Total',
        '(D) Total Lipid-Lipid Interaction Energies'
    )
    
    plt.tight_layout()
    plt.savefig('Decomposed_Interaction_Energies_4Panel.png', dpi=300, 
                facecolor='white', edgecolor='none')
    plt.savefig('Decomposed_Interaction_Energies_4Panel.pdf', 
                facecolor='white', edgecolor='none')
    print("Saved: Decomposed_Interaction_Energies_4Panel.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Create Individual Plots (Larger, More Detailed)
# =============================================================================

def create_individual_plot(data_key, title, filename):
    """
    Create a single detailed plot for one interaction type.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    bar_width = 0.07
    scale_factor = 10000
    
    # Define x positions for each group with more spacing
    group_positions = [0, 1.5, 3.0]
    
    all_bars = []
    all_labels = []
    
    for group_idx, (group, group_pos) in enumerate(zip(ALL_GROUPS, group_positions)):
        # For each PC type
        for pc_idx, pc in enumerate(['DOPC', 'POPC', 'SOPC']):
            # Get indices for this PC type
            indices = [i for i, p in enumerate(group['PC_type']) if p == pc]
            
            # Calculate x positions
            pc_offset = pc_idx * 0.32
            
            for i, idx in enumerate(indices):
                x_pos = group_pos + pc_offset + i * bar_width
                value = group[data_key][idx] / scale_factor
                
                bar = ax.bar(x_pos, value, bar_width * 0.85, 
                            color=PC_COLORS[pc], edgecolor='black', linewidth=0.8)
                
                # Add value label
                if value >= 0:
                    va = 'bottom'
                    y_pos = value + 0.1
                else:
                    va = 'top'
                    y_pos = value - 0.1
                
                ax.text(x_pos, y_pos, f'{value:.1f}', ha='center', va=va, 
                       fontsize=7, rotation=90)
    
    # Customize axes
    ax.set_ylabel('Interaction Energy (×10$^4$ kJ/mol)', fontsize=13)
    ax.set_title(title, fontweight='bold', fontsize=15)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Set x-ticks
    group_centers = [pos + 0.42 for pos in group_positions]
    ax.set_xticks(group_centers)
    ax.set_xticklabels([g['short_name'] for g in ALL_GROUPS], fontsize=12)
    
    # Add vertical separators between groups
    for i in range(len(group_positions) - 1):
        sep_x = (group_positions[i] + 0.95 + group_positions[i+1]) / 2
        ax.axvline(x=sep_x, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    
    # Legend
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=PC_COLORS[pc], 
                                      edgecolor='black', label=pc)
                       for pc in ['DOPC', 'POPC', 'SOPC']]
    ax.legend(handles=legend_elements, loc='best', title='PC Lipid', 
              framealpha=0.9, fontsize=11)
    
    plt.tight_layout()
    plt.savefig(f'{filename}.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig(f'{filename}.pdf', facecolor='white', edgecolor='none')
    print(f"Saved: {filename}.png/pdf")
    plt.show()


def create_all_individual_plots():
    """
    Create all four individual plots.
    """
    plots = [
        ('Secondary_Secondary', '(A) Secondary-Secondary Interactions', 
         'Decomposed_Secondary_Secondary'),
        ('PC_PC', '(B) PC-PC Interactions', 'Decomposed_PC_PC'),
        ('PC_Secondary', '(C) PC-Secondary Interactions', 'Decomposed_PC_Secondary'),
        ('Total', '(D) Total Lipid-Lipid Interactions', 'Decomposed_Total'),
    ]
    
    for data_key, title, filename in plots:
        create_individual_plot(data_key, title, filename)


# =============================================================================
# FUNCTION: Create Stacked Bar Chart (Alternative Visualization)
# =============================================================================

def create_stacked_comparison():
    """
    Create a stacked bar chart showing the contribution of each interaction type.
    """
    fig, ax = plt.subplots(figsize=(16, 8))
    
    bar_width = 0.08
    scale_factor = 10000
    
    # Collect all systems
    all_systems = []
    all_pc_pc = []
    all_pc_sec = []
    all_sec_sec = []
    all_colors = []
    group_boundaries = []
    
    x_pos = 0
    for group_idx, group in enumerate(ALL_GROUPS):
        group_start = len(all_systems)  # index of the group's first system
        for i in range(len(group['systems'])):
            all_systems.append(group['systems'][i])
            all_pc_pc.append(group['PC_PC'][i] / scale_factor)
            all_pc_sec.append(group['PC_Secondary'][i] / scale_factor)
            all_sec_sec.append(group['Secondary_Secondary'][i] / scale_factor)
            all_colors.append(PC_COLORS[group['PC_type'][i]])
            x_pos += 1
        group_boundaries.append((group_start, len(all_systems) - 1, group['short_name']))
        x_pos += 0.5  # Gap between groups
    
    # Create x positions
    x = np.arange(len(all_systems))
    x_adjusted = []
    current = 0
    for i, sys in enumerate(all_systems):
        x_adjusted.append(current)
        current += 1
        # Add gap after each group
        if i == 8 or i == 17:
            current += 0.5
    x_adjusted = np.array(x_adjusted)
    
    # Plot stacked bars
    # Bottom layer: Secondary-Secondary (most negative)
    bars1 = ax.bar(x_adjusted, all_sec_sec, width=0.7, label='Secondary-Secondary', 
                   color='#e74c3c', edgecolor='black', linewidth=0.5)
    
    # Middle layer: PC-Secondary
    bars2 = ax.bar(x_adjusted, all_pc_sec, width=0.7, bottom=all_sec_sec,
                   label='PC-Secondary', color='#f39c12', edgecolor='black', linewidth=0.5)
    
    # Top layer: PC-PC (positive, starts from 0)
    bars3 = ax.bar(x_adjusted, all_pc_pc, width=0.7, label='PC-PC', 
                   color='#3498db', edgecolor='black', linewidth=0.5)
    
    # Add total line
    totals = [s + p + pc for s, p, pc in zip(all_sec_sec, all_pc_sec, all_pc_pc)]
    ax.plot(x_adjusted, totals, 'ko-', markersize=5, linewidth=1.5, label='Total')
    
    # Customize
    ax.set_ylabel('Interaction Energy (×10$^4$ kJ/mol)', fontsize=13)
    ax.set_title('Decomposed Interaction Energies: All 27 Systems', fontweight='bold', fontsize=15)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # X-axis labels
    ax.set_xticks(x_adjusted)
    ax.set_xticklabels(all_systems, rotation=45, ha='right', fontsize=8)
    
    # Add group labels
    for start, end, label in group_boundaries:
        mid = (x_adjusted[start] + x_adjusted[end]) / 2
        ax.text(mid, ax.get_ylim()[1] * 0.95, label, ha='center', va='top', 
                fontsize=11, fontweight='bold', style='italic')
        if end < len(all_systems) - 1:
            ax.axvline(x=(x_adjusted[end] + x_adjusted[min(end+1, len(x_adjusted)-1)]) / 2 + 0.25, 
                      color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    ax.legend(loc='lower left', fontsize=10, framealpha=0.9)
    
    plt.tight_layout()
    plt.savefig('Decomposed_Stacked_All_Systems.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('Decomposed_Stacked_All_Systems.pdf', facecolor='white', edgecolor='none')
    print("Saved: Decomposed_Stacked_All_Systems.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Print Summary
# =============================================================================

def print_summary():
    """
    Print summary of decomposed interaction energies.
    """
    print("\n" + "=" * 85)
    print("DECOMPOSED TOTAL INTERACTION ENERGIES SUMMARY (×10⁴ kJ/mol)")
    print("=" * 85)
    
    for group in ALL_GROUPS:
        print(f"\n{group['name']}")
        print("-" * 75)
        print(f"{'System':<15} {'Sec-Sec':>12} {'PC-Sec':>12} {'PC-PC':>12} {'Total':>12}")
        print("-" * 75)
        
        for i, sys in enumerate(group['systems']):
            ss = group['Secondary_Secondary'][i] / 10000
            ps = group['PC_Secondary'][i] / 10000
            pp = group['PC_PC'][i] / 10000
            total = group['Total'][i] / 10000
            print(f"{sys:<15} {ss:>+12.2f} {ps:>+12.2f} {pp:>+12.2f} {total:>+12.2f}")
        
        # Averages
        print("-" * 75)
        ss_avg = np.mean(group['Secondary_Secondary']) / 10000
        ps_avg = np.mean(group['PC_Secondary']) / 10000
        pp_avg = np.mean(group['PC_PC']) / 10000
        total_avg = np.mean(group['Total']) / 10000
        print(f"{'Average':<15} {ss_avg:>+12.2f} {ps_avg:>+12.2f} {pp_avg:>+12.2f} {total_avg:>+12.2f}")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 85)
    print("Decomposed Total Interaction Energies in Mixed Lipid Bilayers")
    print("=" * 85)
    
    # Print summary
    print_summary()
    
    # Generate plots
    print("\n[1/3] Generating 4-panel figure...")
    create_four_panel_figure()
    
    print("\n[2/3] Generating individual plots...")
    create_all_individual_plots()
    
    print("\n[3/3] Generating stacked comparison plot...")
    create_stacked_comparison()
    
    print("\n" + "=" * 85)
    print("All plots generated successfully!")
    print("=" * 85)
