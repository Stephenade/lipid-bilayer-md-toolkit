#!/usr/bin/env python3
"""
Partial Interaction Energy Comparison in Mixed Lipid Bilayers
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
    'axes.titlesize': 13,
    'xtick.labelsize': 9,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
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
    # PC-PC interactions
    'PC_PC': [9088, 8888, 9065, 13082, 12863, 13180, 11859, 11673, 12003],
    # PC-Secondary interactions
    'PC_Secondary': [-17639, -16906, -17421, -17205, -16623, -17129, -17339, -16673, -17263],
    # Secondary-Secondary interactions
    'Secondary_Secondary': [-38791, -38509, -38620, -38406, -38673, -38709, -39136, -39183, -38821],
    # Total
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
    # PC-PC interactions
    'PC_PC': [17266, 16890, 17440, 20024, 19839, 20082, 19553, 18846, 19798],
    # PC-Secondary interactions
    'PC_Secondary': [-25981, -25007, -25898, -25410, -24762, -24991, -26098, -24432, -25987],
    # Secondary-Secondary interactions
    'Secondary_Secondary': [-110579, -109618, -109475, -110612, -109393, -109515, -110392, -109649, -109210],
    # Total
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
    # PC-PC interactions
    'PC_PC': [8753, 8874, 8808, 13400, 13454, 13575, 12233, 12012, 12120],
    # PC-Secondary interactions
    'PC_Secondary': [-5121, -5140, -5257, -4975, -5004, -5130, -5170, -5074, -5219],
    # Secondary-Secondary interactions
    'Secondary_Secondary': [-14858, -14671, -14633, -14875, -14691, -14666, -14810, -14679, -14647],
    # Total
    'Total': [-11226, -10938, -11082, -6451, -6241, -6222, -7746, -7740, -7746],
}

ALL_GROUPS = [group1_data, group2_data, group3_data]


# =============================================================================
# FUNCTION: Create Grouped Bar Plot
# =============================================================================

def create_grouped_bar_plot(ax, data_key, title, ylabel, scale_factor=10000):
    """
    Create a grouped bar plot for a specific interaction type.
    Groups are arranged by PC lipid type (DOPC, POPC, SOPC) within each group.
    """
    bar_width = 0.25
    group_spacing = 1.0
    
    # Positions for each group
    n_systems_per_pc = 3  # 3 secondary variants per PC type
    n_pc_types = 3  # DOPC, POPC, SOPC
    
    all_positions = []
    all_values = []
    all_colors = []
    all_labels = []
    
    x_offset = 0
    group_centers = []
    group_labels = []
    
    for group_idx, group in enumerate(ALL_GROUPS):
        group_start = x_offset
        
        for pc_idx, pc in enumerate(['DOPC', 'POPC', 'SOPC']):
            # Get indices for this PC type
            indices = [i for i, p in enumerate(group['PC_type']) if p == pc]
            
            for i, idx in enumerate(indices):
                pos = x_offset + i * bar_width
                all_positions.append(pos)
                all_values.append(group[data_key][idx] / scale_factor)
                all_colors.append(PC_COLORS[pc])
                all_labels.append(group['secondary_type'][idx].replace('DO', 'DO-').replace('PO', 'PO-').replace('SO', 'SO-'))
            
            x_offset += n_systems_per_pc * bar_width + 0.15  # Small gap between PC types
        
        group_end = x_offset - 0.15
        group_centers.append((group_start + group_end) / 2)
        group_labels.append(group['short_name'])
        x_offset += group_spacing  # Larger gap between groups
    
    # Create bars
    bars = ax.bar(all_positions, all_values, bar_width * 0.9, 
                  color=all_colors, edgecolor='black', linewidth=0.8)
    
    # Add value labels on bars
    for bar, val in zip(bars, all_values):
        height = bar.get_height()
        if height >= 0:
            va = 'bottom'
            y_pos = height + 0.05
        else:
            va = 'top'
            y_pos = height - 0.05
        ax.text(bar.get_x() + bar.get_width()/2., y_pos,
                f'{val:.1f}', ha='center', va=va, fontsize=7, rotation=90)
    
    # Customize axes
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold', fontsize=13)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Set x-ticks at group centers
    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_labels, fontsize=11)
    
    # Add legend for PC types
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=PC_COLORS[pc], 
                                      edgecolor='black', label=pc)
                       for pc in ['DOPC', 'POPC', 'SOPC']]
    ax.legend(handles=legend_elements, loc='best', title='PC Lipid', 
              framealpha=0.9, fontsize=9)
    
    return ax


# =============================================================================
# FUNCTION: Create 4-Panel Figure
# =============================================================================

def create_four_panel_figure():
    """
    Create a 2x2 figure with all four interaction types.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # (A) Secondary-Secondary interactions
    create_grouped_bar_plot(
        axes[0, 0], 
        'Secondary_Secondary',
        '(A) Secondary-Secondary Interactions (PA-PA or PG-PG)',
        'Interaction Energy (×10$^4$ kJ/mol)'
    )
    
    # (B) PC-PC interactions
    create_grouped_bar_plot(
        axes[0, 1], 
        'PC_PC',
        '(B) PC-PC Interactions',
        'Interaction Energy (×10$^4$ kJ/mol)'
    )
    
    # (C) PC-Secondary interactions
    create_grouped_bar_plot(
        axes[1, 0], 
        'PC_Secondary',
        '(C) PC-Secondary Interactions (PC-PA or PC-PG)',
        'Interaction Energy (×10$^4$ kJ/mol)'
    )
    
    # (D) Total interactions
    create_grouped_bar_plot(
        axes[1, 1], 
        'Total',
        '(D) Total Lipid-Lipid Interactions',
        'Interaction Energy (×10$^4$ kJ/mol)'
    )
    
    plt.tight_layout()
    plt.savefig('Partial_Interaction_Energies_4Panel.png', dpi=300, 
                facecolor='white', edgecolor='none')
    plt.savefig('Partial_Interaction_Energies_4Panel.pdf', 
                facecolor='white', edgecolor='none')
    print("Saved: Partial_Interaction_Energies_4Panel.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Create Individual Plots
# =============================================================================

def create_individual_plots():
    """
    Create separate figures for each interaction type.
    """
    plot_configs = [
        ('Secondary_Secondary', '(A) Secondary-Secondary Interactions (PA-PA or PG-PG)', 
         'Secondary_Secondary_Interactions'),
        ('PC_PC', '(B) PC-PC Interactions', 'PC_PC_Interactions'),
        ('PC_Secondary', '(C) PC-Secondary Interactions (PC-PA or PC-PG)', 
         'PC_Secondary_Interactions'),
        ('Total', '(D) Total Lipid-Lipid Interactions', 'Total_Interactions'),
    ]
    
    for data_key, title, filename in plot_configs:
        fig, ax = plt.subplots(figsize=(14, 6))
        create_grouped_bar_plot(ax, data_key, title, 
                                'Interaction Energy (×10$^4$ kJ/mol)')
        plt.tight_layout()
        plt.savefig(f'{filename}.png', dpi=300, facecolor='white', edgecolor='none')
        plt.savefig(f'{filename}.pdf', facecolor='white', edgecolor='none')
        print(f"Saved: {filename}.png/pdf")
        plt.show()


# =============================================================================
# FUNCTION: Create Detailed Comparison Plot
# =============================================================================

def create_detailed_comparison():
    """
    Create a more detailed comparison showing all systems side by side.
    """
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    
    interaction_types = [
        ('Secondary_Secondary', '(A) Secondary-Secondary Interactions'),
        ('PC_PC', '(B) PC-PC Interactions'),
        ('PC_Secondary', '(C) PC-Secondary Interactions'),
        ('Total', '(D) Total Lipid-Lipid Interactions'),
    ]
    
    for ax, (data_key, title) in zip(axes.flat, interaction_types):
        # Collect all data
        all_systems = []
        all_values = []
        all_colors = []
        group_boundaries = []
        
        current_pos = 0
        for group in ALL_GROUPS:
            group_start = current_pos
            for i, (sys, pc, val) in enumerate(zip(group['systems'], 
                                                    group['PC_type'], 
                                                    group[data_key])):
                all_systems.append(sys)
                all_values.append(val / 10000)  # Scale
                all_colors.append(PC_COLORS[pc])
                current_pos += 1
            group_boundaries.append((group_start, current_pos - 1, group['short_name']))
        
        # Create bars
        x = np.arange(len(all_systems))
        bars = ax.bar(x, all_values, color=all_colors, edgecolor='black', linewidth=0.6)
        
        # Add group separators and labels
        for start, end, label in group_boundaries:
            mid = (start + end) / 2
            # Add vertical line between groups (except at the end)
            if end < len(all_systems) - 1:
                ax.axvline(x=end + 0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)
            # Add group label at top
            y_max = max(all_values) if max(all_values) > 0 else 0
            y_min = min(all_values) if min(all_values) < 0 else 0
            y_pos = y_max + abs(y_max - y_min) * 0.05
            ax.text(mid, y_pos, label, ha='center', va='bottom', fontsize=10, 
                    fontweight='bold', style='italic')
        
        # Customize
        ax.set_xticks(x)
        ax.set_xticklabels(all_systems, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Interaction Energy (×10$^4$ kJ/mol)')
        ax.set_title(title, fontweight='bold', fontsize=12)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        
        # Legend
        legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=PC_COLORS[pc], 
                                          edgecolor='black', label=pc)
                           for pc in ['DOPC', 'POPC', 'SOPC']]
        ax.legend(handles=legend_elements, loc='best', title='PC Lipid', fontsize=8)
    
    plt.tight_layout()
    plt.savefig('Partial_Interactions_Detailed.png', dpi=300, 
                facecolor='white', edgecolor='none')
    plt.savefig('Partial_Interactions_Detailed.pdf', 
                facecolor='white', edgecolor='none')
    print("Saved: Partial_Interactions_Detailed.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Print Summary Statistics
# =============================================================================

def print_summary_statistics():
    """
    Print summary statistics for all interaction types.
    """
    print("\n" + "=" * 80)
    print("PARTIAL INTERACTION ENERGY SUMMARY (All values in kJ/mol)")
    print("=" * 80)
    
    for group in ALL_GROUPS:
        print(f"\n{group['name']}")
        print("-" * 70)
        print(f"{'System':<15} {'PC-PC':>12} {'PC-Secondary':>14} {'Sec-Sec':>12} {'Total':>12}")
        print("-" * 70)
        
        for i, sys in enumerate(group['systems']):
            print(f"{sys:<15} {group['PC_PC'][i]:>+12,} {group['PC_Secondary'][i]:>+14,} "
                  f"{group['Secondary_Secondary'][i]:>+12,} {group['Total'][i]:>+12,}")
        
        # Averages
        print("-" * 70)
        print(f"{'Average':<15} {np.mean(group['PC_PC']):>+12,.0f} "
              f"{np.mean(group['PC_Secondary']):>+14,.0f} "
              f"{np.mean(group['Secondary_Secondary']):>+12,.0f} "
              f"{np.mean(group['Total']):>+12,.0f}")
    
    print("\n" + "=" * 80)
    print("KEY OBSERVATIONS:")
    print("=" * 80)
    print("""
1. SECONDARY-SECONDARY INTERACTIONS:
   - Strongest in 60:40 PG systems (~-110,000 kJ/mol)
   - Moderate in 90:10 PA systems (~-38,500 kJ/mol)
   - Weakest in 80:20 PG systems (~-14,700 kJ/mol)
   - Scales with secondary lipid concentration

2. PC-PC INTERACTIONS:
   - Always repulsive (positive values)
   - POPC-POPC > SOPC-SOPC > DOPC-DOPC
   - Higher in 60:40 systems (fewer PC lipids, more concentrated)

3. PC-SECONDARY INTERACTIONS:
   - Always attractive (negative values)
   - Strongest in 60:40 PG systems (~-25,000 kJ/mol)
   - Moderate in 90:10 PA systems (~-17,000 kJ/mol)
   - Weakest in 80:20 PG systems (~-5,000 kJ/mol)

4. TOTAL INTERACTIONS:
   - All systems have net attractive interactions
   - 60:40 PG systems most cohesive (~-116,000 kJ/mol)
   - 90:10 PA systems moderately cohesive (~-44,000 kJ/mol)
   - 80:20 PG systems least cohesive (~-8,000 kJ/mol)
    """)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("Partial Interaction Energy Comparison in Mixed Lipid Bilayers")
    print("=" * 80)
    
    # Print statistics
    print_summary_statistics()
    
    # Generate plots
    print("\n[1/3] Generating 4-panel figure...")
    create_four_panel_figure()
    
    print("\n[2/3] Generating individual plots...")
    create_individual_plots()
    
    print("\n[3/3] Generating detailed comparison plot...")
    create_detailed_comparison()
    
    print("\n" + "=" * 80)
    print("All plots generated successfully!")
    print("=" * 80)
