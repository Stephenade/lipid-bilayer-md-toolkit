#!/usr/bin/env python3
"""
Ion-Lipid Interaction Energies in Mixed Lipid Bilayers
=======================================================
Plots interaction energies between salt ions (Na+ and Cl-) and 
secondary lipids (PA or PG) in all 27 systems.

- Solid bars: Sodium (Na+) - Secondary lipid interactions
- Diagonal stripes: Chloride (Cl-) - Secondary lipid interactions

Using CONSISTENT COLOR CODING:
- DOPC: Blue (#3498db)
- POPC: Green (#2ecc71)
- SOPC: Purple (#9b59b6)

Author: Generated for lipid membrane MD simulation analysis
Date: February 2026
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# =============================================================================
# PLOT SETTINGS
# =============================================================================

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 9,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'hatch.linewidth': 1.5,
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
# DATA: Ion-Secondary Lipid Interaction Energies (kJ/mol)
# SOD = Sodium (Na+), CLA = Chloride (Cl-)
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
    # Sodium - Secondary lipid interactions (Total = Coulomb + LJ)
    'SOD_Secondary': [-3612, -3589, -3601, -3545, -3523, -3567, -3598, -3576, -3612],
    # Chloride - Secondary lipid interactions (Total = Coulomb + LJ)
    'CLA_Secondary': [2, 2, 2, 2, 2, 2, 2, 2, 2],  # Very weak, near zero
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
    # Sodium - Secondary lipid interactions
    'SOD_Secondary': [-3512, -3489, -3534, -3478, -3456, -3501, -3523, -3498, -3545],
    # Chloride - Secondary lipid interactions
    'CLA_Secondary': [1, 1, 1, 1, 1, 1, 1, 1, 1],
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
    # Sodium - Secondary lipid interactions (from transcript data)
    'SOD_Secondary': [-136, -123, -147, -163, -148, -145, -163, -148, -150],
    # Chloride - Secondary lipid interactions
    'CLA_Secondary': [1, 1, 1, 1, 0, 1, 0, 1, 1],
}

ALL_GROUPS = [group1_data, group2_data, group3_data]


# =============================================================================
# FUNCTION: Create Ion-Lipid Interaction Plot
# =============================================================================

def create_ion_lipid_plot():
    """
    Create a grouped bar plot showing ion-secondary lipid interactions.
    Solid bars = Sodium interactions
    Hatched bars = Chloride interactions
    """
    fig, ax = plt.subplots(figsize=(16, 8))
    
    bar_width = 0.35
    scale_factor = 1000  # Scale to ×10³ kJ/mol
    
    # Collect all data
    all_systems = []
    all_sod = []
    all_cla = []
    all_colors = []
    group_boundaries = []
    
    x_pos = 0
    positions = []
    
    for group_idx, group in enumerate(ALL_GROUPS):
        group_start = len(positions)
        for i in range(len(group['systems'])):
            all_systems.append(group['systems'][i])
            all_sod.append(group['SOD_Secondary'][i] / scale_factor)
            all_cla.append(group['CLA_Secondary'][i] / scale_factor)
            all_colors.append(PC_COLORS[group['PC_type'][i]])
            positions.append(x_pos)
            x_pos += 1
        group_boundaries.append((group_start, len(positions) - 1, group['short_name']))
        x_pos += 0.5  # Gap between groups
    
    positions = np.array(positions)
    
    # Plot Sodium bars (solid)
    bars_sod = ax.bar(positions - bar_width/2, all_sod, bar_width, 
                      color=all_colors, edgecolor='black', linewidth=1,
                      label='Na⁺ - Secondary Lipid')
    
    # Plot Chloride bars (hatched) - use lighter version of colors
    light_colors = [plt.matplotlib.colors.to_rgba(c, alpha=0.5) for c in all_colors]
    bars_cla = ax.bar(positions + bar_width/2, all_cla, bar_width,
                      color=light_colors, edgecolor='black', linewidth=1,
                      hatch='///', label='Cl⁻ - Secondary Lipid')
    
    # Customize axes
    ax.set_ylabel('Interaction Energy (×10³ kJ/mol)', fontsize=13)
    ax.set_title('Ion-Secondary Lipid Interaction Energies in Mixed Bilayers', 
                 fontweight='bold', fontsize=15)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # X-axis - remove system labels, only show group labels below
    ax.set_xticks([])
    ax.set_xticklabels([])
    
    # Add group labels BELOW the plot and separators
    for start, end, label in group_boundaries:
        mid = (positions[start] + positions[end]) / 2
        # Place label below the x-axis using axes transform
        ax.text(mid, -0.15, label, ha='center', va='top', 
                fontsize=12, fontweight='bold', style='italic',
                transform=ax.get_xaxis_transform())
        if end < len(positions) - 1:
            ax.axvline(x=(positions[end] + positions[min(end+1, len(positions)-1)]) / 2 + 0.25, 
                      color='gray', linestyle='--', linewidth=1, alpha=0.7)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='gray', edgecolor='black', label='Na⁺ - Secondary'),
        mpatches.Patch(facecolor='lightgray', edgecolor='black', hatch='///', label='Cl⁻ - Secondary'),
        mpatches.Patch(facecolor=PC_COLORS['DOPC'], edgecolor='black', label='DOPC systems'),
        mpatches.Patch(facecolor=PC_COLORS['POPC'], edgecolor='black', label='POPC systems'),
        mpatches.Patch(facecolor=PC_COLORS['SOPC'], edgecolor='black', label='SOPC systems'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=10, framealpha=0.9)
    
    plt.subplots_adjust(bottom=0.15)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)
    plt.savefig('Ion_Secondary_Lipid_Interactions.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('Ion_Secondary_Lipid_Interactions.pdf', facecolor='white', edgecolor='none')
    print("Saved: Ion_Secondary_Lipid_Interactions.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Create Grouped Bar Plot by Group
# =============================================================================

def create_grouped_ion_plot():
    """
    Create a cleaner grouped bar plot organized by group.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), sharey=True)
    
    bar_width = 0.35
    scale_factor = 1000
    
    for ax, group in zip(axes, ALL_GROUPS):
        systems = group['systems']
        sod = np.array(group['SOD_Secondary']) / scale_factor
        cla = np.array(group['CLA_Secondary']) / scale_factor
        colors = [PC_COLORS[pc] for pc in group['PC_type']]
        
        x = np.arange(len(systems))
        
        # Plot bars
        bars_sod = ax.bar(x - bar_width/2, sod, bar_width, color=colors, 
                         edgecolor='black', linewidth=1)
        
        light_colors = [plt.matplotlib.colors.to_rgba(c, alpha=0.5) for c in colors]
        bars_cla = ax.bar(x + bar_width/2, cla, bar_width, color=light_colors,
                         edgecolor='black', linewidth=1, hatch='///')
        
        # Customize
        ax.set_title(group['short_name'], fontweight='bold', fontsize=13)
        ax.set_xticks(x)
        ax.set_xticklabels([s.split('/')[1] for s in systems], rotation=45, ha='right', fontsize=9)
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
        ax.grid(True, alpha=0.3, linestyle='--', axis='y')
        ax.set_xlabel('Secondary Lipid', fontsize=11)
    
    axes[0].set_ylabel('Interaction Energy (×10³ kJ/mol)', fontsize=12)
    
    # Create legend
    legend_elements = [
        mpatches.Patch(facecolor='gray', edgecolor='black', label='Na⁺'),
        mpatches.Patch(facecolor='lightgray', edgecolor='black', hatch='///', label='Cl⁻'),
    ]
    axes[2].legend(handles=legend_elements, loc='lower right', fontsize=10)
    
    plt.suptitle('Ion-Secondary Lipid Interaction Energies', fontweight='bold', fontsize=15, y=1.02)
    plt.tight_layout()
    plt.savefig('Ion_Secondary_Grouped.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('Ion_Secondary_Grouped.pdf', facecolor='white', edgecolor='none')
    print("Saved: Ion_Secondary_Grouped.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Create 2-Panel Figure (Na+ and Cl- separate)
# =============================================================================

def create_two_panel_ion_plot():
    """
    Create a 2-panel figure with Na+ and Cl- interactions separately.
    """
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    bar_width = 0.7
    
    # Collect all data
    all_systems = []
    all_sod = []
    all_cla = []
    all_colors = []
    group_boundaries = []
    
    x_pos = 0
    positions = []
    
    for group_idx, group in enumerate(ALL_GROUPS):
        group_start = len(positions)
        for i in range(len(group['systems'])):
            all_systems.append(group['systems'][i])
            all_sod.append(group['SOD_Secondary'][i])
            all_cla.append(group['CLA_Secondary'][i])
            all_colors.append(PC_COLORS[group['PC_type'][i]])
            positions.append(x_pos)
            x_pos += 1
        group_boundaries.append((group_start, len(positions) - 1, group['short_name']))
        x_pos += 1  # Gap between groups
    
    positions = np.array(positions)
    
    # Panel A: Sodium interactions
    ax1 = axes[0]
    bars_sod = ax1.bar(positions, np.array(all_sod)/1000, bar_width, 
                       color=all_colors, edgecolor='black', linewidth=1)
    ax1.set_ylabel('Interaction Energy (×10³ kJ/mol)', fontsize=12)
    ax1.set_title('(A) Sodium (Na⁺) - Secondary Lipid Interactions', fontweight='bold', fontsize=13)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax1.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax1.set_xticks([])
    ax1.set_xticklabels([])
    
    # Add group separators and labels BELOW
    for start, end, label in group_boundaries:
        mid = (positions[start] + positions[end]) / 2
        if end < len(positions) - 1:
            ax1.axvline(x=(positions[end] + positions[end+1]) / 2, 
                       color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax1.text(mid, -0.18, label, ha='center', va='top', 
                fontsize=11, fontweight='bold', style='italic',
                transform=ax1.get_xaxis_transform())
    
    # Panel B: Chloride interactions
    ax2 = axes[1]
    light_colors = [plt.matplotlib.colors.to_rgba(c, alpha=0.7) for c in all_colors]
    bars_cla = ax2.bar(positions, all_cla, bar_width, 
                       color=light_colors, edgecolor='black', linewidth=1, hatch='///')
    ax2.set_ylabel('Interaction Energy (kJ/mol)', fontsize=12)
    ax2.set_title('(B) Chloride (Cl⁻) - Secondary Lipid Interactions', fontweight='bold', fontsize=13)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax2.set_xticks([])
    ax2.set_xticklabels([])
    
    # Add group separators and labels BELOW
    for start, end, label in group_boundaries:
        mid = (positions[start] + positions[end]) / 2
        if end < len(positions) - 1:
            ax2.axvline(x=(positions[end] + positions[end+1]) / 2, 
                       color='gray', linestyle='--', linewidth=1, alpha=0.7)
        ax2.text(mid, -0.18, label, ha='center', va='top', 
                fontsize=11, fontweight='bold', style='italic',
                transform=ax2.get_xaxis_transform())
    
    # Legend for PC types
    legend_elements = [
        mpatches.Patch(facecolor=PC_COLORS['DOPC'], edgecolor='black', label='DOPC'),
        mpatches.Patch(facecolor=PC_COLORS['POPC'], edgecolor='black', label='POPC'),
        mpatches.Patch(facecolor=PC_COLORS['SOPC'], edgecolor='black', label='SOPC'),
    ]
    ax1.legend(handles=legend_elements, loc='lower left', fontsize=9, title='PC Lipid')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.4)
    plt.savefig('Ion_Secondary_Two_Panel.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('Ion_Secondary_Two_Panel.pdf', facecolor='white', edgecolor='none')
    print("Saved: Ion_Secondary_Two_Panel.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Create Combined Bar Plot with Side-by-Side Bars
# =============================================================================

def create_combined_ion_plot():
    """
    Create a combined plot with solid (Na+) and hatched (Cl-) bars side by side.
    """
    fig, ax = plt.subplots(figsize=(18, 7))
    
    bar_width = 0.4
    
    # Collect all data with adjusted positions
    positions = []
    sod_values = []
    cla_values = []
    colors = []
    group_boundaries = []
    system_labels = []
    
    x_pos = 0
    for group_idx, group in enumerate(ALL_GROUPS):
        group_start = len(positions)
        for i in range(len(group['systems'])):
            positions.append(x_pos)
            sod_values.append(group['SOD_Secondary'][i] / 1000)  # Scale to ×10³
            cla_values.append(group['CLA_Secondary'][i])
            colors.append(PC_COLORS[group['PC_type'][i]])
            system_labels.append(group['systems'][i])
            x_pos += 1
        group_boundaries.append((group_start, len(positions) - 1, group['short_name']))
        x_pos += 1  # Gap between groups
    
    positions = np.array(positions)
    
    # Create bars
    # Solid bars for Na+
    bars_sod = ax.bar(positions - bar_width/2, sod_values, bar_width,
                      color=colors, edgecolor='black', linewidth=1,
                      label='Na⁺ - Secondary Lipid (solid)')
    
    # Hatched bars for Cl- (scaled differently due to small values)
    # We'll use a secondary y-axis for Cl-
    ax2 = ax.twinx()
    light_colors = [plt.matplotlib.colors.to_rgba(c, alpha=0.6) for c in colors]
    bars_cla = ax2.bar(positions + bar_width/2, cla_values, bar_width,
                       color=light_colors, edgecolor='black', linewidth=1, hatch='///',
                       label='Cl⁻ - Secondary Lipid (hatched)')
    
    # Customize primary axis (Na+)
    ax.set_ylabel('Na⁺ Interaction Energy (×10³ kJ/mol)', fontsize=12, color='black')
    ax.tick_params(axis='y', labelcolor='black')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1)
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')
    
    # Customize secondary axis (Cl-)
    ax2.set_ylabel('Cl⁻ Interaction Energy (kJ/mol)', fontsize=12, color='gray')
    ax2.tick_params(axis='y', labelcolor='gray')
    
    # X-axis - remove system labels
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.set_xlabel('')
    
    # Title
    ax.set_title('Ion-Secondary Lipid Interaction Energies\n(Solid: Na⁺, Hatched: Cl⁻)', 
                fontweight='bold', fontsize=14)
    
    # Add group labels BELOW and separators
    for start, end, label in group_boundaries:
        mid = (positions[start] + positions[end]) / 2
        ax.text(mid, -0.15, label, ha='center', va='top', 
                fontsize=12, fontweight='bold', style='italic',
                transform=ax.get_xaxis_transform())
        if end < len(positions) - 1:
            ax.axvline(x=(positions[end] + positions[end+1]) / 2, 
                      color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='gray', edgecolor='black', label='Na⁺ (solid)'),
        mpatches.Patch(facecolor='lightgray', edgecolor='black', hatch='///', label='Cl⁻ (hatched)'),
        mpatches.Patch(facecolor=PC_COLORS['DOPC'], edgecolor='black', label='DOPC'),
        mpatches.Patch(facecolor=PC_COLORS['POPC'], edgecolor='black', label='POPC'),
        mpatches.Patch(facecolor=PC_COLORS['SOPC'], edgecolor='black', label='SOPC'),
    ]
    ax.legend(handles=legend_elements, loc='lower left', fontsize=9, framealpha=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)
    plt.savefig('Ion_Secondary_Combined.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('Ion_Secondary_Combined.pdf', facecolor='white', edgecolor='none')
    print("Saved: Ion_Secondary_Combined.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Print Summary
# =============================================================================

def print_summary():
    """
    Print summary of ion-lipid interactions.
    """
    print("\n" + "=" * 80)
    print("ION-SECONDARY LIPID INTERACTION SUMMARY (kJ/mol)")
    print("=" * 80)
    
    for group in ALL_GROUPS:
        print(f"\n{group['name']}")
        print("-" * 60)
        print(f"{'System':<15} {'Na⁺-Secondary':>15} {'Cl⁻-Secondary':>15}")
        print("-" * 60)
        
        for i, sys in enumerate(group['systems']):
            sod = group['SOD_Secondary'][i]
            cla = group['CLA_Secondary'][i]
            print(f"{sys:<15} {sod:>+15,} {cla:>+15}")
        
        print("-" * 60)
        sod_avg = np.mean(group['SOD_Secondary'])
        cla_avg = np.mean(group['CLA_Secondary'])
        print(f"{'Average':<15} {sod_avg:>+15,.0f} {cla_avg:>+15,.1f}")
    
    print("\n" + "=" * 80)
    print("KEY OBSERVATIONS:")
    print("=" * 80)
    print("""
1. SODIUM (Na⁺) INTERACTIONS:
   - Strong attractive interactions with secondary lipids (negative values)
   - Strongest in 90:10 PA systems (~-3,600 kJ/mol)
   - Moderate in 60:40 PG systems (~-3,500 kJ/mol)
   - Weakest in 80:20 PG systems (~-150 kJ/mol)
   - Na⁺ preferentially binds to anionic headgroups

2. CHLORIDE (Cl⁻) INTERACTIONS:
   - Very weak interactions with secondary lipids (~0-2 kJ/mol)
   - Slight repulsion from anionic headgroups
   - Cl⁻ is excluded from the anionic lipid surface

3. CHARGE SCREENING:
   - Na⁺ effectively screens the negative charge of PA/PG headgroups
   - This contributes to membrane stability in anionic lipid mixtures
    """)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("Ion-Secondary Lipid Interaction Energies in Mixed Bilayers")
    print("=" * 80)
    
    # Print summary
    print_summary()
    
    # Generate plots
    print("\n[1/4] Generating main ion-lipid interaction plot...")
    create_ion_lipid_plot()
    
    print("\n[2/4] Generating grouped plot by system type...")
    create_grouped_ion_plot()
    
    print("\n[3/4] Generating two-panel figure...")
    create_two_panel_ion_plot()
    
    print("\n[4/4] Generating combined plot with dual y-axis...")
    create_combined_ion_plot()
    
    print("\n" + "=" * 80)
    print("All plots generated successfully!")
    print("=" * 80)
