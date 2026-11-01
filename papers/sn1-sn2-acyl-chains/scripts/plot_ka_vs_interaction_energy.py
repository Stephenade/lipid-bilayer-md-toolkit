#!/usr/bin/env python3
"""
K_A vs Total Lipid-Lipid Interaction Energy Analysis
=====================================================
Generates publication-quality figures comparing:
1. K_A vs Total Interaction Energy for each group separately
2. K_A vs Total Interaction Energy for all 27 systems combined

Using CONSISTENT COLOR CODING from previous K_A and Order Parameter plots:
- DOPC: Blue (#3498db)
- POPC: Green (#2ecc71)
- SOPC: Purple (#9b59b6)

Author: Generated for lipid membrane MD simulation analysis
Date: February 2026
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.lines import Line2D

# =============================================================================
# PLOT SETTINGS
# =============================================================================

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 11,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 11,
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
# CONSISTENT COLOR SCHEME (from K_A and Order Parameter plots)
# =============================================================================

# PC Lipids - CONSISTENT with previous plots
PC_COLORS = {
    'DOPC': '#3498db',  # Blue
    'POPC': '#2ecc71',  # Green (was #27ae60 in some scripts, using #2ecc71 for consistency)
    'SOPC': '#9b59b6',  # Purple
}

# Anionic Lipids - CONSISTENT with previous plots
ANIONIC_COLORS = {
    'DOPA': '#e74c3c',  # Red
    'POPA': '#f39c12',  # Orange
    'SOPA': '#7f8c8d',  # Dark Gray
    'DOPG': '#e74c3c',  # Red
    'POPG': '#f39c12',  # Orange
    'SOPG': '#7f8c8d',  # Dark Gray
}

# Marker styles for groups
GROUP_MARKERS = {
    'Group 1: PC/PA (80:20)': 'o',
    'Group 2: PC/PG (65:35)': 's',
    'Group 3: PC/PG (95:5)': '^',
}

# =============================================================================
# DATA: K_A and Total Interaction Energies for All 27 Systems
# =============================================================================

# -----------------------------------------------------------------------------
# K_A, APL and error are read from the published dataset so the figures always
# match it: Adepoju & Hung, Zenodo, https://doi.org/10.5281/zenodo.23850654
# (a copy of the CSV is kept in ../data/).
# -----------------------------------------------------------------------------
import csv as _csv
from pathlib import Path as _Path

DATASET_CSV = _Path(__file__).resolve().parent.parent / "data" / "KA_PC_PA_PG_Mixtures.csv"

def _load_dataset(path=DATASET_CSV):
    table = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in _csv.DictReader(fh):
            key = (row["PC_lipid"], row["Secondary_lipid"], int(row["PC_mol_percent"]))
            table[key] = {"KA": float(row["KA_mN_per_m"]),
                          "APL": float(row["APL_A2"]),
                          "error": float(row["KA_error_mN_per_m"])}
    return table

_DATASET = _load_dataset()

# Group 1: PC/PA (80/20)
group1_data = {
    'name': 'Group 1: PC/PA (80:20)',
    'short_name': '80:20 PC/PA',
    'systems': [
        'DOPC_DOPA', 'DOPC_POPA', 'DOPC_SOPA',
        'POPC_DOPA', 'POPC_POPA', 'POPC_SOPA',
        'SOPC_DOPA', 'SOPC_POPA', 'SOPC_SOPA'
    ],
    'KA': None,  # filled from the dataset below (80:20)
    'total_energy': [-47342, -46527, -46976, -42529, -42433, -42658, -44616, -44183, -44081],  # kJ/mol
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
}

# Group 2: PC/PG (65/35)
group2_data = {
    'name': 'Group 2: PC/PG (65:35)',
    'short_name': '65:35 PC/PG',
    'systems': [
        'DOPC_DOPG', 'DOPC_POPG', 'DOPC_SOPG',
        'POPC_DOPG', 'POPC_POPG', 'POPC_SOPG',
        'SOPC_DOPG', 'SOPC_POPG', 'SOPC_SOPG'
    ],
    'KA': None,  # filled from the dataset below (65:35)
    'total_energy': [-119294, -117735, -117933, -115997, -114315, -114424, -116937, -115235, -115400],  # kJ/mol
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
}

# Group 3: PC/PG (95/5)
group3_data = {
    'name': 'Group 3: PC/PG (95:5)',
    'short_name': '95:5 PC/PG',
    'systems': [
        'DOPC_DOPG', 'DOPC_POPG', 'DOPC_SOPG',
        'POPC_DOPG', 'POPC_POPG', 'POPC_SOPG',
        'SOPC_DOPG', 'SOPC_POPG', 'SOPC_SOPG'
    ],
    'KA': None,  # filled from the dataset below (95:5)
    'total_energy': [-11226, -10938, -11082, -6451, -6241, -6222, -7746, -7740, -7746],  # kJ/mol
    'PC_type': ['DOPC', 'DOPC', 'DOPC', 'POPC', 'POPC', 'POPC', 'SOPC', 'SOPC', 'SOPC'],
}


# Fill K_A for every system from the published dataset
for _g, _pct in ((group1_data, 80), (group2_data, 65), (group3_data, 95)):
    _g['KA'] = [_DATASET[(*_sys.split('_'), _pct)]['KA'] for _sys in _g['systems']]


# =============================================================================
# FUNCTION: Plot Single Group
# =============================================================================

def plot_single_group(group_data, ax, show_labels=True):
    """
    Create a scatter plot of K_A vs Total Interaction Energy for a single group.
    Points are colored by PC type using consistent color scheme.
    X-axis scaled by 10^4 for cleaner display.
    """
    # Convert to numpy arrays
    ka = np.array(group_data['KA'])
    energy = np.array(group_data['total_energy'])
    energy_scaled = energy / 10000  # Scale by 10^4
    pc_types = group_data['PC_type']
    systems = group_data['systems']
    
    # Plot each point colored by PC type
    for pc in ['DOPC', 'POPC', 'SOPC']:
        mask = [p == pc for p in pc_types]
        ax.scatter(np.array(energy_scaled)[mask], np.array(ka)[mask], 
                   c=PC_COLORS[pc], s=120, edgecolors='black', 
                   linewidth=1, zorder=5, label=pc)
    
    # Add system labels if requested
    if show_labels:
        for k, e, sys in zip(ka, energy_scaled, systems):
            short_name = sys.split('_')[1]  # Get anionic lipid name
            ax.annotate(short_name, (e, k), xytext=(5, 5), textcoords='offset points', 
                        fontsize=8, alpha=0.7)
    
    # Linear regression (on scaled data for plotting)
    slope, intercept, r_value, p_value, std_err = stats.linregress(energy_scaled, ka)
    x_line = np.linspace(min(energy_scaled) - abs(min(energy_scaled))*0.02, 
                         max(energy_scaled) + abs(max(energy_scaled))*0.02, 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'k--', linewidth=1.5, alpha=0.7)
    
    # Add regression info
    eq_text = f'R² = {r_value**2:.3f}\np = {p_value:.2e}'
    ax.text(0.05, 0.95, eq_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Labels and title
    ax.set_xlabel('Total Lipid-Lipid Interaction Energy (×10$^4$ kJ/mol)')
    ax.set_ylabel('K$_A$ (mN/m)')
    ax.set_title(group_data['name'], fontweight='bold')
    ax.legend(loc='lower right', title='PC Lipid')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Return original slope (per kJ/mol, not per 10^4 kJ/mol)
    return slope/10000, intercept, r_value, p_value


# =============================================================================
# FUNCTION: Plot All Three Groups Separately (Figure 1)
# =============================================================================

def plot_all_groups_separate():
    """
    Create a figure with 3 subplots, one for each group.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    groups = [group1_data, group2_data, group3_data]
    
    print("\n" + "=" * 70)
    print("CORRELATION ANALYSIS: K_A vs Total Lipid-Lipid Interaction Energy")
    print("=" * 70)
    
    for ax, group in zip(axes, groups):
        slope, intercept, r_value, p_value = plot_single_group(group, ax)
        
        print(f"\n{group['name']}")
        print("-" * 50)
        print(f"  Slope: {slope:.6f} (mN/m) / (kJ/mol)")
        print(f"  Intercept: {intercept:.2f} mN/m")
        print(f"  R²: {r_value**2:.4f}")
        print(f"  p-value: {p_value:.2e}")
        if p_value < 0.05:
            print(f"  *** Statistically significant (p < 0.05) ***")
    
    plt.tight_layout()
    plt.savefig('KA_vs_Energy_Three_Groups.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('KA_vs_Energy_Three_Groups.pdf', facecolor='white', edgecolor='none')
    print("\nSaved: KA_vs_Energy_Three_Groups.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Plot All 27 Systems Combined (Figure 2)
# =============================================================================

def plot_all_systems_combined():
    """
    Create a single plot with all 27 systems, colored by PC type and shaped by group.
    X-axis scaled by 10^4 for cleaner display.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    groups = [group1_data, group2_data, group3_data]
    
    # Collect all data
    all_ka = []
    all_energy = []
    
    # Plot each group with different markers
    for group in groups:
        ka = np.array(group['KA'])
        energy = np.array(group['total_energy'])
        energy_scaled = energy / 10000  # Scale by 10^4
        pc_types = group['PC_type']
        group_name = group['name']
        marker = GROUP_MARKERS[group_name]
        
        all_ka.extend(ka)
        all_energy.extend(energy)
        
        # Plot each point
        for k, e, pc in zip(ka, energy_scaled, pc_types):
            ax.scatter(e, k, c=PC_COLORS[pc], marker=marker, s=150, 
                       edgecolors='black', linewidth=1, zorder=5)
    
    # Overall linear regression (on scaled data)
    all_ka = np.array(all_ka)
    all_energy = np.array(all_energy)
    all_energy_scaled = all_energy / 10000
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_energy_scaled, all_ka)
    
    x_line = np.linspace(min(all_energy_scaled) - 0.5, max(all_energy_scaled) + 0.3, 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'k--', linewidth=2, alpha=0.7, label=f'Linear fit (R² = {r_value**2:.3f})')
    
    # Add regression equation (convert slope back to per kJ/mol for display)
    slope_original = slope / 10000
    eq_text = f'K$_A$ = {slope_original:.5f} × E + {intercept:.1f}\nR² = {r_value**2:.3f}, p = {p_value:.2e}'
    ax.text(0.02, 0.98, eq_text, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
    
    # Create custom legend
    # PC type legend (colors) - using consistent colors
    pc_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=PC_COLORS[pc], 
                         markersize=12, markeredgecolor='black', label=pc) 
                  for pc in ['DOPC', 'POPC', 'SOPC']]
    
    # Group legend (markers)
    group_handles = [Line2D([0], [0], marker=GROUP_MARKERS[g['name']], color='w', 
                            markerfacecolor='gray', markersize=12, markeredgecolor='black', 
                            label=g['short_name']) for g in groups]
    
    # Create two legends
    legend1 = ax.legend(handles=pc_handles, loc='lower left', title='PC Lipid', framealpha=0.9)
    ax.add_artist(legend1)
    ax.legend(handles=group_handles, loc='lower right', title='Group', framealpha=0.9)
    
    # Labels and title
    ax.set_xlabel('Total Lipid-Lipid Interaction Energy (×10$^4$ kJ/mol)', fontsize=14)
    ax.set_ylabel('Area Compressibility Modulus, K$_A$ (mN/m)', fontsize=14)
    ax.set_title('K$_A$ vs Total Interaction Energy: All 27 Systems', fontsize=16, fontweight='bold')
    
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Set axis limits with padding
    ax.set_xlim(min(all_energy_scaled) - 0.8, max(all_energy_scaled) + 0.3)
    ax.set_ylim(min(all_ka) - 10, max(all_ka) + 10)
    
    plt.tight_layout()
    plt.savefig('KA_vs_Energy_All_27_Systems.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('KA_vs_Energy_All_27_Systems.pdf', facecolor='white', edgecolor='none')
    
    print(f"\n{'=' * 50}")
    print("ALL 27 SYSTEMS COMBINED")
    print("=" * 50)
    print(f"  Slope: {slope_original:.6f} (mN/m) / (kJ/mol)")
    print(f"  Intercept: {intercept:.2f} mN/m")
    print(f"  R²: {r_value**2:.4f}")
    print(f"  p-value: {p_value:.2e}")
    
    print("\nSaved: KA_vs_Energy_All_27_Systems.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Plot 4-Panel Detailed Figure (Figure 3)
# =============================================================================

def plot_four_panel():
    """
    Create a 2x2 figure with individual groups and combined plot.
    X-axis scaled by 10^4 for cleaner display.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    groups = [group1_data, group2_data, group3_data]
    
    # Top row and bottom left: Individual groups
    plot_single_group(group1_data, axes[0, 0], show_labels=True)
    plot_single_group(group2_data, axes[0, 1], show_labels=True)
    plot_single_group(group3_data, axes[1, 0], show_labels=True)
    
    # Bottom right: Combined plot
    ax = axes[1, 1]
    
    all_ka = []
    all_energy = []
    
    for group in groups:
        ka = np.array(group['KA'])
        energy = np.array(group['total_energy'])
        energy_scaled = energy / 10000  # Scale by 10^4
        pc_types = group['PC_type']
        group_name = group['name']
        marker = GROUP_MARKERS[group_name]
        
        all_ka.extend(ka)
        all_energy.extend(energy)
        
        for k, e, pc in zip(ka, energy_scaled, pc_types):
            ax.scatter(e, k, c=PC_COLORS[pc], marker=marker, s=100, 
                       edgecolors='black', linewidth=1, zorder=5)
    
    # Regression (on scaled data)
    all_ka = np.array(all_ka)
    all_energy = np.array(all_energy)
    all_energy_scaled = all_energy / 10000
    slope, intercept, r_value, p_value, std_err = stats.linregress(all_energy_scaled, all_ka)
    x_line = np.linspace(min(all_energy_scaled) - 0.5, max(all_energy_scaled) + 0.3, 100)
    y_line = slope * x_line + intercept
    ax.plot(x_line, y_line, 'k--', linewidth=2, alpha=0.7)
    
    eq_text = f'R² = {r_value**2:.3f}\np = {p_value:.2e}'
    ax.text(0.05, 0.95, eq_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9))
    
    ax.set_xlabel('Total Interaction Energy (×10$^4$ kJ/mol)')
    ax.set_ylabel('K$_A$ (mN/m)')
    ax.set_title('All 27 Systems Combined', fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Legend
    pc_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=PC_COLORS[pc], 
                         markersize=8, markeredgecolor='black', label=pc) 
                  for pc in ['DOPC', 'POPC', 'SOPC']]
    group_handles = [Line2D([0], [0], marker=GROUP_MARKERS[g['name']], color='w', 
                            markerfacecolor='gray', markersize=8, markeredgecolor='black', 
                            label=g['short_name']) for g in groups]
    
    legend1 = ax.legend(handles=pc_handles, loc='lower left', title='PC', fontsize=8)
    ax.add_artist(legend1)
    ax.legend(handles=group_handles, loc='lower right', title='Group', fontsize=8)
    
    # Add panel labels (A, B, C, D)
    for i, (ax_i, label) in enumerate(zip(axes.flat, ['A', 'B', 'C', 'D'])):
        ax_i.text(-0.12, 1.05, label, transform=ax_i.transAxes, fontsize=16, 
                fontweight='bold', va='top')
    
    plt.tight_layout()
    plt.savefig('KA_vs_Energy_4Panel.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('KA_vs_Energy_4Panel.pdf', facecolor='white', edgecolor='none')
    print("\nSaved: KA_vs_Energy_4Panel.png/pdf")
    plt.show()


# =============================================================================
# FUNCTION: Print Interpretation
# =============================================================================

def print_interpretation():
    """
    Print the regression statistics computed from the data above.
    (Values are calculated at run time, so they always reflect the K_A in the
    dataset CSV and the energies in this script.)
    """
    print("\n" + "=" * 70)
    print("REGRESSION SUMMARY: K_A vs total lipid-lipid interaction energy")
    print("=" * 70)
    all_e, all_ka = [], []
    for g in (group1_data, group2_data, group3_data):
        e = np.array(g['total_energy'], dtype=float)
        ka = np.array(g['KA'], dtype=float)
        all_e.extend(e); all_ka.extend(ka)
        r = stats.linregress(e, ka)
        print(f"  {g['name']:<26} slope = {r.slope:+.6f} (mN/m)/(kJ/mol)   "
              f"R^2 = {r.rvalue**2:.3f}   p = {r.pvalue:.2e}")
    r = stats.linregress(np.array(all_e), np.array(all_ka))
    print(f"  {'All 27 systems':<26} slope = {r.slope:+.6f} (mN/m)/(kJ/mol)   "
          f"R^2 = {r.rvalue**2:.3f}   p = {r.pvalue:.2e}")
    print("\nColour coding: DOPC blue (#3498db), POPC green (#2ecc71), SOPC purple (#9b59b6)")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("K_A vs Lipid-Lipid Interaction Energy Analysis")
    print("(Using consistent color coding from K_A and Order Parameter plots)")
    print("=" * 70)
    
    # Figure 1: Three groups separately
    print("\n[1/3] Generating three-group comparison plot...")
    plot_all_groups_separate()
    
    # Figure 2: All 27 systems combined
    print("\n[2/3] Generating combined 27-system plot...")
    plot_all_systems_combined()
    
    # Figure 3: Four-panel detailed figure
    print("\n[3/3] Generating four-panel detailed figure...")
    plot_four_panel()
    
    # Print interpretation
    print_interpretation()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
