
#!/usr/bin/env python3
"""
Area Compressibility Modulus (K_A) Comparison Plots
====================================================
This script generates publication-quality figures comparing K_A across:
1. Group 1: PC/PA systems (80:20 ratio)
2. Group 2: PC/PG systems (65:35 ratio)
3. Group 3: PC/PG systems (95:5 ratio)
4. All 27 systems combined

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
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# =============================================================================
# DATA: Area Compressibility Modulus (K_A) and Area Per Lipid (APL)
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

def _pc_block(pc, anionic, pc_percent):
    rows = [_DATASET[(pc, a, pc_percent)] for a in anionic]
    return {k: [r[k] for r in rows] for k in ("KA", "APL", "error")}

# Group 1: PC/PA Systems (80:20 ratio)
group1_data = {
    'name': 'Group 1: PC/PA (80:20)',
    'anionic_lipids': ['DOPA', 'POPA', 'SOPA'],
    'DOPC': _pc_block('DOPC', ['DOPA', 'POPA', 'SOPA'], 80),
    'POPC': _pc_block('POPC', ['DOPA', 'POPA', 'SOPA'], 80),
    'SOPC': _pc_block('SOPC', ['DOPA', 'POPA', 'SOPA'], 80),
}

# Group 2: PC/PG Systems (65:35 ratio)
group2_data = {
    'name': 'Group 2: PC/PG (65:35)',
    'anionic_lipids': ['DOPG', 'POPG', 'SOPG'],
    'DOPC': _pc_block('DOPC', ['DOPG', 'POPG', 'SOPG'], 65),
    'POPC': _pc_block('POPC', ['DOPG', 'POPG', 'SOPG'], 65),
    'SOPC': _pc_block('SOPC', ['DOPG', 'POPG', 'SOPG'], 65),
}

# Group 3: PC/PG Systems (95:5 ratio)
group3_data = {
    'name': 'Group 3: PC/PG (95:5)',
    'anionic_lipids': ['DOPG', 'POPG', 'SOPG'],
    'DOPC': _pc_block('DOPC', ['DOPG', 'POPG', 'SOPG'], 95),
    'POPC': _pc_block('POPC', ['DOPG', 'POPG', 'SOPG'], 95),
    'SOPC': _pc_block('SOPC', ['DOPG', 'POPG', 'SOPG'], 95),
}

# Reference value
PURE_DOPC_KA = 245  # mN/m (literature value)

# Colors
PC_COLORS = {
    'DOPC': '#3498db',  # Blue
    'POPC': '#2ecc71',  # Green
    'SOPC': '#9b59b6',  # Purple
}


# =============================================================================
# PLOT 1: Group 1 - PC/PA Systems (80:20)
# =============================================================================

def plot_group1():
    """Plot K_A comparison for Group 1: PC/PA (80:20) systems."""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    data = group1_data
    anionic_lipids = data['anionic_lipids']
    x = np.arange(len(anionic_lipids))
    width = 0.25
    
    # Create bars
    bars1 = ax.bar(x - width, data['DOPC']['KA'], width, 
                   yerr=data['DOPC']['error'],
                   label='DOPC mixtures', color=PC_COLORS['DOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars2 = ax.bar(x, data['POPC']['KA'], width,
                   yerr=data['POPC']['error'],
                   label='POPC mixtures', color=PC_COLORS['POPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars3 = ax.bar(x + width, data['SOPC']['KA'], width,
                   yerr=data['SOPC']['error'],
                   label='SOPC mixtures', color=PC_COLORS['SOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    # Reference line for pure DOPC
    ax.axhline(y=PURE_DOPC_KA, color='#e74c3c', linestyle='--', linewidth=2.5,
               label=f'Pure DOPC ({PURE_DOPC_KA} mN/m)')
    
    # Add value labels on bars
    def add_bar_labels(bars, errors):
        for bar, err in zip(bars, errors):
            height = bar.get_height()
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + err + 3),
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    add_bar_labels(bars1, data['DOPC']['error'])
    add_bar_labels(bars2, data['POPC']['error'])
    add_bar_labels(bars3, data['SOPC']['error'])
    
    # Formatting
    ax.set_xlabel('Anionic Lipid (PA)', fontsize=13, fontweight='bold')
    ax.set_ylabel('$K_A$ (mN/m)', fontsize=13, fontweight='bold')
    ax.set_title('Group 1: Area Compressibility Modulus ($K_A$)\nPC/PA Systems (80:20 Ratio)', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(anionic_lipids, fontsize=12, fontweight='bold')
    ax.set_ylim(0, 280)
    ax.set_xlim(-0.5, 2.5)
    
    # Grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    # Legend
    ax.legend(loc='upper right', fontsize=11, framealpha=1, edgecolor='black')
    
    # Add summary statistics box
    dopc_avg = np.mean(data['DOPC']['KA'])
    popc_avg = np.mean(data['POPC']['KA'])
    sopc_avg = np.mean(data['SOPC']['KA'])
    
    textstr = f'Average $K_A$ (mN/m):\nDOPC: {dopc_avg:.0f}\nPOPC: {popc_avg:.0f}\nSOPC: {sopc_avg:.0f}'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('KA_Group1_PC_PA_80_20.png', dpi=300, facecolor='white')
    plt.savefig('KA_Group1_PC_PA_80_20.pdf', facecolor='white')
    print("Saved: KA_Group1_PC_PA_80_20.png/pdf")
    plt.show()


# =============================================================================
# PLOT 2: Group 2 - PC/PG Systems (65:35)
# =============================================================================

def plot_group2():
    """Plot K_A comparison for Group 2: PC/PG (65:35) systems."""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    data = group2_data
    anionic_lipids = data['anionic_lipids']
    x = np.arange(len(anionic_lipids))
    width = 0.25
    
    # Create bars
    bars1 = ax.bar(x - width, data['DOPC']['KA'], width, 
                   yerr=data['DOPC']['error'],
                   label='DOPC mixtures', color=PC_COLORS['DOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars2 = ax.bar(x, data['POPC']['KA'], width,
                   yerr=data['POPC']['error'],
                   label='POPC mixtures', color=PC_COLORS['POPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars3 = ax.bar(x + width, data['SOPC']['KA'], width,
                   yerr=data['SOPC']['error'],
                   label='SOPC mixtures', color=PC_COLORS['SOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    # Reference line for pure DOPC
    ax.axhline(y=PURE_DOPC_KA, color='#e74c3c', linestyle='--', linewidth=2.5,
               label=f'Pure DOPC ({PURE_DOPC_KA} mN/m)')
    
    # Add value labels on bars
    def add_bar_labels(bars, errors):
        for bar, err in zip(bars, errors):
            height = bar.get_height()
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + err + 3),
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    add_bar_labels(bars1, data['DOPC']['error'])
    add_bar_labels(bars2, data['POPC']['error'])
    add_bar_labels(bars3, data['SOPC']['error'])
    
    # Formatting
    ax.set_xlabel('Anionic Lipid (PG)', fontsize=13, fontweight='bold')
    ax.set_ylabel('$K_A$ (mN/m)', fontsize=13, fontweight='bold')
    ax.set_title('Group 2: Area Compressibility Modulus ($K_A$)\nPC/PG Systems (65:35 Ratio)', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(anionic_lipids, fontsize=12, fontweight='bold')
    ax.set_ylim(0, 280)
    ax.set_xlim(-0.5, 2.5)
    
    # Grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    # Legend
    ax.legend(loc='upper right', fontsize=11, framealpha=1, edgecolor='black')
    
    # Add summary statistics box
    dopc_avg = np.mean(data['DOPC']['KA'])
    popc_avg = np.mean(data['POPC']['KA'])
    sopc_avg = np.mean(data['SOPC']['KA'])
    
    textstr = f'Average $K_A$ (mN/m):\nDOPC: {dopc_avg:.0f}\nPOPC: {popc_avg:.0f}\nSOPC: {sopc_avg:.0f}'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('KA_Group2_PC_PG_65_35.png', dpi=300, facecolor='white')
    plt.savefig('KA_Group2_PC_PG_65_35.pdf', facecolor='white')
    print("Saved: KA_Group2_PC_PG_65_35.png/pdf")
    plt.show()


# =============================================================================
# PLOT 3: Group 3 - PC/PG Systems (95:5)
# =============================================================================

def plot_group3():
    """Plot K_A comparison for Group 3: PC/PG (95:5) systems."""
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    data = group3_data
    anionic_lipids = data['anionic_lipids']
    x = np.arange(len(anionic_lipids))
    width = 0.25
    
    # Create bars
    bars1 = ax.bar(x - width, data['DOPC']['KA'], width, 
                   yerr=data['DOPC']['error'],
                   label='DOPC mixtures', color=PC_COLORS['DOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars2 = ax.bar(x, data['POPC']['KA'], width,
                   yerr=data['POPC']['error'],
                   label='POPC mixtures', color=PC_COLORS['POPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    bars3 = ax.bar(x + width, data['SOPC']['KA'], width,
                   yerr=data['SOPC']['error'],
                   label='SOPC mixtures', color=PC_COLORS['SOPC'],
                   edgecolor='black', linewidth=1.2,
                   capsize=4, error_kw={'elinewidth': 1.5, 'capthick': 1.5})
    
    # Reference line for pure DOPC
    ax.axhline(y=PURE_DOPC_KA, color='#e74c3c', linestyle='--', linewidth=2.5,
               label=f'Pure DOPC ({PURE_DOPC_KA} mN/m)')
    
    # Add value labels on bars
    def add_bar_labels(bars, errors):
        for bar, err in zip(bars, errors):
            height = bar.get_height()
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + err + 3),
                       ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    add_bar_labels(bars1, data['DOPC']['error'])
    add_bar_labels(bars2, data['POPC']['error'])
    add_bar_labels(bars3, data['SOPC']['error'])
    
    # Formatting
    ax.set_xlabel('Anionic Lipid (PG)', fontsize=13, fontweight='bold')
    ax.set_ylabel('$K_A$ (mN/m)', fontsize=13, fontweight='bold')
    ax.set_title('Group 3: Area Compressibility Modulus ($K_A$)\nPC/PG Systems (95:5 Ratio)', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(anionic_lipids, fontsize=12, fontweight='bold')
    ax.set_ylim(0, 280)
    ax.set_xlim(-0.5, 2.5)
    
    # Grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    # Legend
    ax.legend(loc='upper right', fontsize=11, framealpha=1, edgecolor='black')
    
    # Add summary statistics box
    dopc_avg = np.mean(data['DOPC']['KA'])
    popc_avg = np.mean(data['POPC']['KA'])
    sopc_avg = np.mean(data['SOPC']['KA'])
    
    textstr = f'Average $K_A$ (mN/m):\nDOPC: {dopc_avg:.0f}\nPOPC: {popc_avg:.0f}\nSOPC: {sopc_avg:.0f}'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray')
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig('KA_Group3_PC_PG_95_5.png', dpi=300, facecolor='white')
    plt.savefig('KA_Group3_PC_PG_95_5.pdf', facecolor='white')
    print("Saved: KA_Group3_PC_PG_95_5.png/pdf")
    plt.show()


# =============================================================================
# PLOT 4: All 27 Systems Combined
# =============================================================================

def plot_all_27_systems():
    """Plot K_A comparison for all 27 systems in one figure."""
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # X-axis labels
    x_labels = ['DOPA\n(80:20)', 'POPA\n(80:20)', 'SOPA\n(80:20)', 
                'DOPG\n(65:35)', 'POPG\n(65:35)', 'SOPG\n(65:35)',
                'DOPG\n(95:5)', 'POPG\n(95:5)', 'SOPG\n(95:5)']
    
    # Combine K_A data from all groups
    dopc_ka = group1_data['DOPC']['KA'] + group2_data['DOPC']['KA'] + group3_data['DOPC']['KA']
    dopc_err = group1_data['DOPC']['error'] + group2_data['DOPC']['error'] + group3_data['DOPC']['error']
    
    popc_ka = group1_data['POPC']['KA'] + group2_data['POPC']['KA'] + group3_data['POPC']['KA']
    popc_err = group1_data['POPC']['error'] + group2_data['POPC']['error'] + group3_data['POPC']['error']
    
    sopc_ka = group1_data['SOPC']['KA'] + group2_data['SOPC']['KA'] + group3_data['SOPC']['KA']
    sopc_err = group1_data['SOPC']['error'] + group2_data['SOPC']['error'] + group3_data['SOPC']['error']
    
    x = np.arange(len(x_labels))
    width = 0.25
    
    # Create bars
    bars1 = ax.bar(x - width, dopc_ka, width, 
                   yerr=dopc_err,
                   label='DOPC mixtures', color=PC_COLORS['DOPC'],
                   edgecolor='black', linewidth=1,
                   capsize=3, error_kw={'elinewidth': 1.2, 'capthick': 1.2})
    
    bars2 = ax.bar(x, popc_ka, width,
                   yerr=popc_err,
                   label='POPC mixtures', color=PC_COLORS['POPC'],
                   edgecolor='black', linewidth=1,
                   capsize=3, error_kw={'elinewidth': 1.2, 'capthick': 1.2})
    
    bars3 = ax.bar(x + width, sopc_ka, width,
                   yerr=sopc_err,
                   label='SOPC mixtures', color=PC_COLORS['SOPC'],
                   edgecolor='black', linewidth=1,
                   capsize=3, error_kw={'elinewidth': 1.2, 'capthick': 1.2})
    
    # Reference line for pure DOPC
    ax.axhline(y=PURE_DOPC_KA, color='#e74c3c', linestyle='--', linewidth=2.5,
               label=f'Pure DOPC ({PURE_DOPC_KA} mN/m)')
    
    # Add value labels on bars (rotated for space)
    def add_bar_labels(bars, errors):
        for bar, err in zip(bars, errors):
            height = bar.get_height()
            ax.annotate(f'{int(height)}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + err + 2),
                       ha='center', va='bottom', fontsize=8, fontweight='bold',
                       rotation=90)
    
    add_bar_labels(bars1, dopc_err)
    add_bar_labels(bars2, popc_err)
    add_bar_labels(bars3, sopc_err)
    
    # Add vertical separators between groups
    ax.axvline(x=2.5, color='gray', linestyle='-', linewidth=2, alpha=0.5)
    ax.axvline(x=5.5, color='gray', linestyle='-', linewidth=2, alpha=0.5)
    
    # Add group labels at top
    ax.text(1, 270, 'Group 1\nPC/PA (80:20)', ha='center', fontsize=11, 
            fontweight='bold', color='#2c3e50')
    ax.text(4, 270, 'Group 2\nPC/PG (65:35)', ha='center', fontsize=11, 
            fontweight='bold', color='#2c3e50')
    ax.text(7, 270, 'Group 3\nPC/PG (95:5)', ha='center', fontsize=11, 
            fontweight='bold', color='#2c3e50')
    
    # Add shaded regions for groups
    ax.axvspan(-0.5, 2.5, alpha=0.08, color='blue')
    ax.axvspan(2.5, 5.5, alpha=0.08, color='green')
    ax.axvspan(5.5, 8.5, alpha=0.08, color='purple')
    
    # Formatting
    ax.set_xlabel('Anionic Lipid (Concentration)', fontsize=13, fontweight='bold')
    ax.set_ylabel('$K_A$ (mN/m)', fontsize=13, fontweight='bold')
    ax.set_title('Area Compressibility Modulus ($K_A$) Comparison\nAll 27 Lipid Membrane Systems', 
                 fontsize=15, fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10, fontweight='bold')
    ax.set_ylim(0, 290)
    ax.set_xlim(-0.5, 8.5)
    
    # Grid
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)
    
    # Spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    # Legend
    ax.legend(loc='upper right', fontsize=11, framealpha=1, edgecolor='black')
    
    plt.tight_layout()
    plt.savefig('KA_All_27_Systems.png', dpi=300, facecolor='white')
    plt.savefig('KA_All_27_Systems.pdf', facecolor='white')
    print("Saved: KA_All_27_Systems.png/pdf")
    plt.show()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Area Compressibility Modulus (K_A) Comparison Plots")
    print("=" * 70)
    print("\nGenerating plots...\n")
    
    # Plot 1: Group 1
    print("1. Plotting Group 1: PC/PA Systems (80:20)")
    plot_group1()
    
    # Plot 2: Group 2
    print("\n2. Plotting Group 2: PC/PG Systems (65:35)")
    plot_group2()
    
    # Plot 3: Group 3
    print("\n3. Plotting Group 3: PC/PG Systems (95:5)")
    plot_group3()
    
    # Plot 4: All 27 Systems
    print("\n4. Plotting All 27 Systems Combined")
    plot_all_27_systems()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
    print("\nOutput files:")
    print("  1. KA_Group1_PC_PA_80_20.png/pdf")
    print("  2. KA_Group2_PC_PG_65_35.png/pdf")
    print("  3. KA_Group3_PC_PG_95_5.png/pdf")
    print("  4. KA_All_27_Systems.png/pdf")
    print("=" * 70)
