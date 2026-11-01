#!/usr/bin/env python3
"""
K_A Comparison Plot for 17 CG Systems - Publication Quality
============================================================
Compares the area compressibility modulus (K_A) of 17 coarse-grained
binary lipid mixtures with pure DOPC as reference.

Selected systems for AA simulations:
- DOPC-DOPA (80:20) → Group 1
- DOPC-DOPG (65:35) → Group 2
- DOPC-DOPG (95:5)  → Group 3

Color scheme matches AA analysis:
- DOPA systems: Red tones
- DOPG systems: Blue tones

Author: Generated for lipid membrane MD simulation analysis
Date: February 2026
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# =============================================================================
# PLOT SETTINGS - Publication Quality
# =============================================================================

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# =============================================================================
# COLOR SCHEME - Matching AA Analysis
# =============================================================================

# DOPA systems - Red/Orange tones
DOPA_COLOR = '#e74c3c'           # Red
DOPA_SELECTED_COLOR = '#c0392b'  # Darker red for selected

# DOPG systems - Blue tones  
DOPG_COLOR = '#3498db'           # Blue
DOPG_SELECTED_COLOR = '#2980b9'  # Darker blue for selected

# Reference
PURE_DOPC_COLOR = '#2c3e50'      # Dark gray

# =============================================================================
# DATA: K_A VALUES FROM CG SIMULATIONS
# =============================================================================

pure_dopc_ka = 183.8631644

# DOPC-DOPA systems
dopa_systems = {
    '60:40': 192.6517788,
    '65:35': 189.8053147,
    '70:30': 192.6517788,
    '75:25': 192.1922541,
    '80:20': 178.0213748,   # SELECTED FOR AA
    '85:15': 190.433088,
    '90:10': 183.9681571,
    '95:5': 199.6700376,
}

# DOPC-DOPG systems
dopg_systems = {
    '60:40': 193.2027203,
    '65:35': 179.8043954,   # SELECTED FOR AA
    '70:30': 196.4631106,
    '75:25': 204.841617,
    '80:20': 190.888853,
    '85:15': 200.9016407,
    '90:10': 192.5373272,
    '95:5': 178.5696663,    # SELECTED FOR AA
}

# Selected systems for AA study
SELECTED_SYSTEMS = {
    'DOPA': ['80:20'],
    'DOPG': ['65:35', '95:5'],
}


# =============================================================================
# PLOT: Publication Quality - All Systems with Selected Highlighted
# =============================================================================

def create_publication_plot():
    """
    Create publication-quality bar plot with selected systems highlighted.
    """
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Prepare data
    systems = []
    values = []
    colors = []
    edge_colors = []
    edge_widths = []
    hatches = []
    
    # Pure DOPC
    systems.append('Pure\nDOPC')
    values.append(pure_dopc_ka)
    colors.append(PURE_DOPC_COLOR)
    edge_colors.append('black')
    edge_widths.append(1.5)
    hatches.append('')
    
    # Add DOPA systems
    for ratio, ka in dopa_systems.items():
        systems.append(f'{ratio}')
        values.append(ka)
        if ratio in SELECTED_SYSTEMS['DOPA']:
            colors.append(DOPA_SELECTED_COLOR)
            edge_colors.append('#f1c40f')  # Gold border
            edge_widths.append(4)
            hatches.append('')
        else:
            colors.append(DOPA_COLOR)
            edge_colors.append('black')
            edge_widths.append(1.5)
            hatches.append('')
    
    # Add DOPG systems
    for ratio, ka in dopg_systems.items():
        systems.append(f'{ratio}')
        values.append(ka)
        if ratio in SELECTED_SYSTEMS['DOPG']:
            colors.append(DOPG_SELECTED_COLOR)
            edge_colors.append('#f1c40f')  # Gold border
            edge_widths.append(4)
            hatches.append('')
        else:
            colors.append(DOPG_COLOR)
            edge_colors.append('black')
            edge_widths.append(1.5)
            hatches.append('')
    
    x = np.arange(len(systems))
    
    # Create bars
    bars = ax.bar(x, values, width=0.75, color=colors, edgecolor=edge_colors, 
                  linewidth=edge_widths)
    
    # Apply edge widths individually
    for bar, ec, ew in zip(bars, edge_colors, edge_widths):
        bar.set_edgecolor(ec)
        bar.set_linewidth(ew)
    
    # Add pure DOPC reference line
    ax.axhline(y=pure_dopc_ka, color='black', linestyle='--', linewidth=2, 
               alpha=0.7, zorder=1)
    ax.text(len(systems) - 0.5, pure_dopc_ka + 1, f'Pure DOPC = {pure_dopc_ka:.1f} mN/m', 
            ha='right', va='bottom', fontsize=11, style='italic')
    
    # Add vertical separators between groups
    ax.axvline(x=0.5, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.axvline(x=8.5, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    
    # Add group labels at top
    ax.text(0, 222, 'Ref', ha='center', fontsize=12, fontweight='bold', color=PURE_DOPC_COLOR)
    ax.text(4.5, 222, 'DOPC/DOPA', ha='center', fontsize=14, fontweight='bold', color=DOPA_COLOR)
    ax.text(13, 222, 'DOPC/DOPG', ha='center', fontsize=14, fontweight='bold', color=DOPG_COLOR)
    
    # Mark selected systems with stars
    selected_indices = [5, 10, 16]  # 80:20 DOPA, 65:35 DOPG, 95:5 DOPG
    for idx in selected_indices:
        ax.annotate('★', xy=(idx, values[idx] + 2), ha='center', fontsize=16, 
                    color='#f1c40f', fontweight='bold')
    
    # Customize axes
    ax.set_ylabel('Area Compressibility Modulus, $K_A$ (mN/m)', fontsize=14, fontweight='bold')
    ax.set_xlabel('PC : Secondary Lipid Ratio', fontsize=14, fontweight='bold')
    ax.set_title('Coarse-Grained MD Screening: Selection of Systems for All-Atom Simulations',
                 fontweight='bold', fontsize=16, pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(systems, fontsize=10)
    ax.set_ylim(165, 225)
    ax.set_xlim(-0.6, len(systems) - 0.4)
    
    # Grid
    ax.yaxis.grid(True, alpha=0.3, linestyle='--', zorder=0)
    ax.set_axisbelow(True)
    
    # Legend
    legend_elements = [
        Patch(facecolor=PURE_DOPC_COLOR, edgecolor='black', linewidth=1.5, label='Pure DOPC'),
        Patch(facecolor=DOPA_COLOR, edgecolor='black', linewidth=1.5, label='DOPC/DOPA'),
        Patch(facecolor=DOPG_COLOR, edgecolor='black', linewidth=1.5, label='DOPC/DOPG'),
        Patch(facecolor='white', edgecolor='#f1c40f', linewidth=4, label='Selected for AA'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=11, 
              framealpha=0.95, edgecolor='black')
    
    # Add annotation box for selected systems
    textstr = 'Selected for AA Simulations:\n• DOPC/DOPA 80:20 → Group 1\n• DOPC/DOPG 65:35 → Group 2\n• DOPC/DOPG 95:5 → Group 3'
    props = dict(boxstyle='round,pad=0.5', facecolor='#fef9e7', edgecolor='#f1c40f', linewidth=2)
    ax.text(0.98, 0.02, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    plt.savefig('CG_KA_Publication_Quality.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('CG_KA_Publication_Quality.pdf', facecolor='white', edgecolor='none')
    print("Saved: CG_KA_Publication_Quality.png/pdf")
    plt.show()
    plt.close()


# =============================================================================
# PLOT: Grouped Comparison - High Quality
# =============================================================================

def create_grouped_comparison():
    """
    Create high-quality grouped bar plot.
    """
    fig, ax = plt.subplots(figsize=(14, 7))
    
    ratios = list(dopa_systems.keys())
    dopa_values = list(dopa_systems.values())
    dopg_values = list(dopg_systems.values())
    
    x = np.arange(len(ratios))
    width = 0.38
    
    # Create bars
    bars1 = ax.bar(x - width/2, dopa_values, width, label='DOPC/DOPA', 
                   color=DOPA_COLOR, edgecolor='black', linewidth=1.5)
    bars2 = ax.bar(x + width/2, dopg_values, width, label='DOPC/DOPG', 
                   color=DOPG_COLOR, edgecolor='black', linewidth=1.5)
    
    # Highlight selected systems
    # DOPA 80:20 (index 4)
    bars1[4].set_edgecolor('#f1c40f')
    bars1[4].set_linewidth(4)
    bars1[4].set_color(DOPA_SELECTED_COLOR)
    
    # DOPG 65:35 (index 1) and 95:5 (index 7)
    bars2[1].set_edgecolor('#f1c40f')
    bars2[1].set_linewidth(4)
    bars2[1].set_color(DOPG_SELECTED_COLOR)
    
    bars2[7].set_edgecolor('#f1c40f')
    bars2[7].set_linewidth(4)
    bars2[7].set_color(DOPG_SELECTED_COLOR)
    
    # Add pure DOPC reference line
    ax.axhline(y=pure_dopc_ka, color='black', linestyle='--', linewidth=2.5, 
               label=f'Pure DOPC ({pure_dopc_ka:.1f} mN/m)')
    
    # Add stars for selected
    ax.annotate('★', xy=(4 - width/2, dopa_values[4] + 2), ha='center', fontsize=14, color='#f1c40f')
    ax.annotate('★', xy=(1 + width/2, dopg_values[1] + 2), ha='center', fontsize=14, color='#f1c40f')
    ax.annotate('★', xy=(7 + width/2, dopg_values[7] + 2), ha='center', fontsize=14, color='#f1c40f')
    
    # Customize
    ax.set_xlabel('PC : Secondary Lipid Ratio', fontsize=14, fontweight='bold')
    ax.set_ylabel('Area Compressibility Modulus, $K_A$ (mN/m)', fontsize=14, fontweight='bold')
    ax.set_title('CG Simulations: $K_A$ Comparison of DOPC/DOPA vs DOPC/DOPG',
                 fontweight='bold', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(ratios, fontsize=12)
    ax.legend(loc='upper right', fontsize=11, framealpha=0.95, edgecolor='black')
    ax.yaxis.grid(True, alpha=0.3, linestyle='--', zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylim(165, 215)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('CG_KA_Grouped_Publication.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('CG_KA_Grouped_Publication.pdf', facecolor='white', edgecolor='none')
    print("Saved: CG_KA_Grouped_Publication.png/pdf")
    plt.show()
    plt.close()


# =============================================================================
# PLOT: Line Plot - High Quality
# =============================================================================

def create_line_plot():
    """
    Create high-quality line plot showing K_A vs secondary lipid concentration.
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    ratios = ['60:40', '65:35', '70:30', '75:25', '80:20', '85:15', '90:10', '95:5']
    secondary_percent = [40, 35, 30, 25, 20, 15, 10, 5]
    
    dopa_values = [dopa_systems[r] for r in ratios]
    dopg_values = [dopg_systems[r] for r in ratios]
    
    # Plot lines
    ax.plot(secondary_percent, dopa_values, 'o-', color=DOPA_COLOR, linewidth=2.5, 
            markersize=10, markeredgecolor='black', markeredgewidth=1.5, label='DOPC/DOPA')
    ax.plot(secondary_percent, dopg_values, 's-', color=DOPG_COLOR, linewidth=2.5, 
            markersize=10, markeredgecolor='black', markeredgewidth=1.5, label='DOPC/DOPG')
    
    # Highlight selected points
    # DOPA 80:20 (20%)
    ax.plot(20, dopa_systems['80:20'], 'o', color=DOPA_SELECTED_COLOR, markersize=18,
            markeredgecolor='#f1c40f', markeredgewidth=4, zorder=10)
    # DOPG 65:35 (35%)
    ax.plot(35, dopg_systems['65:35'], 's', color=DOPG_SELECTED_COLOR, markersize=18,
            markeredgecolor='#f1c40f', markeredgewidth=4, zorder=10)
    # DOPG 95:5 (5%)
    ax.plot(5, dopg_systems['95:5'], 's', color=DOPG_SELECTED_COLOR, markersize=18,
            markeredgecolor='#f1c40f', markeredgewidth=4, zorder=10)
    
    # Add pure DOPC reference line
    ax.axhline(y=pure_dopc_ka, color='black', linestyle='--', linewidth=2.5, 
               label=f'Pure DOPC ({pure_dopc_ka:.1f} mN/m)')
    
    # Add annotations for selected
    ax.annotate('Group 1\n(80:20 PA)', xy=(20, dopa_systems['80:20']), 
                xytext=(25, 170), fontsize=10, ha='left',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef9e7', edgecolor='#f1c40f'))
    
    ax.annotate('Group 2\n(65:35 PG)', xy=(35, dopg_systems['65:35']), 
                xytext=(38, 172), fontsize=10, ha='left',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef9e7', edgecolor='#f1c40f'))
    
    ax.annotate('Group 3\n(95:5 PG)', xy=(5, dopg_systems['95:5']), 
                xytext=(8, 172), fontsize=10, ha='left',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5),
                bbox=dict(boxstyle='round,pad=0.3', facecolor='#fef9e7', edgecolor='#f1c40f'))
    
    # Customize
    ax.set_xlabel('Secondary Lipid Concentration (%)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Area Compressibility Modulus, $K_A$ (mN/m)', fontsize=14, fontweight='bold')
    ax.set_title('CG Simulations: $K_A$ vs Secondary Lipid Concentration',
                 fontweight='bold', fontsize=16)
    ax.legend(loc='upper center', fontsize=11, framealpha=0.95, edgecolor='black', ncol=3)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim(0, 45)
    ax.set_ylim(165, 210)
    ax.invert_xaxis()
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('CG_KA_Line_Publication.png', dpi=300, facecolor='white', edgecolor='none')
    plt.savefig('CG_KA_Line_Publication.pdf', facecolor='white', edgecolor='none')
    print("Saved: CG_KA_Line_Publication.png/pdf")
    plt.show()
    plt.close()


# =============================================================================
# PRINT SUMMARY
# =============================================================================

def print_summary():
    print("=" * 70)
    print("CG SIMULATIONS: K_A SUMMARY")
    print("=" * 70)
    print(f"\nPure DOPC Reference: {pure_dopc_ka:.2f} mN/m")
    
    print("\n" + "-" * 70)
    print("SELECTED SYSTEMS FOR ALL-ATOM SIMULATIONS:")
    print("-" * 70)
    print(f"  • DOPC/DOPA 80:20 → K_A = {dopa_systems['80:20']:.2f} mN/m (ΔK_A = {dopa_systems['80:20'] - pure_dopc_ka:+.2f})")
    print(f"  • DOPC/DOPG 65:35 → K_A = {dopg_systems['65:35']:.2f} mN/m (ΔK_A = {dopg_systems['65:35'] - pure_dopc_ka:+.2f})")
    print(f"  • DOPC/DOPG 95:5  → K_A = {dopg_systems['95:5']:.2f} mN/m (ΔK_A = {dopg_systems['95:5'] - pure_dopc_ka:+.2f})")
    print("=" * 70)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print_summary()
    
    print("\n[1/3] Creating publication-quality bar plot...")
    create_publication_plot()
    
    print("\n[2/3] Creating grouped comparison plot...")
    create_grouped_comparison()
    
    print("\n[3/3] Creating line plot...")
    create_line_plot()
    
    print("\n" + "=" * 70)
    print("All publication-quality plots generated!")
    print("=" * 70)
