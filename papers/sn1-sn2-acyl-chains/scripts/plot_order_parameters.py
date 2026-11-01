#!/usr/bin/env python3
"""
Order Parameter (S_CD) Comparison Plots
========================================
This script generates publication-quality figures comparing order parameters across:
1. Group 1: PC/PA systems (80:20 ratio)
2. Group 2: PC/PG systems (65:35 ratio)
3. Group 3: PC/PG systems (95:5 ratio)
4. All 27 systems combined
5. Chain type comparisons (18:1, 16:0, 18:0)

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
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
})

# =============================================================================
# DATA: Order Parameters (S_CD) for All 27 Systems
# =============================================================================

# Carbon positions
carbons_18 = np.arange(2, 17)  # C2-C16 for 18-carbon chains
carbons_16 = np.arange(2, 15)  # C2-C14 for 16-carbon chains (palmitoyl terminates at C16)

# -----------------------------------------------------------------------------
# GROUP 1: PC/PA Systems (80:20)
# -----------------------------------------------------------------------------

group1_data = {
    'name': 'Group 1: PC/PA (80:20)',
    'systems': {
        # DOPC/DOPA - all 18:1 chains
        'DOPC_DOPA': {
            'DOPC_sn2': [0.195, 0.186, 0.193, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'DOPA_sn2': [0.192, 0.183, 0.193, 0.168, 0.162, 0.095, 0.081, 0.078, 0.074, 0.099, 0.099, 0.101, 0.087, 0.074, 0.056],
            'DOPA_sn1': [0.143, 0.167, 0.165, 0.160, 0.145, 0.103, 0.071, 0.070, 0.066, 0.089, 0.087, 0.089, 0.074, 0.065, 0.050],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'DOPA_sn2': '18:1', 'DOPA_sn1': '18:1'},
        },
        # DOPC/POPA - DOPC 18:1/18:1, POPA 18:1/16:0
        'DOPC_POPA': {
            'DOPC_sn2': [0.195, 0.186, 0.193, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'POPA_sn2': [0.189, 0.180, 0.188, 0.166, 0.159, 0.096, 0.082, 0.079, 0.074, 0.099, 0.100, 0.102, 0.089, 0.078, 0.058],
            'POPA_sn1_16': [0.148, 0.180, 0.178, 0.189, 0.180, 0.177, 0.161, 0.152, 0.133, 0.119, 0.100, 0.083, 0.062],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'POPA_sn2': '18:1', 'POPA_sn1_16': '16:0'},
        },
        # DOPC/SOPA - DOPC 18:1/18:1, SOPA 18:1/18:0
        'DOPC_SOPA': {
            'DOPC_sn2': [0.194, 0.186, 0.192, 0.169, 0.163, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'DOPC_sn1': [0.144, 0.167, 0.165, 0.163, 0.148, 0.107, 0.073, 0.071, 0.067, 0.091, 0.088, 0.090, 0.075, 0.067, 0.050],
            'SOPA_sn2': [0.186, 0.178, 0.186, 0.165, 0.158, 0.097, 0.083, 0.080, 0.075, 0.100, 0.101, 0.103, 0.090, 0.079, 0.059],
            'SOPA_sn1': [0.147, 0.182, 0.181, 0.193, 0.185, 0.184, 0.170, 0.161, 0.143, 0.130, 0.113, 0.099, 0.083, 0.070, 0.052],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'SOPA_sn2': '18:1', 'SOPA_sn1': '18:0'},
        },
        # POPC/DOPA - POPC 18:1/16:0, DOPA 18:1/18:1
        'POPC_DOPA': {
            'POPC_sn2': [0.197, 0.188, 0.193, 0.169, 0.162, 0.097, 0.083, 0.080, 0.075, 0.101, 0.102, 0.104, 0.091, 0.081, 0.061],
            'POPC_sn1_16': [0.151, 0.185, 0.183, 0.194, 0.184, 0.181, 0.165, 0.157, 0.138, 0.124, 0.104, 0.087, 0.065],
            'DOPA_sn2': [0.185, 0.177, 0.184, 0.162, 0.156, 0.094, 0.081, 0.078, 0.073, 0.097, 0.098, 0.100, 0.087, 0.076, 0.057],
            'DOPA_sn1': [0.140, 0.163, 0.161, 0.157, 0.143, 0.104, 0.071, 0.070, 0.065, 0.088, 0.086, 0.088, 0.074, 0.065, 0.049],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'DOPA_sn2': '18:1', 'DOPA_sn1': '18:1'},
        },
        # POPC/POPA - POPC 18:1/16:0, POPA 18:1/16:0
        'POPC_POPA': {
            'POPC_sn2': [0.197, 0.187, 0.192, 0.168, 0.161, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.060],
            'POPC_sn1_16': [0.151, 0.184, 0.182, 0.193, 0.183, 0.180, 0.164, 0.156, 0.137, 0.123, 0.103, 0.086, 0.064],
            'POPA_sn2': [0.183, 0.175, 0.182, 0.160, 0.154, 0.093, 0.080, 0.078, 0.073, 0.097, 0.098, 0.100, 0.087, 0.077, 0.058],
            'POPA_sn1_16': [0.145, 0.177, 0.175, 0.186, 0.177, 0.174, 0.159, 0.150, 0.132, 0.118, 0.099, 0.083, 0.062],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'POPA_sn2': '18:1', 'POPA_sn1_16': '16:0'},
        },
        # POPC/SOPA - POPC 18:1/16:0, SOPA 18:1/18:0
        'POPC_SOPA': {
            'POPC_sn2': [0.197, 0.187, 0.192, 0.169, 0.162, 0.097, 0.083, 0.080, 0.075, 0.101, 0.102, 0.104, 0.091, 0.081, 0.061],
            'POPC_sn1_16': [0.152, 0.185, 0.183, 0.194, 0.184, 0.181, 0.165, 0.157, 0.138, 0.124, 0.104, 0.087, 0.065],
            'SOPA_sn2': [0.183, 0.175, 0.183, 0.161, 0.155, 0.094, 0.081, 0.079, 0.074, 0.098, 0.099, 0.101, 0.088, 0.078, 0.059],
            'SOPA_sn1': [0.145, 0.179, 0.178, 0.190, 0.182, 0.181, 0.167, 0.158, 0.141, 0.128, 0.111, 0.097, 0.082, 0.069, 0.052],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'SOPA_sn2': '18:1', 'SOPA_sn1': '18:0'},
        },
        # SOPC/DOPA - SOPC 18:1/18:0, DOPA 18:1/18:1
        'SOPC_DOPA': {
            'SOPC_sn2': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'SOPC_sn1': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'DOPA_sn2': [0.179, 0.175, 0.182, 0.161, 0.155, 0.094, 0.082, 0.080, 0.077, 0.103, 0.103, 0.106, 0.096, 0.085, 0.064],
            'DOPA_sn1': [0.139, 0.164, 0.160, 0.157, 0.144, 0.106, 0.074, 0.068, 0.066, 0.090, 0.090, 0.095, 0.081, 0.073, 0.055],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'DOPA_sn2': '18:1', 'DOPA_sn1': '18:1'},
        },
        # SOPC/POPA - SOPC 18:1/18:0, POPA 18:1/16:0
        'SOPC_POPA': {
            'SOPC_sn2': [0.200, 0.190, 0.195, 0.172, 0.165, 0.099, 0.086, 0.082, 0.079, 0.107, 0.109, 0.112, 0.100, 0.090, 0.069],
            'SOPC_sn1': [0.151, 0.185, 0.184, 0.196, 0.188, 0.188, 0.174, 0.165, 0.148, 0.136, 0.117, 0.104, 0.087, 0.073, 0.054],
            'POPA_sn2': [0.180, 0.173, 0.181, 0.160, 0.154, 0.094, 0.082, 0.079, 0.075, 0.101, 0.102, 0.105, 0.093, 0.083, 0.063],
            'POPA_sn1_16': [0.143, 0.175, 0.173, 0.184, 0.175, 0.173, 0.158, 0.149, 0.131, 0.118, 0.100, 0.084, 0.063],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'POPA_sn2': '18:1', 'POPA_sn1_16': '16:0'},
        },
        # SOPC/SOPA - SOPC 18:1/18:0, SOPA 18:1/18:0
        'SOPC_SOPA': {
            'SOPC_sn2': [0.200, 0.191, 0.196, 0.173, 0.166, 0.100, 0.086, 0.083, 0.078, 0.106, 0.108, 0.111, 0.098, 0.088, 0.068],
            'SOPC_sn1': [0.152, 0.186, 0.186, 0.198, 0.190, 0.188, 0.175, 0.166, 0.148, 0.135, 0.117, 0.103, 0.086, 0.072, 0.054],
            'SOPA_sn2': [0.187, 0.181, 0.189, 0.170, 0.163, 0.099, 0.085, 0.080, 0.075, 0.101, 0.102, 0.103, 0.089, 0.078, 0.058],
            'SOPA_sn1': [0.145, 0.180, 0.180, 0.193, 0.186, 0.184, 0.170, 0.159, 0.140, 0.123, 0.106, 0.091, 0.074, 0.062, 0.046],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'SOPA_sn2': '18:1', 'SOPA_sn1': '18:0'},
        },
    }
}

# -----------------------------------------------------------------------------
# GROUP 2: PC/PG Systems (65:35)
# -----------------------------------------------------------------------------

group2_data = {
    'name': 'Group 2: PC/PG (65:35)',
    'systems': {
        # DOPC/DOPG - all 18:1 chains
        'DOPC_DOPG': {
            'DOPC_sn2': [0.194, 0.185, 0.192, 0.170, 0.163, 0.098, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'DOPG_sn2': [0.191, 0.182, 0.190, 0.168, 0.161, 0.096, 0.081, 0.078, 0.074, 0.099, 0.099, 0.101, 0.087, 0.076, 0.057],
            'DOPG_sn1': [0.142, 0.165, 0.163, 0.161, 0.146, 0.105, 0.072, 0.070, 0.066, 0.089, 0.087, 0.089, 0.075, 0.066, 0.049],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # DOPC/POPG
        'DOPC_POPG': {
            'DOPC_sn2': [0.194, 0.186, 0.192, 0.170, 0.163, 0.097, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'POPG_sn2': [0.188, 0.180, 0.187, 0.165, 0.158, 0.095, 0.081, 0.079, 0.074, 0.099, 0.100, 0.102, 0.089, 0.078, 0.059],
            'POPG_sn1_16': [0.147, 0.179, 0.177, 0.188, 0.179, 0.176, 0.161, 0.152, 0.134, 0.120, 0.101, 0.085, 0.063],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # DOPC/SOPG
        'DOPC_SOPG': {
            'DOPC_sn2': [0.194, 0.185, 0.191, 0.169, 0.162, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'SOPG_sn2': [0.186, 0.178, 0.185, 0.163, 0.156, 0.094, 0.081, 0.078, 0.073, 0.098, 0.099, 0.101, 0.088, 0.078, 0.058],
            'SOPG_sn1': [0.147, 0.183, 0.182, 0.194, 0.186, 0.185, 0.171, 0.162, 0.145, 0.132, 0.115, 0.101, 0.085, 0.072, 0.054],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
        # POPC/DOPG
        'POPC_DOPG': {
            'POPC_sn2': [0.197, 0.188, 0.193, 0.169, 0.162, 0.097, 0.084, 0.080, 0.076, 0.102, 0.103, 0.105, 0.092, 0.082, 0.062],
            'POPC_sn1_16': [0.152, 0.186, 0.184, 0.195, 0.185, 0.182, 0.166, 0.158, 0.139, 0.125, 0.105, 0.088, 0.066],
            'DOPG_sn2': [0.184, 0.176, 0.183, 0.162, 0.155, 0.094, 0.081, 0.078, 0.073, 0.097, 0.098, 0.100, 0.087, 0.077, 0.058],
            'DOPG_sn1': [0.139, 0.162, 0.160, 0.157, 0.143, 0.103, 0.071, 0.069, 0.065, 0.088, 0.086, 0.088, 0.074, 0.066, 0.049],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # POPC/POPG
        'POPC_POPG': {
            'POPC_sn2': [0.198, 0.189, 0.194, 0.170, 0.163, 0.098, 0.084, 0.080, 0.076, 0.102, 0.103, 0.106, 0.093, 0.083, 0.063],
            'POPC_sn1_16': [0.153, 0.187, 0.185, 0.197, 0.187, 0.184, 0.168, 0.163, 0.143, 0.128, 0.108, 0.090, 0.067],
            'POPG_sn2': [0.185, 0.177, 0.184, 0.162, 0.156, 0.094, 0.081, 0.077, 0.073, 0.097, 0.098, 0.100, 0.088, 0.078, 0.059],
            'POPG_sn1_16': [0.146, 0.178, 0.176, 0.186, 0.177, 0.174, 0.159, 0.147, 0.129, 0.116, 0.098, 0.082, 0.061],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # POPC/SOPG
        'POPC_SOPG': {
            'POPC_sn2': [0.198, 0.189, 0.194, 0.170, 0.164, 0.098, 0.084, 0.080, 0.076, 0.102, 0.104, 0.106, 0.093, 0.083, 0.063],
            'POPC_sn1_16': [0.154, 0.188, 0.186, 0.198, 0.188, 0.186, 0.170, 0.166, 0.147, 0.132, 0.111, 0.092, 0.068],
            'SOPG_sn2': [0.183, 0.175, 0.183, 0.162, 0.156, 0.095, 0.082, 0.080, 0.075, 0.100, 0.101, 0.103, 0.091, 0.081, 0.061],
            'SOPG_sn1': [0.146, 0.181, 0.180, 0.192, 0.184, 0.183, 0.169, 0.166, 0.149, 0.137, 0.119, 0.105, 0.088, 0.074, 0.055],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
        # SOPC/DOPG
        'SOPC_DOPG': {
            'SOPC_sn2': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'SOPC_sn1': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'DOPG_sn2': [0.182, 0.175, 0.182, 0.163, 0.156, 0.095, 0.083, 0.080, 0.076, 0.102, 0.103, 0.106, 0.094, 0.084, 0.063],
            'DOPG_sn1': [0.139, 0.164, 0.161, 0.158, 0.145, 0.107, 0.075, 0.074, 0.070, 0.094, 0.092, 0.095, 0.082, 0.074, 0.056],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # SOPC/POPG
        'SOPC_POPG': {
            'SOPC_sn2': [0.200, 0.191, 0.196, 0.174, 0.167, 0.100, 0.087, 0.083, 0.079, 0.107, 0.109, 0.112, 0.100, 0.090, 0.069],
            'SOPC_sn1': [0.152, 0.186, 0.186, 0.199, 0.190, 0.189, 0.176, 0.170, 0.152, 0.140, 0.121, 0.107, 0.090, 0.076, 0.056],
            'POPG_sn2': [0.185, 0.177, 0.184, 0.164, 0.158, 0.096, 0.083, 0.082, 0.077, 0.103, 0.104, 0.107, 0.095, 0.085, 0.064],
            'POPG_sn1_16': [0.148, 0.181, 0.179, 0.190, 0.181, 0.179, 0.164, 0.155, 0.137, 0.123, 0.104, 0.087, 0.065],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # SOPC/SOPG
        'SOPC_SOPG': {
            'SOPC_sn2': [0.200, 0.190, 0.196, 0.172, 0.166, 0.100, 0.086, 0.083, 0.079, 0.107, 0.108, 0.114, 0.100, 0.090, 0.070],
            'SOPC_sn1': [0.152, 0.186, 0.186, 0.197, 0.189, 0.188, 0.175, 0.168, 0.149, 0.137, 0.118, 0.109, 0.088, 0.074, 0.058],
            'SOPG_sn2': [0.192, 0.182, 0.190, 0.168, 0.162, 0.099, 0.085, 0.082, 0.075, 0.102, 0.102, 0.103, 0.089, 0.077, 0.058],
            'SOPG_sn1': [0.146, 0.183, 0.184, 0.193, 0.184, 0.183, 0.169, 0.159, 0.140, 0.122, 0.104, 0.091, 0.073, 0.061, 0.045],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
    }
}

# -----------------------------------------------------------------------------
# GROUP 3: PC/PG Systems (95:5)
# -----------------------------------------------------------------------------

group3_data = {
    'name': 'Group 3: PC/PG (95:5)',
    'systems': {
        # DOPC/DOPG
        'DOPC_DOPG': {
            'DOPC_sn2': [0.194, 0.185, 0.192, 0.170, 0.163, 0.097, 0.082, 0.080, 0.076, 0.101, 0.101, 0.103, 0.088, 0.077, 0.059],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'DOPG_sn2': [0.189, 0.180, 0.188, 0.166, 0.159, 0.095, 0.080, 0.076, 0.072, 0.097, 0.097, 0.099, 0.086, 0.075, 0.056],
            'DOPG_sn1': [0.140, 0.163, 0.161, 0.158, 0.144, 0.104, 0.071, 0.070, 0.065, 0.088, 0.086, 0.088, 0.074, 0.065, 0.049],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # DOPC/POPG
        'DOPC_POPG': {
            'DOPC_sn2': [0.194, 0.186, 0.192, 0.170, 0.163, 0.097, 0.081, 0.081, 0.076, 0.102, 0.102, 0.104, 0.090, 0.080, 0.061],
            'DOPC_sn1': [0.144, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.071, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'POPG_sn2': [0.193, 0.186, 0.191, 0.171, 0.163, 0.098, 0.084, 0.078, 0.074, 0.102, 0.105, 0.107, 0.096, 0.086, 0.067],
            'POPG_sn1_16': [0.150, 0.183, 0.185, 0.193, 0.184, 0.182, 0.167, 0.156, 0.138, 0.125, 0.106, 0.091, 0.069],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # DOPC/SOPG
        'DOPC_SOPG': {
            'DOPC_sn2': [0.194, 0.185, 0.191, 0.169, 0.162, 0.097, 0.081, 0.079, 0.075, 0.100, 0.100, 0.102, 0.088, 0.077, 0.058],
            'DOPC_sn1': [0.145, 0.168, 0.166, 0.164, 0.149, 0.108, 0.074, 0.072, 0.068, 0.092, 0.089, 0.091, 0.076, 0.068, 0.051],
            'SOPG_sn2': [0.186, 0.182, 0.189, 0.173, 0.162, 0.101, 0.086, 0.081, 0.074, 0.102, 0.105, 0.107, 0.096, 0.085, 0.065],
            'SOPG_sn1': [0.147, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'chains': {'DOPC_sn2': '18:1', 'DOPC_sn1': '18:1', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
        # POPC/DOPG
        'POPC_DOPG': {
            'POPC_sn2': [0.197, 0.187, 0.191, 0.167, 0.160, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.061],
            'POPC_sn1_16': [0.151, 0.184, 0.182, 0.192, 0.182, 0.179, 0.163, 0.153, 0.134, 0.120, 0.101, 0.085, 0.063],
            'DOPG_sn2': [0.182, 0.173, 0.179, 0.164, 0.155, 0.093, 0.080, 0.078, 0.073, 0.098, 0.098, 0.100, 0.088, 0.078, 0.055],
            'DOPG_sn1': [0.137, 0.159, 0.156, 0.151, 0.138, 0.099, 0.068, 0.067, 0.063, 0.085, 0.083, 0.086, 0.072, 0.064, 0.049],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # POPC/POPG
        'POPC_POPG': {
            'POPC_sn2': [0.197, 0.186, 0.191, 0.167, 0.160, 0.096, 0.083, 0.079, 0.074, 0.100, 0.101, 0.103, 0.090, 0.080, 0.061],
            'POPC_sn1_16': [0.151, 0.184, 0.182, 0.193, 0.183, 0.180, 0.164, 0.154, 0.135, 0.121, 0.101, 0.085, 0.063],
            'POPG_sn2': [0.182, 0.174, 0.179, 0.160, 0.153, 0.093, 0.081, 0.077, 0.072, 0.096, 0.097, 0.098, 0.087, 0.076, 0.058],
            'POPG_sn1_16': [0.143, 0.175, 0.173, 0.183, 0.174, 0.171, 0.156, 0.147, 0.129, 0.115, 0.097, 0.081, 0.061],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # POPC/SOPG
        'POPC_SOPG': {
            'POPC_sn2': [0.197, 0.187, 0.192, 0.168, 0.161, 0.096, 0.083, 0.079, 0.075, 0.101, 0.102, 0.103, 0.091, 0.081, 0.061],
            'POPC_sn1_16': [0.152, 0.185, 0.184, 0.194, 0.184, 0.181, 0.165, 0.156, 0.138, 0.124, 0.104, 0.085, 0.061],
            'SOPG_sn2': [0.183, 0.175, 0.183, 0.162, 0.155, 0.094, 0.081, 0.079, 0.074, 0.099, 0.099, 0.101, 0.088, 0.078, 0.059],
            'SOPG_sn1': [0.139, 0.182, 0.180, 0.191, 0.183, 0.182, 0.167, 0.156, 0.138, 0.125, 0.108, 0.089, 0.074, 0.061, 0.044],
            'chains': {'POPC_sn2': '18:1', 'POPC_sn1_16': '16:0', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
        # SOPC/DOPG
        'SOPC_DOPG': {
            'SOPC_sn2': [0.201, 0.191, 0.196, 0.173, 0.166, 0.100, 0.087, 0.083, 0.080, 0.108, 0.110, 0.113, 0.101, 0.091, 0.070],
            'SOPC_sn1': [0.151, 0.185, 0.185, 0.197, 0.189, 0.189, 0.175, 0.166, 0.149, 0.137, 0.118, 0.105, 0.088, 0.074, 0.055],
            'DOPG_sn2': [0.179, 0.175, 0.182, 0.161, 0.155, 0.094, 0.082, 0.080, 0.077, 0.103, 0.103, 0.106, 0.096, 0.085, 0.064],
            'DOPG_sn1': [0.139, 0.164, 0.160, 0.157, 0.144, 0.106, 0.074, 0.068, 0.066, 0.090, 0.090, 0.095, 0.081, 0.073, 0.055],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'DOPG_sn2': '18:1', 'DOPG_sn1': '18:1'},
        },
        # SOPC/POPG
        'SOPC_POPG': {
            'SOPC_sn2': [0.201, 0.191, 0.196, 0.172, 0.165, 0.099, 0.087, 0.083, 0.080, 0.108, 0.111, 0.113, 0.102, 0.091, 0.070],
            'SOPC_sn1': [0.150, 0.184, 0.184, 0.196, 0.189, 0.189, 0.175, 0.167, 0.150, 0.138, 0.120, 0.107, 0.089, 0.074, 0.055],
            'POPG_sn2': [0.193, 0.186, 0.191, 0.171, 0.163, 0.098, 0.084, 0.078, 0.074, 0.102, 0.105, 0.107, 0.096, 0.086, 0.067],
            'POPG_sn1_16': [0.150, 0.183, 0.185, 0.193, 0.184, 0.182, 0.167, 0.156, 0.138, 0.125, 0.106, 0.091, 0.069],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'POPG_sn2': '18:1', 'POPG_sn1_16': '16:0'},
        },
        # SOPC/SOPG
        'SOPC_SOPG': {
            'SOPC_sn2': [0.200, 0.190, 0.196, 0.172, 0.165, 0.100, 0.088, 0.084, 0.081, 0.109, 0.112, 0.115, 0.103, 0.092, 0.071],
            'SOPC_sn1': [0.151, 0.186, 0.185, 0.197, 0.189, 0.189, 0.176, 0.167, 0.150, 0.139, 0.120, 0.107, 0.089, 0.075, 0.056],
            'SOPG_sn2': [0.187, 0.182, 0.189, 0.170, 0.162, 0.101, 0.086, 0.078, 0.074, 0.102, 0.105, 0.108, 0.096, 0.085, 0.065],
            'SOPG_sn1': [0.144, 0.180, 0.180, 0.191, 0.184, 0.184, 0.171, 0.162, 0.144, 0.131, 0.114, 0.101, 0.085, 0.072, 0.055],
            'chains': {'SOPC_sn2': '18:1', 'SOPC_sn1': '18:0', 'SOPG_sn2': '18:1', 'SOPG_sn1': '18:0'},
        },
    }
}

# Colors
PC_COLORS = {
    'DOPC': '#3498db',
    'POPC': '#2ecc71', 
    'SOPC': '#9b59b6',
}

CHAIN_COLORS = {
    '18:1': '#e74c3c',  # Red for unsaturated
    '16:0': '#3498db',  # Blue for saturated 16
    '18:0': '#2c3e50',  # Dark for saturated 18
}


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_group_sn2_comparison(group_data, filename_prefix):
    """
    Plot sn-2 chain order parameters for all 9 systems in a group.
    """
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True, sharey=True)
    
    systems = list(group_data['systems'].keys())
    
    for idx, (ax, sys_name) in enumerate(zip(axes.flatten(), systems)):
        sys_data = group_data['systems'][sys_name]
        
        # Get PC and anionic lipid names
        pc_name, anionic_name = sys_name.split('_')
        
        # Plot PC sn-2
        pc_sn2_key = f'{pc_name}_sn2'
        ax.plot(carbons_18, sys_data[pc_sn2_key], 'o-', color=PC_COLORS[pc_name],
                linewidth=2, markersize=5, label=f'{pc_name} sn-2')
        
        # Plot anionic sn-2
        anionic_sn2_key = f'{anionic_name}_sn2'
        ax.plot(carbons_18, sys_data[anionic_sn2_key], 's--', color='#e74c3c',
                linewidth=2, markersize=5, label=f'{anionic_name} sn-2', alpha=0.8)
        
        # Highlight double bond region
        ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        
        ax.set_title(f'{pc_name}/{anionic_name}', fontsize=11, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0, 0.22)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8)
    
    # Common labels
    fig.text(0.5, 0.02, 'Carbon Position', ha='center', fontsize=13, fontweight='bold')
    fig.text(0.02, 0.5, '$|S_{CD}|$', va='center', rotation='vertical', fontsize=13, fontweight='bold')
    fig.suptitle(f'{group_data["name"]}\nsn-2 Chain Order Parameters (All 18:1)', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.06, bottom=0.06, top=0.92)
    plt.savefig(f'{filename_prefix}_sn2.png', dpi=300, facecolor='white')
    plt.savefig(f'{filename_prefix}_sn2.pdf', facecolor='white')
    print(f"Saved: {filename_prefix}_sn2.png/pdf")
    plt.show()


def plot_group_sn1_comparison(group_data, filename_prefix):
    """
    Plot sn-1 chain order parameters for all 9 systems in a group.
    """
    fig, axes = plt.subplots(3, 3, figsize=(15, 12), sharex=True, sharey=True)
    
    systems = list(group_data['systems'].keys())
    
    for idx, (ax, sys_name) in enumerate(zip(axes.flatten(), systems)):
        sys_data = group_data['systems'][sys_name]
        chains_info = sys_data['chains']
        
        # Get PC and anionic lipid names
        pc_name, anionic_name = sys_name.split('_')
        
        # Find and plot PC sn-1
        pc_sn1_keys = [k for k in sys_data.keys() if k.startswith(pc_name) and 'sn1' in k]
        for key in pc_sn1_keys:
            chain_type = chains_info[key]
            carbons = carbons_16 if '16' in key else carbons_18
            label = f'{pc_name} sn-1 ({chain_type})'
            ax.plot(carbons, sys_data[key], 'o-', color=PC_COLORS[pc_name],
                    linewidth=2, markersize=5, label=label)
        
        # Find and plot anionic sn-1
        anionic_sn1_keys = [k for k in sys_data.keys() if k.startswith(anionic_name) and 'sn1' in k]
        for key in anionic_sn1_keys:
            chain_type = chains_info[key]
            carbons = carbons_16 if '16' in key else carbons_18
            label = f'{anionic_name} sn-1 ({chain_type})'
            ax.plot(carbons, sys_data[key], 's--', color='#e74c3c',
                    linewidth=2, markersize=5, label=label, alpha=0.8)
        
        # Highlight double bond region (only for 18:1 chains)
        ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        
        ax.set_title(f'{pc_name}/{anionic_name}', fontsize=11, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0, 0.22)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=7)
    
    # Common labels
    fig.text(0.5, 0.02, 'Carbon Position', ha='center', fontsize=13, fontweight='bold')
    fig.text(0.02, 0.5, '$|S_{CD}|$', va='center', rotation='vertical', fontsize=13, fontweight='bold')
    fig.suptitle(f'{group_data["name"]}\nsn-1 Chain Order Parameters', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.06, bottom=0.06, top=0.92)
    plt.savefig(f'{filename_prefix}_sn1.png', dpi=300, facecolor='white')
    plt.savefig(f'{filename_prefix}_sn1.pdf', facecolor='white')
    print(f"Saved: {filename_prefix}_sn1.png/pdf")
    plt.show()


def plot_group_pc_comparison(group_data, filename_prefix):
    """
    Plot PC lipid sn-2 chains grouped by PC type (DOPC, POPC, SOPC).
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    pc_lipids = ['DOPC', 'POPC', 'SOPC']
    # Anionic lipid colors - consistent scheme:
    # DO- variants (DOPA, DOPG) = Red
    # PO- variants (POPA, POPG) = Orange  
    # SO- variants (SOPA, SOPG) = Dark gray
    anionic_colors = {'DOPA': '#e74c3c', 'POPA': '#f39c12', 'SOPA': '#7f8c8d',
                      'DOPG': '#e74c3c', 'POPG': '#f39c12', 'SOPG': '#7f8c8d'}
    
    for ax, pc in zip(axes, pc_lipids):
        # Find all systems with this PC lipid
        for sys_name, sys_data in group_data['systems'].items():
            if sys_name.startswith(pc):
                anionic = sys_name.split('_')[1]
                pc_sn2_key = f'{pc}_sn2'
                
                ax.plot(carbons_18, sys_data[pc_sn2_key], 'o-', 
                        color=anionic_colors.get(anionic, '#333'),
                        linewidth=2, markersize=5, label=f'{pc}/{anionic}')
        
        ax.axvspan(8.5, 10.5, alpha=0.15, color='yellow')
        ax.set_xlabel('Carbon Position', fontsize=11, fontweight='bold')
        ax.set_title(f'{pc} sn-2 (18:1)', fontsize=12, fontweight='bold')
        ax.set_xlim(1, 17)
        ax.set_ylim(0, 0.22)
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=9)
    
    axes[0].set_ylabel('$|S_{CD}|$', fontsize=12, fontweight='bold')
    
    fig.suptitle(f'{group_data["name"]}\nPC Lipid sn-2 Chain Comparison', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(f'{filename_prefix}_PC_comparison.png', dpi=300, facecolor='white')
    plt.savefig(f'{filename_prefix}_PC_comparison.pdf', facecolor='white')
    print(f"Saved: {filename_prefix}_PC_comparison.png/pdf")
    plt.show()


def plot_all_systems_chain_types():
    """
    Compare order parameters by chain type across all 27 systems.
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    
    # Collect all chains by type
    chains_18_1_sn2 = []
    chains_18_1_sn1 = []
    chains_16_0 = []
    chains_18_0 = []
    
    for group_data in [group1_data, group2_data, group3_data]:
        for sys_name, sys_data in group_data['systems'].items():
            chains_info = sys_data['chains']
            for key, chain_type in chains_info.items():
                if 'sn2' in key and chain_type == '18:1':
                    chains_18_1_sn2.append(sys_data[key])
                elif 'sn1' in key:
                    if chain_type == '18:1':
                        chains_18_1_sn1.append(sys_data[key])
                    elif chain_type == '16:0':
                        chains_16_0.append(sys_data[key])
                    elif chain_type == '18:0':
                        chains_18_0.append(sys_data[key])
    
    # Calculate averages
    avg_18_1_sn2 = np.mean(chains_18_1_sn2, axis=0)
    std_18_1_sn2 = np.std(chains_18_1_sn2, axis=0)
    
    avg_18_1_sn1 = np.mean(chains_18_1_sn1, axis=0)
    std_18_1_sn1 = np.std(chains_18_1_sn1, axis=0)
    
    avg_16_0 = np.mean(chains_16_0, axis=0)
    std_16_0 = np.std(chains_16_0, axis=0)
    
    avg_18_0 = np.mean(chains_18_0, axis=0)
    std_18_0 = np.std(chains_18_0, axis=0)
    
    # Plot 1: 18:1 chains
    ax1 = axes[0]
    ax1.errorbar(carbons_18, avg_18_1_sn2, yerr=std_18_1_sn2, fmt='o-', color='#e74c3c',
                 linewidth=2, markersize=6, capsize=3, label='18:1 sn-2')
    ax1.errorbar(carbons_18, avg_18_1_sn1, yerr=std_18_1_sn1, fmt='s--', color='#c0392b',
                 linewidth=2, markersize=6, capsize=3, label='18:1 sn-1', alpha=0.8)
    ax1.axvspan(8.5, 10.5, alpha=0.2, color='yellow', label='Double bond')
    ax1.set_xlabel('Carbon Position', fontsize=11, fontweight='bold')
    ax1.set_ylabel('$|S_{CD}|$', fontsize=11, fontweight='bold')
    ax1.set_title('18:1 (Oleoyl) Chains\nAll 27 Systems', fontsize=12, fontweight='bold')
    ax1.set_xlim(1, 17)
    ax1.set_ylim(0, 0.22)
    ax1.legend(loc='upper right', fontsize=9)
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Plot 2: 16:0 chains
    ax2 = axes[1]
    ax2.errorbar(carbons_16, avg_16_0, yerr=std_16_0, fmt='^-', color='#3498db',
                 linewidth=2, markersize=6, capsize=3, label='16:0 sn-1')
    ax2.set_xlabel('Carbon Position', fontsize=11, fontweight='bold')
    ax2.set_title('16:0 (Palmitoyl) Chains\nAll 27 Systems', fontsize=12, fontweight='bold')
    ax2.set_xlim(1, 17)
    ax2.set_ylim(0, 0.22)
    ax2.legend(loc='upper right', fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # Plot 3: 18:0 chains
    ax3 = axes[2]
    ax3.errorbar(carbons_18, avg_18_0, yerr=std_18_0, fmt='D-', color='#2c3e50',
                 linewidth=2, markersize=6, capsize=3, label='18:0 sn-1')
    ax3.set_xlabel('Carbon Position', fontsize=11, fontweight='bold')
    ax3.set_title('18:0 (Stearoyl) Chains\nAll 27 Systems', fontsize=12, fontweight='bold')
    ax3.set_xlim(1, 17)
    ax3.set_ylim(0, 0.22)
    ax3.legend(loc='upper right', fontsize=9)
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    fig.suptitle('Order Parameters by Chain Type - Average Across All 27 Systems', 
                 fontsize=14, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig('Order_All_Systems_Chain_Types.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_Systems_Chain_Types.pdf', facecolor='white')
    print("Saved: Order_All_Systems_Chain_Types.png/pdf")
    plt.show()


def plot_all_systems_comparison():
    """
    Create a comprehensive comparison of all 27 systems.
    """
    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    
    groups = [group1_data, group2_data, group3_data]
    group_names = ['PA (90:10)', 'PG (60:40)', 'PG (90:10)']
    pc_lipids = ['DOPC', 'POPC', 'SOPC']
    
    for row, (group_data, group_name) in enumerate(zip(groups, group_names)):
        for col, pc in enumerate(pc_lipids):
            ax = axes[row, col]
            
            # Find all systems with this PC lipid
            for sys_name, sys_data in group_data['systems'].items():
                if sys_name.startswith(pc):
                    anionic = sys_name.split('_')[1]
                    pc_sn2_key = f'{pc}_sn2'
                    
                    # Color based on anionic lipid type
                    # DO- variants = Red, PO- variants = Orange, SO- variants = Dark gray
                    if 'DO' in anionic:
                        color = '#e74c3c'  # Red
                        marker = 'o'
                    elif 'PO' in anionic:
                        color = '#f39c12'  # Orange
                        marker = 's'
                    else:  # SO
                        color = '#7f8c8d'  # Dark gray
                        marker = '^'
                    
                    ax.plot(carbons_18, sys_data[pc_sn2_key], f'{marker}-', 
                            color=color, linewidth=1.5, markersize=4, 
                            label=anionic, alpha=0.9)
            
            ax.axvspan(8.5, 10.5, alpha=0.1, color='yellow')
            ax.set_xlim(1, 17)
            ax.set_ylim(0.04, 0.22)
            ax.grid(True, linestyle='--', alpha=0.4)
            
            if row == 0:
                ax.set_title(f'{pc}', fontsize=12, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{group_name}\n$|S_{{CD}}|$', fontsize=10, fontweight='bold')
            if row == 2:
                ax.set_xlabel('Carbon Position', fontsize=10, fontweight='bold')
            
            ax.legend(loc='upper right', fontsize=7, ncol=1)
    
    fig.suptitle('Order Parameters (sn-2 Chains) - All 27 Systems\nRows: Groups | Columns: PC Lipids', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.savefig('Order_All_27_Systems_Grid.png', dpi=300, facecolor='white')
    plt.savefig('Order_All_27_Systems_Grid.pdf', facecolor='white')
    print("Saved: Order_All_27_Systems_Grid.png/pdf")
    plt.show()


def plot_chain_type_summary():
    """
    Create a summary plot comparing all chain types on one figure.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Collect and average chains
    chains_18_1_sn2 = []
    chains_18_1_sn1 = []
    chains_16_0 = []
    chains_18_0 = []
    
    for group_data in [group1_data, group2_data, group3_data]:
        for sys_name, sys_data in group_data['systems'].items():
            chains_info = sys_data['chains']
            for key, chain_type in chains_info.items():
                if 'sn2' in key and chain_type == '18:1':
                    chains_18_1_sn2.append(sys_data[key])
                elif 'sn1' in key:
                    if chain_type == '18:1':
                        chains_18_1_sn1.append(sys_data[key])
                    elif chain_type == '16:0':
                        chains_16_0.append(sys_data[key])
                    elif chain_type == '18:0':
                        chains_18_0.append(sys_data[key])
    
    # Calculate averages
    avg_18_1_sn2 = np.mean(chains_18_1_sn2, axis=0)
    avg_18_1_sn1 = np.mean(chains_18_1_sn1, axis=0)
    avg_16_0 = np.mean(chains_16_0, axis=0)
    avg_18_0 = np.mean(chains_18_0, axis=0)
    
    # Plot all chain types
    ax.plot(carbons_18, avg_18_1_sn2, 'o-', color='#e74c3c', linewidth=2.5, 
            markersize=8, label='18:1 sn-2 (oleoyl)')
    ax.plot(carbons_18, avg_18_1_sn1, 's--', color='#c0392b', linewidth=2.5, 
            markersize=8, label='18:1 sn-1 (oleoyl)', alpha=0.8)
    ax.plot(carbons_16, avg_16_0, '^-', color='#3498db', linewidth=2.5, 
            markersize=8, label='16:0 sn-1 (palmitoyl)')
    ax.plot(carbons_18, avg_18_0, 'D-', color='#2c3e50', linewidth=2.5, 
            markersize=8, label='18:0 sn-1 (stearoyl)')
    
    # Highlight double bond region
    ax.axvspan(8.5, 10.5, alpha=0.2, color='yellow', label='Double bond region')
    
    # Annotations
    ax.annotate('Saturated chains:\nmonotonic decrease', 
                xy=(6, 0.185), xytext=(3, 0.21),
                fontsize=10, ha='center',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    ax.annotate('Unsaturated chains:\ndouble bond dip', 
                xy=(9.5, 0.072), xytext=(13, 0.05),
                fontsize=10, ha='center',
                arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
    
    ax.set_xlabel('Carbon Position', fontsize=13, fontweight='bold')
    ax.set_ylabel('$|S_{CD}|$', fontsize=13, fontweight='bold')
    ax.set_title('Order Parameter Profiles by Chain Type\n(Averaged Across All 27 Systems)', 
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1, 17)
    ax.set_ylim(0, 0.23)
    ax.set_xticks(range(2, 17))
    ax.legend(loc='upper right', fontsize=11, framealpha=1)
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('Order_Chain_Type_Summary.png', dpi=300, facecolor='white')
    plt.savefig('Order_Chain_Type_Summary.pdf', facecolor='white')
    print("Saved: Order_Chain_Type_Summary.png/pdf")
    plt.show()


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Order Parameter (S_CD) Comparison Plots")
    print("=" * 70)
    print("\nGenerating plots...\n")
    
    # Group 1 plots
    print("-" * 50)
    print("GROUP 1: PC/PA Systems (80:20)")
    print("-" * 50)
    print("1a. Plotting Group 1 sn-2 chains...")
    plot_group_sn2_comparison(group1_data, 'Order_Group1_PA')
    print("1b. Plotting Group 1 sn-1 chains...")
    plot_group_sn1_comparison(group1_data, 'Order_Group1_PA')
    print("1c. Plotting Group 1 PC comparison...")
    plot_group_pc_comparison(group1_data, 'Order_Group1_PA')
    
    # Group 2 plots
    print("\n" + "-" * 50)
    print("GROUP 2: PC/PG Systems (65:35)")
    print("-" * 50)
    print("2a. Plotting Group 2 sn-2 chains...")
    plot_group_sn2_comparison(group2_data, 'Order_Group2_PG_65_35')
    print("2b. Plotting Group 2 sn-1 chains...")
    plot_group_sn1_comparison(group2_data, 'Order_Group2_PG_65_35')
    print("2c. Plotting Group 2 PC comparison...")
    plot_group_pc_comparison(group2_data, 'Order_Group2_PG_65_35')
    
    # Group 3 plots
    print("\n" + "-" * 50)
    print("GROUP 3: PC/PG Systems (95:5)")
    print("-" * 50)
    print("3a. Plotting Group 3 sn-2 chains...")
    plot_group_sn2_comparison(group3_data, 'Order_Group3_PG_95_5')
    print("3b. Plotting Group 3 sn-1 chains...")
    plot_group_sn1_comparison(group3_data, 'Order_Group3_PG_95_5')
    print("3c. Plotting Group 3 PC comparison...")
    plot_group_pc_comparison(group3_data, 'Order_Group3_PG_95_5')
    
    # All systems combined
    print("\n" + "-" * 50)
    print("ALL 27 SYSTEMS COMBINED")
    print("-" * 50)
    print("4a. Plotting chain type comparison...")
    plot_all_systems_chain_types()
    print("4b. Plotting all systems grid...")
    plot_all_systems_comparison()
    print("4c. Plotting chain type summary...")
    plot_chain_type_summary()
    
    print("\n" + "=" * 70)
    print("All plots generated successfully!")
    print("=" * 70)
    print("\nOutput files:")
    print("\n  GROUP 1 (PC/PA 80:20):")
    print("    - Order_Group1_PA_sn2.png/pdf")
    print("    - Order_Group1_PA_sn1.png/pdf")
    print("    - Order_Group1_PA_PC_comparison.png/pdf")
    print("\n  GROUP 2 (PC/PG 65:35):")
    print("    - Order_Group2_PG_65_35_sn2.png/pdf")
    print("    - Order_Group2_PG_65_35_sn1.png/pdf")
    print("    - Order_Group2_PG_65_35_PC_comparison.png/pdf")
    print("\n  GROUP 3 (PC/PG 90:10):")
    print("    - Order_Group3_PG_95_5_sn2.png/pdf")
    print("    - Order_Group3_PG_95_5_sn1.png/pdf")
    print("    - Order_Group3_PG_95_5_PC_comparison.png/pdf")
    print("\n  ALL SYSTEMS:")
    print("    - Order_All_Systems_Chain_Types.png/pdf")
    print("    - Order_All_27_Systems_Grid.png/pdf")
    print("    - Order_Chain_Type_Summary.png/pdf")
    print("=" * 70)
