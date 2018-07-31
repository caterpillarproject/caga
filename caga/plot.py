from __future__ import absolute_import, division

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

__all__ = ["mass_history"]

def mass_history(gamtree):
    """
    Plot halo mass of a gamma tree, colored by branch
    """
    br_halo_ID, br_age, br_z, br_t_merge, br_ID_merge, \
    br_m_halo, br_r_vir, br_is_prim, redshifts, times, \
        tree_trunk_ID = [x for x in gamtree.all_vars]
    nb_redshifts = len(redshifts)
    
    # Set figure frame and font properties
    fig = plt.figure(figsize=(6.5,4.5))
    matplotlib.rcParams.update({'font.size': 13})

    # Plot the evolution of the virial mass of each branch
    # Non star-forming branches will be shown as thin gray lines
    for i_z in range(nb_redshifts):
        for i_b in range(len(br_halo_ID[i_z])):
            plt.plot(br_z[i_z][i_b], br_m_halo[i_z][i_b], '.')

    # Set labels
    plt.xlabel('Redshift', fontsize=14)
    plt.ylabel('M$_\mathrm{vir}$ [M$_\odot$]', fontsize=14)
    #plt.legend(fontsize=13, frameon=False)
    plt.yscale('log')

    # Adjust the figure
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(right=0.95)
    plt.subplots_adjust(left=0.12)
    plt.subplots_adjust(bottom=0.15)
    return fig
