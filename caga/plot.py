from __future__ import absolute_import, division

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from . import caga, calc

__all__ = ["mass_history","metallicity_evolution"]

def mass_history(gamtree):
    """
    Plot halo mass of a gamma tree, colored by branch
    """
    br_halo_ID, br_age, br_z, br_t_merge, br_ID_merge, \
    br_m_halo, br_r_vir, br_is_prim, redshifts, times, \
    tree_trunk_ID, br_is_sub = [x for x in gamtree.all_vars]
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

def metallicity_evolution(g, tree_trunk_ID=None):
    """
    Plot metallicity evolution of each branch
    Input: g = gamma calculation
    Optional: tree_trunk_ID if you want to mark a particular branch
    """
    # Set figure frame and font properties
    fig = plt.figure(figsize=(6.5,4.5))
    matplotlib.rcParams.update({'font.size': 13})
    
    # For each star-forming branch ..
    for i_z in range(len(g.redshifts)):
        for i_b in range(len(g.br_halo_ID[i_z])):
            if g.br_is_SF[i_z][i_b]:
                
                # Copy the OMEGA (star-forming) calculation
                o = g.galaxy_inst[i_z][i_b].inner
                
                # Check if there is star formation.  Even if this is
                # technically a star-forming branch, it is possible
                # that other criteria such as gas density threshold
                # does not allow star formation.
                if sum(o.history.m_locked) > 0.0:
                
                    # Copy the times and metallicities
                    t_temp = np.array(o.history.age) + g.times[i_z]
                    Z_temp = o.history.metallicity
                
                    # If the galaxy eventually merges ..
                    if not g.br_halo_ID[i_z][i_b][-1] == tree_trunk_ID:
                    
                        # Copy the timestep at which the galaxy (branch) merged
                        # Please see the comment below for why the OMEGA+ calculation
                        # was not simply stopped when the galaxy merged..
                        i_t_merger = o.i_t_merger
                        
                        # Plot the evolution of metallicity from the formation
                        # of the branch until it merges
                        plt.plot(t_temp[:i_t_merger+1], Z_temp[:i_t_merger+1])
    
                    # Plot the complete evolution if this is the root/trunk
                    else:
                        plt.plot(t_temp, Z_temp, linewidth=5, label='Root of the merger tree')
                    
    # Set labels
    plt.yscale('log')
    plt.xlabel('Time [yr]', fontsize=14)
    plt.ylabel('Metallicity (mass fraction)', fontsize=14)
    plt.legend(fontsize=13, frameon=False)
    plt.ylim(5e-6,2e-2)
    
    # Adjust the figure
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(right=0.95)
    plt.subplots_adjust(left=0.12)
    plt.subplots_adjust(bottom=0.15)
    
    return fig

def sfh_evolution(g, tree_trunk_ID=None):
    """
    Plot sfh evolution of each branch
    Input: g = gamma calculation
    Optional: tree_trunk_ID if you want to mark a particular branch
    """
    fig = plt.figure(figsize=(6.5,4.5))
    matplotlib.rcParams.update({'font.size': 13})
    
    # For each star-forming branch ..
    for i_z in range(len(g.redshifts)):
        for i_b in range(len(g.br_halo_ID[i_z])):
            if g.br_is_SF[i_z][i_b]:
                
                # Copy the OMEGA (star-forming) calculation
                o = g.galaxy_inst[i_z][i_b].inner
                
                # Check if there is star formation
                if sum(o.history.m_locked) > 0.0:
                
                    # Copy the times and star formation rates
                    t_temp = np.array(o.history.age) + g.times[i_z]
                    sfh_temp = o.history.sfr_abs
                
                    # If the galaxy eventually merges ..
                    if not g.br_halo_ID[i_z][i_b][-1] == tree_trunk_ID:
                    
                        # Plot the evolution of metallicity for from the formation
                        # of the branch until it merges
                        i_t_merger = o.i_t_merger
                        plt.plot(t_temp[:i_t_merger+1], sfh_temp[:i_t_merger+1])
    
                    # Plot the complete evolution if this is the root/trunk
                    else:
                        plt.plot(t_temp, sfh_temp, linewidth=5, label='Root of the merger tree')
                    
    # Set labels
    plt.yscale('log')
    plt.xlabel('Time [yr]', fontsize=14)
    plt.ylabel('Star formation rate [M$_\odot$ yr$^{-1}$]', fontsize=14)
    plt.legend(fontsize=13, frameon=False)
    plt.ylim(1e-4,20)
    
    # Adjust the figure
    plt.subplots_adjust(top=0.95)
    plt.subplots_adjust(right=0.95)
    plt.subplots_adjust(left=0.12)
    plt.subplots_adjust(bottom=0.15)
    
    return fig

def mstar_evolution(g, reduction_factor=1.0):
    """
    Plot total mstar(z)
    """
    # saving history of total stellar mass
    mstar_hist = calc.mstar_evolution(g)
    #
    ## For each star-forming branch ..
    #mstar_tot = 0
    #for i_z in range(len(g.redshifts)):
    #    for i_b in range(len(g.br_halo_ID[i_z])):
    #        if g.br_is_SF[i_z][i_b]:
    #            
    #            # Copy the OMEGA (star-forming) calculation
    #            o = g.galaxy_inst[i_z][i_b].inner
    #            
    #            mstar_tot += sum(o.history.m_locked)
    #    
    #    mstar_hist[i_z] = mstar_tot
    
    # Set figure frame and font properties
    fig = plt.figure(figsize=(6,4))
    matplotlib.rcParams.update({'font.size': 13})
    
    plt.plot(g.redshifts,mstar_hist*reduction_factor)
    plt.yscale('log')
    plt.xlabel('redshift')
    plt.ylabel('Total $M_{star}$ [M$_\odot$]')
    return fig, (g.redshifts, mstar_hist)

def metallicity_distribution(g, sigma_gauss=0.1,
                             Fe_H_min=-6.0, Fe_H_max=1.0, d_Fe_H=0.05, 
                             plotxLow=-3, plotxHigh=0.5,logScale=False,
                             plotyLow=None, plotyHigh=None):
    fig = plt.figure(figsize=(6,4))
    mdf_x, mdf_all_norm = calc.mdf(g, Fe_H_min=Fe_H_min, Fe_H_max=Fe_H_max, d_Fe_H=d_Fe_H)
    mdf_all_gauss_norm = caga.convolve_gauss(mdf_x, mdf_all_norm, sigma_gauss)
    plt.plot(mdf_x,mdf_all_norm, label="Raw MDF")
    plt.plot(mdf_x,mdf_all_gauss_norm, label="Smoothed MDF")
    plt.xlim(plotxLow,plotxHigh)
    plt.ylim(plotyLow,plotyHigh)
    if logScale:
        plt.gca().set_yscale('log')
    matplotlib.rcParams.update({'font.size': 14.0})
    plt.xlabel('[Fe/H]', fontsize=14)
    plt.ylabel('MDF', fontsize=14)
    
    return fig, (mdf_x, mdf_all_norm, mdf_all_gauss_norm)
