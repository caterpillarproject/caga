from __future__ import absolute_import, division

"""
Operations on gamma trees
"""

import numpy as np
from . import caga

def root_mstar(g):
    return mstar_evolution(g)[0]
def mstar_evolution(g):
    """
    Find mstar(z) (same order as g.redshifts)
    """
    mstar_hist = np.zeros(len(g.redshifts))
    mstar_tot = 0
    for i_z in range(len(g.redshifts)):
        for i_b in range(len(g.br_halo_ID[i_z])):
            if g.br_is_SF[i_z][i_b]:
                # Copy the OMEGA (star-forming) calculation
                o = g.galaxy_inst[i_z][i_b].inner
                mstar_tot += sum(o.history.m_locked)
        mstar_hist[i_z] = mstar_tot
    return mstar_hist

def root_FeH_mean(g):
    """ Convenience function to compute MDF and find mean """
    return caga.find_distribution_mean(*mdf(g))
def root_FeH_std(g):
    """ Convenience function to compute MDF and find stdev """
    return caga.find_distribution_std(*mdf(g))
def root_FeH_medscat(g,p=[.16,.5,.84]):
    """ Convenience function to compute MDF and find medscat """
    x,y=mdf(g)
    return caga.find_distribution_percentile(x,y,p)
def mdf(g, Fe_H_min=-6.0, Fe_H_max=2.0, d_Fe_H=0.05):
    """
    Loop through GAMMA tree to find the total metallicity distribution function
    """
    def get_Fe_H(inst):
        # Find the Fe index
        i_Fe = inst.history.elements.index('Fe')
        i_H = inst.history.elements.index('H')
        # Define the age and Fe_H array
        age = []
        Fe_H = []
        m_locked = []
        # For each timestep ..
        for i_t in range(inst.nb_timesteps):
            # If there are metals ..
            if inst.history.ism_elem_yield[i_t][i_Fe] > 0.0:
                # Calculate the metallicity 
                m_Fe_H_ratio = inst.history.ism_elem_yield[i_t][i_Fe] / inst.history.ism_elem_yield[i_t][i_H]
                Fe_H.append( np.log10(m_Fe_H_ratio) - np.log10((10**(7.50-12))*56.0) )
                # Copy the m_locked and age
                age.append(inst.history.age[i_t])
                m_locked.append(inst.history.m_locked[i_t])
        return age, Fe_H, m_locked    
    
    nb_Fe_H = int((Fe_H_max - Fe_H_min) / d_Fe_H) + 1
    mdf_x = np.zeros(nb_Fe_H)
    mdf_y = np.zeros(nb_Fe_H)
    for i_Fe in range(nb_Fe_H):
        mdf_x[i_Fe] = Fe_H_min + i_Fe * d_Fe_H
    # Create the MDF array of all galaxies
    mdf_inst = []
    for i_zz in range(len(g.redshifts)):
        mdf_inst.append([])
        for i_br in range(len(g.br_age[i_zz])):
            mdf_inst[i_zz].append(np.zeros(nb_Fe_H))
    # Return the [Fe/H] value
    # For all progenitor galaxies ..
    for i_zz in range(len(g.redshifts)):
        for i_br in range(len(g.galaxy_inst[i_zz])):
    
            # Verify if the galaxy formed stars ..
            go_for_it = False
            if len(g.br_is_SF) == 0:
                go_for_it = True
            elif g.br_is_SF[i_zz][i_br]:
                go_for_it = True
            if go_for_it:
                
                # Copy the galaxy instance
                inst = g.galaxy_inst[i_zz][i_br].inner
    
                # Get the Age - [Fe/H] relation and the stellar mass formed
                age, Fe_H, m_locked = get_Fe_H(inst)
    
                # For each timestep ..
                for i_t in range(len(age)-1):
    
                    if not np.isnan(Fe_H[i_t]) and not np.isnan(Fe_H[i_t+1]) and \
                       not np.isinf(Fe_H[i_t]) and not np.isinf(Fe_H[i_t+1]) and \
                        m_locked[i_t] > 0.0:
    
                        # Calculate the star formation rate for this timestep
                        sfr_inst = m_locked[i_t] / (age[i_t+1] - age[i_t])
                            
                        # If we do not need to interpolate ..
                        if Fe_H[i_t+1] == Fe_H[i_t]:
                            
                            # Find the lower MDF_x boundary
                            i_low = 0
                            while mdf_x[i_low+1] < Fe_H[i_t]:
                                i_low += 1
    
                            # Add the stellar mass in the MDF bin
                            mdf_inst[i_zz][i_br][i_low] += sfr_inst * (age[i_t+1] - age[i_t])
    
                            # Go to the next MDF bin
                            i_low += 1
    
                        # If we need to interpolate
                        else:
                            
                            # Calculate the interpolation coefficients
                            # t = a * [Fe/H] + b
                            aa = (age[i_t+1] - age[i_t]) / (Fe_H[i_t+1] - Fe_H[i_t])
                            bb = age[i_t] - aa * Fe_H[i_t]
    
                            # Calculate the maximum and minimum  Fe_H values
                            # So it doesn't matter when [Fe/H] is increasing or decreasing
                            min_Fe_inst = min(Fe_H[i_t], Fe_H[i_t+1])
                            max_Fe_inst = max(Fe_H[i_t], Fe_H[i_t+1])
    
                            # Find the lower MDF_x boundary
                            i_low = 0
                            while mdf_x[i_low+1] < min_Fe_inst:
                                i_low += 1
    
                            # While some MDF bins still need to be treated ..
                            while mdf_x[i_low] < max_Fe_inst:
    
                                # Get the lower and upper [Fe/H] boundaries covered
                                # by the OMEGA timesteps within the currend MDF bin
                                Fe_low = max(mdf_x[i_low], min_Fe_inst)
                                Fe_up = min(mdf_x[i_low+1], max_Fe_inst)
    
                                # Calculate the OMEGA time interval spent on the current MDF bin
                                t_low = aa * Fe_low + bb
                                t_up = aa * Fe_up + bb
                                dt_spent = abs(t_up - t_low)
    
                                # Add the stellar mass in the MDF bin
                                mdf_inst[i_zz][i_br][i_low] += sfr_inst * dt_spent
    
                                # Go to the next MDF bin
                                i_low += 1
    # Sum all MDFs
    mdf_all = np.zeros(nb_Fe_H)
    for i_zz in range(len(g.redshifts)):
        for i_br in range(len(g.galaxy_inst[i_zz])):
            mdf_all += mdf_inst[i_zz][i_br]
    mdf_all_norm = caga.normalize_distribution(mdf_x,mdf_all)
    
    return mdf_x, mdf_all_norm

def mdf_accreted(g_destroyed, Fe_H_min=-6.0, Fe_H_max=2.0, d_Fe_H=0.05, sigma_gauss=0.1):
    """ Creates the MDF of the stellar halo from all destroyed halos """
    # g_destroyed = array of gamma objects, one for each destroyed halo

    n_destroyed = len(g_destroyed)

    mdf_x, _ = mdf(g_destroyed[0], Fe_H_min, Fe_H_max, d_Fe_H)
    m = len(mdf_x)

    indivMDFs = np.zeros((n_destroyed,m))
    totalMDF = np.zeros(m)

    for i in range(0,n_destroyed):
        mstar_hist = mstar_evolution(g_destroyed[i])
        if mstar_hist[-1] != 0:
            _, mdf_all_norm = mdf(g_destroyed[i], Fe_H_min, Fe_H_max, d_Fe_H)
            mdf_smooth = caga.convolve_gauss(mdf_x, mdf_all_norm, sigma_gauss)
            indivMDFs[i] = mdf_smooth
            totalMDF += mdf_smooth*mstar_hist[-1]
        
    totalMDF = caga.normalize_distribution(mdf_x,totalMDF)

    return mdf_x, totalMDF
