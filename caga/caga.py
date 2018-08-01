from __future__ import absolute_import, division

"""
Modules for using Caterpillar-GAMMA data
"""

import numpy as np
import itertools

from NuPyCEE import omega, sygma, stellab
from JINAPyCEE import gamma, omega_plus

from . import plot


#### Define some element information
# Use stellab to get solar values
_stellab = stellab.stellab()
_asp09 = _stellab.sol_norm[_stellab.paths_norm.index("Asplund_et_al_2009")]
el_to_solar = dict([_elemnorm for _elemnorm in _asp09])
# Mass number of the most abundant isotopes
el_norm = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg',
           'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr',
           'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br',
           'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag',
           'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce',
           'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
           'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl',
           'Pb', 'Bi', 'Th', 'U']
mn_norm = [1.0000200000000001, 3.999834, 6.924099999999999, 9.0, 10.801,
           12.011061999999999, 14.00229, 16.004379000000004, 19.0, 20.13891,
           23.0, 24.3202, 27.0, 28.108575891424106, 31.0, 32.094199999999994,
           35.4844, 36.30859999999999, 39.13588999999999, 40.115629999999996,
           45.0, 47.9183, 50.9975, 52.05541, 55.0, 55.90993, 59.0, 58.759575,
           63.6166, 65.4682, 69.79784000000001, 72.69049999999999, 75.0, 79.0421,
           79.9862, 83.88084, 85.58312, 87.71136299999999, 89.0, 91.31840000000001,
           93.0, 96.05436999999998, 101.1598, 103.0, 106.5111, 107.96322,
           112.508, 114.9142, 118.8077, 121.85580000000002, 127.69839999999999,
           127.0, 131.28810000000001, 133.0, 137.42162, 138.99909000000002,
           140.20986, 141.0, 144.33351666483335, 150.4481, 152.0438, 157.3281,
           159.0, 162.57152999999997, 165.0, 167.32707000000002, 169.0,
           173.11618838116186, 175.02820514102572, 178.54163, 180.99988000000002,
           183.89069999999998, 186.28676, 190.28879711202887, 192.254,
           195.11347999999998, 197.0, 200.6297, 204.40952, 207.34285, 209.0,
           232.0, 237.27134]
el_to_mass = dict(zip(el_norm, mn_norm))

class gamma_tree(object):
    """
    A class that just stores the output of GAMMA merger tree data
    """
    def __init__(self, br_halo_ID, br_age, br_z, br_t_merge, br_ID_merge,
                 br_m_halo, br_r_vir, br_is_prim, redshifts, times,
                 tree_trunk_ID, br_is_sub):
        """
        Input and store gamma merger tree info.
        The constructor just takes all the variables needed.
        """
        self.br_halo_ID = br_halo_ID
        self.br_age = br_age
        self.br_z = br_z
        self.br_t_merge = br_t_merge
        self.br_ID_merge = br_ID_merge
        self.br_m_halo = br_m_halo
        self.br_r_vir = br_r_vir
        self.br_is_prim = br_is_prim
        self.redshifts = redshifts
        self.times = times
        self.tree_trunk_ID = tree_trunk_ID
        self.br_is_sub = br_is_sub
    
        mfilt = list(filter(None, self.br_m_halo))
        self.Mpeak = np.max(list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(mfilt)))))
        self.m_DM_0 = mfilt[-1][-1][-1]
        assert self.tree_trunk_ID == list(filter(None, self.br_halo_ID))[-1][-1][-1], "tree_trunk_ID does not match!"
        
    @classmethod
    def load(cls, path):
        """ Load a gamma tree (from a .npy save) """
        assert path.endswith(".npy")
        # These trees were generated on blender with python 2
        gtree = np.load(path, encoding='latin1')
        return cls(*gtree)
    
    @property
    def all_vars(self):
        """
        Shortcut for getting the full variable set.
        e.g. *gtree.all_vars
        """
        return self.br_halo_ID, self.br_age, self.br_z, self.br_t_merge, \
            self.br_ID_merge, self.br_m_halo, self.br_r_vir, self.br_is_prim, \
            self.redshifts, self.times, self.tree_trunk_ID, self.br_is_sub

    @property
    def kwargs(self):
        """
        Create kwargs to pass to GAMMA with merger tree data.
        """
        kwargs = {"redshifts":self.redshifts, "times":self.times, "br_halo_ID":self.br_halo_ID,
                  "br_age":self.br_age, "br_z":self.br_z, "br_t_merge":self.br_t_merge, 
                  "br_ID_merge":self.br_ID_merge, "br_m_halo":self.br_m_halo,
                  "br_r_vir":self.br_r_vir, "br_is_prim":self.br_is_prim,
                  "tree_trunk_ID":self.tree_trunk_ID,
                  #"br_is_sub":self.br_is_sub,
                  "m_DM_0":self.m_DM_0}
        return kwargs

    def plot_mass_history(self):
        """ Create a plot with the mass history """
        return plot.mass_history(self)
        
def precompute_ssps():
    """
    Pre-calculate the ejecta of simple stellar populations.
    You can ignore the warning.
    """
    o_for_SSPs = omega.omega(special_timesteps=2, pre_calculate_SSPs=True)
    # Copy the SSPs array
    SSPs_in = [o_for_SSPs.ej_SSP, o_for_SSPs.ej_SSP_coef, o_for_SSPs.dt_ssp, o_for_SSPs.t_ssp]
    return SSPs_in
    
def br_is_SF_thresh(gt, m_vir_thresh=0.):
    """
    Fill the input array that will inform GAMMA whether or not a given branch
    will form stars, with a simple mass threshold.
    By default, the threshold is 0.0, i.e. all halos form stars.
    """
    br_is_SF = []
    for i_z in range(len(gt.br_m_halo)):
        br_is_SF.append([])
        for i_b in range(len(gt.br_m_halo[i_z])):
            if max(gt.br_m_halo[i_z][i_b]) >= m_vir_thresh:
                br_is_SF[i_z].append(True)
            else:
                br_is_SF[i_z].append(False)
    return br_is_SF

def generate_kwargs(gt, m_vir_thresh):
    kwargs = gt.kwargs
    br_is_SF = br_is_SF_thresh(gt, m_vir_thresh)
    SSPs_in = precompute_ssps()
    kwargs.update({"br_is_SF":br_is_SF,
                   "pre_calculate_SSPs":True,
                   "SSPs_in":SSPs_in})
    return kwargs

def get_root_Mstar(g):
    """
    """
def get_root_ZH(g):
    """
    """
def get_root_FeH(g):
    """
    """
def get_root_FeH_distr(g):
    """
    """
def get_root_Mpeak(g):
    """
    """

def plot_lumfn(Mstars,MoverL=2.0):
    logLbins = np.arange(1,12,.2)
    logLs = np.log10(Mstars/MoverL)
    h,x = np.histogram(logLs, logLbins)
    
def mass_to_spectro(M1,M2,e1,e2):
    """
    Given two masses M1, M2 and elements e1, e2,
    return compute [e1/e2]
    """
    assert (e1 in el_to_mass) and (e2 in el_to_mass), (e1, e2)
    assert (e1 in el_to_solar) and (e2 in el_to_solar), (e1, e2)
    N1 = M1/el_to_mass[e1]
    N2 = M2/el_to_mass[e2]
    sol1 = el_to_solar[e1]
    sol2 = el_to_solar[e2]
    return np.log10(N1) - np.log10(N2) - sol1 + sol2

### Utility functions for MDF
def _findBin(x,xRange):
    ind = 0
    n = len(xRange)
    if x > xRange[0]:
        while ind < n and xRange[ind] < x :
            ind += 1
    return ind
def _areaUnderCurve(y,x,a,b):
    aIndex = _findBin(a,x)
    bIndex = _findBin(b,x)
    return np.trapz(y[aIndex:bIndex+1],x[aIndex:bIndex+1])

def compute_mdf(g, Fe_H_min=-6.0, Fe_H_max=1.0, d_Fe_H=0.05):
    def get_Fe_H(inst):
        # Find the Fe index
        i_Fe = inst.history.elements.index('Fe')
        # Define the age and Fe_H array
        age = []
        Fe_H = []
        m_locked = []
        # For each timestep ..
        for i_t in range(inst.nb_timesteps):
            # If there are metals ..
            if inst.ymgal[i_t][2] > 0.0:
                # Calculate the metallicity 
                m_Fe_H_ratio = inst.history.ism_elem_yield[i_t][i_Fe] / inst.history.ism_elem_yield[i_t][0]
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
    mdf_all_norm = mdf_all/_areaUnderCurve(mdf_all,mdf_x,Fe_H_min,Fe_H_max)
    
    return mdf_x, mdf_all_norm

def convolve_mdf(mdf_x, mdf_all, sigma_gauss):
    Fe_H_min = np.min(mdf_x)
    d_Fe_H = np.max(np.diff(mdf_x))
    Fe_H_max = np.max(mdf_x)+d_Fe_H
    # Convolve with Gaussian
    mdf_all_gauss = np.zeros(len(mdf_x))
    for i_x in range(len(mdf_x)):
        for i_x_scan in range(len(mdf_x)):
            mdf_all_gauss[i_x_scan] += mdf_all[i_x]*np.exp((-1*(mdf_x[i_x_scan] - mdf_x[i_x])**2) / (2.0*sigma_gauss**2))
    mdf_all_gauss_norm = mdf_all_gauss/_areaUnderCurve(mdf_all_gauss,mdf_x,Fe_H_min,Fe_H_max)
    return mdf_all_gauss_norm
