from __future__ import absolute_import, division

"""
Modules for using Caterpillar-GAMMA data
"""

import numpy as np

from NuPyCEE import omega, sygma
from JINAPyCEE import gamma, omega_plus

from . import plot

class gamma_tree(object):
    """
    A class that just stores the output of GAMMA merger tree data
    """
    def __init__(self, br_halo_ID, br_age, br_z, br_t_merge, br_ID_merge,
                 br_m_halo, br_r_vir, br_is_prim, redshifts, times,
                 tree_trunk_ID, ):
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
    
        self.m_DM_0 = list(filter(None, self.br_m_halo))[-1][-1][-1]
        
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
            self.redshifts, self.times, self.tree_trunk_ID

    @property
    def kwargs(self):
        """
        Create kwargs to pass to GAMMA with merger tree data.
        Does not include 
        """
        kwargs = {"redshifts":self.redshifts, "times":self.times, "br_halo_ID":self.br_halo_ID,
                  "br_age":self.br_age, "br_z":self.br_z, "br_t_merge":self.br_t_merge, 
                  "br_ID_merge":self.br_ID_merge, "br_m_halo":self.br_m_halo,
                  "br_r_vir":self.br_r_vir, "br_is_prim":br_is_prim,
                  "tree_trunk_ID":self.tree_trunk_ID,
                  "m_DM_0":self.m_DM_0}
        #"br_is_SF":br_is_SF, 
        #"pre_calculate_SSPs":True,"SSPs_in":SSPs_in, "print_off":True}
        return kwargs

    def plot_mass_history(self):
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
    

def get_root_Mstar(gam):
    """
    """
def get_root_ZH(gam):
    """
    """
def get_root_FeH(gam):
    """
    """
def get_root_FeH_distr(gam):
    """
    """
def get_root_Mpeak(gam):
    """
    """

