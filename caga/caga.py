from __future__ import absolute_import, division

"""
Modules for using Caterpillar-GAMMA data
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from NuPyCEE import omega, sygma
from JINAPyCEE import gamma, omega_plus

## hardcode some index labels for GAMMA tree data
all_ixs = np.arange(11)
ix_halo_ID, ix_age, ix_z, ix_t_merge, ix_ID_merge, ix_m_halo, ix_r_vir, \
    ix_is_prim, ix_redshifts, ix_times, ix_tree_trunk_ID = all_ixs

def load_gamma_tree(path):
    """ Load a gamma tree (from a .npy save) """
    assert path.endswith(".npy")
    gtree = np.load(path, encoding='latin1')
    assert len(all_ixs) == len(gtree)
    return gtree

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

def precompute_ssps():
    """
    Pre-calculate the ejecta of simple stellar populations.
    You can ignore the warning.
    """
    o_for_SSPs = omega.omega(special_timesteps=2, pre_calculate_SSPs=True)
    # Copy the SSPs array
    SSPs_in = [o_for_SSPs.ej_SSP, o_for_SSPs.ej_SSP_coef, o_for_SSPs.dt_ssp, o_for_SSPs.t_ssp]
    return SSPs_in
    
