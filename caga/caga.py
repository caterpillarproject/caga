from __future__ import absolute_import, division

"""
Modules for using Caterpillar-GAMMA data
"""

import numpy as np
import itertools
from scipy import interpolate, optimize

from NuPyCEE import omega, sygma, stellab
from JINAPyCEE import gamma, omega_plus

from . import plot, calc


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

def mass_fraction_to_spectro(MH,el):
    """
    Given mass ratio M/H for element el, compute [X/H]
    """
    assert (el in el_to_mass) and (el in el_to_solar), el
    NH = MH/el_to_mass[el]
    solar = el_to_solar[el]
    return np.log10(NH) - solar + 12

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

def generate_kwargs(gt, m_vir_thresh, SSPs_in=None):
    kwargs = gt.kwargs
    br_is_SF = br_is_SF_thresh(gt, m_vir_thresh)
    if SSPs_in is None:
        SSPs_in = precompute_ssps()
    kwargs.update({"br_is_SF":br_is_SF,
                   "pre_calculate_SSPs":True,
                   "SSPs_in":SSPs_in})
    return kwargs

def plot_lumfn(Mstars,MoverL=2.0):
    """
    Not implemented properly
    """
    logLbins = np.arange(1,12,.2)
    logLs = np.log10(Mstars/MoverL)
    h,x = np.histogram(logLs, logLbins)
    raise NotImplementedError
    
### Utility functions for MDF
def normalize_distribution(x,y):
    def _findBin(x,xRange):
        """ Find which bin x is in"""
        ind = 0
        n = len(xRange)
        if x > xRange[0]:
            while ind < n and xRange[ind] < x :
                ind += 1
        return ind
    def _areaUnderCurve(y,x,a,b):
        """ Use to normalize distributions """
        aIndex = _findBin(a,x)
        bIndex = _findBin(b,x)
        return np.trapz(y[aIndex:bIndex+1],x[aIndex:bIndex+1])
    x_min = np.min(x)
    d_x = np.max(np.diff(x))
    x_max = np.max(x)+d_x
    y_norm = y/_areaUnderCurve(y,x,x_min,x_max)
    return y_norm
    
def convolve_gauss(x, y, sigma_gauss):
    """
    Simple convolution to smooth a distribution
    This can probably be replaced with scipy or numpy
    """
    # Convolve with Gaussian
    y_gauss = np.zeros(len(x))
    for i_x in range(len(x)):
        for i_x_scan in range(len(x)):
            y_gauss[i_x_scan] += y[i_x]*np.exp((-1*(x[i_x_scan] - x[i_x])**2) / (2.0*sigma_gauss**2))
    y_gauss_norm = normalize_distribution(x,y_gauss)
    return y_gauss_norm

def get_cdf(x,y, return_arrays=False, **kwargs):
    """
    Compute cumulative distribution function and return interpolating function
    Assumes x is left bin edge!!!
    """
    assert len(x) == len(y)
    assert np.all(np.diff(x) > 0), np.diff(x)
    ycum = np.array([0]+list(y.cumsum()))
    ycum = ycum/ycum[-1]
    dx = np.max(np.diff(x))
    xcum = np.array(list(x)+[x[-1]+dx])
    cdf = interpolate.interp1d(xcum, ycum, assume_sorted=True, **kwargs)
    if return_arrays:
        return cdf, xcum, ycum
    return cdf

def find_distribution_percentile(x,y,p):
    """
    Given a PDF x,y find the CDF and linearly interpolate to find where the percentile is p.
    Requires x to be sorted in increasing order with no duplicates.
    """
    p = np.ravel(p)
    assert np.all(np.logical_and(p >= 0, p <= 1))
    cdf, xcum, ycum = get_cdf(x,y,return_arrays=True)
    p_guess = [xcum[np.argmin(np.abs(ycum-p_))] for p_ in p]
    output = [optimize.brentq(lambda x: cdf(x)-p_, xcum[0], xcum[-1], xtol=1e-4) for p_ in p]
    return np.array(output)

def find_distribution_median(x,y):
    return find_distribution_percentile(x,y,0.5)[0]

def get_xmid(x):
    assert np.all(np.diff(x) > 0), np.diff(x)
    dx = np.max(np.diff(x))
    xall = np.array(list(x)+[x[-1]+dx])
    return (xall[1:]+xall[:-1])/2.

def find_distribution_mean(x,y):
    assert len(x) == len(y)
    xmid = get_xmid(x)
    ynorm = y.sum()
    return np.sum(xmid*y/ynorm)

def find_distribution_std(x,y):
    E_X = find_distribution_mean(x,y)
    xmid = get_xmid(x)
    ynorm = y.sum()
    E_X2 = np.sum(xmid*xmid*y/ynorm)
    return np.sqrt(E_X2 - E_X**2)
