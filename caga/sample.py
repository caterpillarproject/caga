from __future__ import absolute_import, division

"""
Tools to take stellar population samples from chemical evolution tracks
"""

import numpy as np
import itertools
from . import caga, calc

default_elems = ["C","N","O","Na","Mg","Si","Ca","Ti","Mn","Fe","Co","Ni","Zn","Sr","Ba","Eu"]

def sample_gamma(g, Nsample, elems=default_elems, convert_to_XH=True):
    """
    Take a GAMMA tree and sample stars from chemical evolution track from an omega run and sample stars along it.
    Require a minimum of 1 star per branch, so will not get back an array of length Nsample.
    """
    ## Create output array.
    ## Add 2 extra indices for i_z and i_b
    Ncols = 1+len(elems) + 2
    
    ## Create mass in each branch and find total mass
    mstar_tot = 0.
    br_mstar = [[] for i_z in range(len(g.redshifts))]
    for i_z in range(len(g.redshifts)):
        for i_b in range(len(g.br_halo_ID[i_z])):
            if g.br_is_SF[i_z][i_b]: # only sample from star forming branches
                o = g.galaxy_inst[i_z][i_b].inner
                mstar = np.sum(o.history.m_locked)
                mstar_tot += mstar
                br_mstar[i_z].append(mstar)
            else:
                br_mstar[i_z].append(0.)
    
    ## Find the number of stars to sample on each branch
    br_N_stars = [[max(int(round(Nsample*mstar/mstar_tot,0)),1) for mstar in z_mstar] for z_mstar in br_mstar]
    total_N_stars = np.sum(list(itertools.chain(*br_N_stars)))
    ## TODO: What to do if Nsample is not high enough?
    output = np.zeros((total_N_stars, Ncols))
    
    all_samples = []
    all_iz = []
    all_ib = []
    for i_z in range(len(g.redshifts)):
        for i_b in range(len(g.br_halo_ID[i_z])):
            N_stars = br_N_stars[i_z][i_b]
            o = g.galaxy_inst[i_z][i_b].inner
            samples = sample_omega(o, N_stars, t0=g.times[i_z], elems=elems, convert_to_XH=convert_to_XH)
            all_samples.append(samples)
            all_iz += [i_z]*N_stars
            all_ib += [i_b]*N_stars
    np.concatenate(all_samples,out=output[:,0:-2])
    output[:,-2] = all_iz
    output[:,-1] = all_ib
    return output
    
def sample_omega(om, Nsample, t0=0.0, elems=default_elems, convert_to_XH=True):
    """
    Take the chemical evolution track from an omega run and sample stars along it.
    Default output columns:
      t_form, element mass ratios (M_X/M_H) for C, N, O, Na, Mg, Si, Ca, Ti, Mn, Fe, Co, Ni, Zn, Sr, Ba, Eu
    
    Parameters
    ----------
    om : `NuPyCEE.omega.omega`
        An omega model to sample

    Nsample : int
        Number of stars to sample from the track
    
    t0 : float, optional
        Value to add to the age of the stars (default 0.0)
    
    elems : list, optional
        Elements to pull out.
    
    convert_to_XH : bool, optional
        if True (default), converts mass ratio arrays to [X/H]
    
    Returns
    -------
    numpy array of size (Nsample, 1+len(elems))
    """
    
    ## Get element indices to grab from the array
    ## This will throw an error if you give it bad elements
    elem_indices = np.array([om.history.elements.index(elem) for elem in ["H"]+elems])
    ## Pull out ism_elem_yield and divide by hydrogen
    ism_elem_yield = np.array(om.history.ism_elem_yield)[:,elem_indices]
    ism_elem_yield = ((ism_elem_yield[:,1:]).T/ism_elem_yield[:,0]).T
    
    ## Create output
    ## columns: t_form, elem1, elem2, ...
    output = np.zeros((Nsample, len(elem_indices)))
    
    ## Get random samples in the time bins, weighted by m_locked
    pmf = om.history.m_locked/np.sum(om.history.m_locked)
    indices = sample_pmf(pmf, Nsample)
    
    ## Pull out the time of formation
    output[:,0] = np.array(om.history.age)[indices] + t0
    ## Pull out elements
    if convert_to_XH:
        for j, el in enumerate(elems):
            output[:,1+j] = caga.mass_fraction_to_spectro(ism_elem_yield[indices,j], el)
    else:
        output[:,1:] = ism_elem_yield[indices,:]
        
    return output

def fill_nonfinite(x, val):
    """ Take an array x and fill ~np.isfinite(x) with val """
    x = x.copy()
    x[~np.isfinite(x)] = val
    return x

def drop_nonfinite(x):
    """ Take an array x and keep finite """
    return x[np.isfinite(x)]

def sample_pmf(pmf, N):
    """ Get N indices from the given PMF """
    return np.random.choice(np.arange(len(pmf)), size=N, p=pmf)

