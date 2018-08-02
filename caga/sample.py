from __future__ import absolute_import, division

"""
Tools to take samples from chemical evolution tracks
"""

import numpy as np
from . import caga, calc

default_elems = ["C","N","O","Na","Mg","Si","Ca","Ti","Mn","Fe","Co","Ni","Zn","Sr","Ba","Eu"]

def sample_omega(om, Nsample, t0=0.0, elems=default_elems, convert_to_XH=True):
    """
    Take the chemical evolution track from an omega run and sample stars along it.
    Output is Nsample x N_values_to_sample array
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

def sample_pmf(pmf, N):
    """ Get N indices from the given PMF """
    return np.random.choice(np.arange(len(pmf)), size=N, p=pmf)

