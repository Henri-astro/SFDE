import numpy as np

G = 0.00449     # gravitational constant in pc^3 / (Myr^2 Msun)


def ComputeCrossingTime( Mini, r ):
    """Computes the crossing time of a SC with mass Mini
    Mini: the initial mass of the cluster [Msun]
    r: the cluster half-mass radius
    returns the crossing time [Myr]"""

    return pow( 2.0 * r, 1.5 ) / np.sqrt( G * Mini )


def ComputeRelaxTime( Mini, N, r ):
    """Computes the relaxtion time of a SC with mass Mini
    Mini: the initial mass of the cluster [Msun]
    N: the number of stars in the cluster
    r: the cluster half-mass radius
    returns the relaxation time [Myr]"""
    
    return 0.1 * N / np.log( N ) * ComputeCrossingTime( Mini, r )


def MassToRadius( mass ):
    """uses the Marks et al. (2012) mass-radius relationship to compute the half mass radius from the initial cluster mass
    mass: the mass of the cluster [Msun]
    returns the initial half-mass radius [pc]"""
    
    return 0.10 * pow( mass, 0.13 )


def ComputeDens( Mini ):
    """computes the dencity from the initial SC mass according to Marks et al. (2012)
    Mini: the initial cluster mass [Msun]
    returns the initial dencity [Msun/pc^3]"""
    
    return pow( 10.0, 0.61 * np.log10( Mini ) + 2.08 )


def ComputeMini( Dens ):
    """computes the initial mass from the initial SC densitz
    Dens: the initial dencity [Msun/pc^3]
    returns the initial cluster mass [Msun]"""
    
    return pow( Dens, 1.0/0.61 ) * pow( 10.0, -2.08/0.61 )
