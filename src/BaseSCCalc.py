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


def ComputeTidalRadius( mass, Rperi ):
    """computes the tidal radius of a GC
    mass: cluster mass [Msun]
    Rperi: the distance of the pericentre from the galactic centre [kpc]
    returns the tidal radius [pc]"""
    
    return pow( 4.44e-8 * mass, 1.0/3.0 ) * pow( Rperi * 1000.0, 2.0/3.0 )


def ComputeKingRadius( r, sfe ):
    """a function to compute the King radius
    r: the half-mass radius of the cluster [pc]
    returns the King radius [pc]"""
    
    return np.sqrt( 3.0/sfe ) * r


def ComputeConcentrationParametre( mass, Rperi, sfe ):
    """computes the concentration of an embedded SC
    mass: the mass of the SC [Msun]
    Rperi: the pericentre distance from the SC to the galactic centre [kpc]
    sfe: the star formation efficiency
    returns the concentration parametre"""
    
    r_t = ComputeTidalRadius( mass, Rperi )
    r_0 = ComputeKingRadius( MassToRadius( mass ), sfe )
    
    return np.log10( r_t/r_0 )


def ComputeKingConcentrationParameter( c ):
    """computes the King concentration parameter from the concentration parameter
    c: the concentration parameter of the SC
    returns the King concentration parameter W0"""
    
    return 4.38 * c
