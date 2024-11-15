import numpy as np
import pytest

from src.BaseSCCalc import *


def test_ComputeCrossingTime():
    """tests that the crossing time is computed correctly"""
    
    epsilon = 1e-5
    
    Mini = 1.0 / G
    
    Res = ComputeCrossingTime( Mini, 3.0 )
    
    assert pytest.approx( Res, epsilon ) == pow( 2.0 * 3.0, 1.5 )


def test_ComputeRelaxTime():
    """tests that the relaxation time is computed correctly"""
    
    epsilon = 1e-5
    
    Mini = 1.0 / G
    N = 10
    
    Res = ComputeRelaxTime( Mini, N, 3.0 )
    
    assert pytest.approx( Res, epsilon ) == pow( 2.0 * 3.0, 1.5 ) / np.log( 10.0 )


def test_MassToRadius():
    """tests the mass to radius function"""
    
    epsilon = 1e-5
    
    mass = 1.0
    
    Res = MassToRadius( mass )
    
    assert pytest.approx( Res, epsilon ) == 0.1
    
    
def test_ComputeDens():
    """tests the transformation from mass to density and vice versa"""
    
    epsilon = 1e-5
    
    Mini = 10.0
    
    Res = ComputeDens( Mini )
    
    assert pytest.approx( Res, epsilon ) == pow( 10.0, 0.61 + 2.08 )
    
    for Mini in [1.0,10.0,20.0,100.0,150.0]:
        Res = ComputeMini( ComputeDens( Mini ))
        
        assert pytest.approx( Res, epsilon ) == Mini


def test_ComputeTidalRadius():
    """tests that the tidal radius is computed correctly"""
    
    #values for NGC 104 according to Baumgardts catalogue
    M = 8.53e5
    R = 7.52
    r_t = 124.63
    
    assert pytest.approx( ComputeTidalRadius( M, R ), 5 ) == r_t
    
    #values for NGC 288 according to Baumgardts catalogue
    M = 9.62e4
    R = 12.21
    r_t = 95.16
    
    assert pytest.approx( ComputeTidalRadius( M, R ), 2e-1 ) == r_t


def test_ComputeKingRadius():
    """tests that the King radius is computed correctly"""
    
    rh = 1.0
    sfe = 0.3
    r0 = 3.162278
    
    assert pytest.approx( ComputeKingRadius( rh, sfe ), 1e-4 ) == r0


def test_ComputeConcentrationParametre():
    """tests that the concentration parameter is computed correctly"""
    
    #values for NGC 288 according to Baumgardts catalogue (apart from sfe)
    M = 9.62e4
    Rperi = 12.21
    sfe = 0.3
    
    assert pytest.approx( ComputeConcentrationParametre( M, Rperi, sfe ), 1e-1 ) == 1.8


def test_ComputeKingConcentrationParameter():
    """tests that the King concentration parameter is computed correctly"""
    
    c = 1.0
    W0 = 4.38
    
    assert pytest.approx( ComputeKingConcentrationParameter( c ) ) == W0
