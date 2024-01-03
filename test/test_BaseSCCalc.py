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
