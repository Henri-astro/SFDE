# general libraries
import pytest
import numpy as np

# own libraries
from src.IMFGenerator import cIMFGenerator
from src.massfunction import cMassFunction
from src.BaseSCCalc import ComputeDens


def test_ComputeAlpha():
    """tests the compute alpha function"""
    
    IMFGen = cIMFGenerator( 0.0 )
    
    alpha = IMFGen.ComputeAlpha( 1e4 )
    
    assert 3 == len( alpha )
    assert 1.3 == pytest.approx( alpha[0] )
    assert 2.3 == pytest.approx( alpha[1] )
    assert 2.3 == pytest.approx( alpha[2] )
    
    alpha = IMFGen.ComputeAlpha( 1e5 )
    
    assert 3 == len( alpha )
    assert 1.3 == pytest.approx( alpha[0] )
    assert 2.3 == pytest.approx( alpha[1] )
    assert 2.293133 == pytest.approx( alpha[2] )
    
    
    IMFGen = cIMFGenerator( -2.0 )
    
    alpha = IMFGen.ComputeAlpha( 1e5 )
    
    assert 3 == len( alpha )
    assert 1.3 - 0.885654 == pytest.approx( alpha[0] )
    assert 2.3 - 0.885654 == pytest.approx( alpha[1] )
    assert 2.178333 == pytest.approx( alpha[2] )
    
    
def test_CheckUpperEnd():
    """tests the CheckUpperEnd function"""
    
    MF = cMassFunction( 100.0, [0.08,150.0], [1.3] )
    IMFGen = cIMFGenerator( 0.0 )
    
    assert -1.0 == pytest.approx( IMFGen.CheckUpperEnd( MF ))
    
    MF = cMassFunction( 2.0, [0.08, 1.0], [1.0] )
    
    assert np.log( 150.0 ) / 0.92 - 1.0 == pytest.approx( IMFGen.CheckUpperEnd( MF ))
    
    MF = cMassFunction( 3.0, [0.08, 1.0, 2.0], [1.0, 2.0] )
    
    assert ( pow( 2.0, -1.0 ) - pow( 150.0, -1.0 ) ) / ( 0.92 + np.log( 2.0 ) ) - 1.0 == pytest.approx( IMFGen.CheckUpperEnd( MF ))
    
    MF = cMassFunction( 2.0, [0.08, 0.5, 1.0], [1.3, 2.3] )
    
    assert ( pow( 1.0, -1.3 ) - pow( 150.0, -1.3 ) ) * 0.5 / ( 1.02081201235762066279 * 1.3 ) - 1.0 == pytest.approx( IMFGen.CheckUpperEnd( MF ))


def test_ComputeMF():
    """checks that the mass function is computed correctly"""
    
    #create test case
    FeH = -2.0
    Mini = 1e5
    
    IMFGen = cIMFGenerator( FeH )
    
    MF = IMFGen.ComputeMF( Mini )
    
    assert 0.08 == pytest.approx( MF.Getbounds()[0] )
    assert 0.5 == pytest.approx( MF.Getbounds()[1] )
    assert 1.0 == pytest.approx( MF.Getbounds()[2] )
    assert 147.6862129 == pytest.approx( MF.Getbounds()[3] )
    
    assert 0.414346 == pytest.approx( MF.Getalphas()[0] )
    assert 1.414346 == pytest.approx( MF.Getalphas()[1] )
    assert 2.178333 == pytest.approx( MF.Getalphas()[2] )
    
    assert Mini == pytest.approx( MF.GetMtot() )
