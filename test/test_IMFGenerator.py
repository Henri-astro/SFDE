# general libraries
import pytest
import numpy as np

# own libraries
from src.IMFGenerator import cIMFGenerator
from src.massfunction import cMassFunction
from src.BaseSCCalc import ComputeDens


def test_ComputeAlpha():
    """tests the compute alpha function"""
    
    IMFGen = cIMFGenerator( 0.0, 20.0 )
    
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
    
    
    IMFGen = cIMFGenerator( -2.0, 20.0 )
    
    alpha = IMFGen.ComputeAlpha( 1e5 )
    
    assert 3 == len( alpha )
    assert 1.3 - 0.885654 == pytest.approx( alpha[0] )
    assert 2.3 - 0.885654 == pytest.approx( alpha[1] )
    assert 2.178333 == pytest.approx( alpha[2] )
    
    
def test_CheckUpperEnd():
    """tests the CheckUpperEnd function"""
    
    MF = cMassFunction( 100.0, [0.08,150.0], [1.3] )
    IMFGen = cIMFGenerator( 0.0, 20.0 )
    
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
    ZH = -2.0
    Mini = 1e5
    
    IMFGen = cIMFGenerator( ZH, 20.0 )
    
    MF = IMFGen.ComputeMF( Mini )
    
    assert 0.08 == pytest.approx( MF.Getbounds()[0] )
    assert 0.5 == pytest.approx( MF.Getbounds()[1] )
    assert 1.0 == pytest.approx( MF.Getbounds()[2] )
    assert 147.6862129 == pytest.approx( MF.Getbounds()[3] )
    
    assert 0.414346 == pytest.approx( MF.Getalphas()[0] )
    assert 1.414346 == pytest.approx( MF.Getalphas()[1] )
    assert 2.178333 == pytest.approx( MF.Getalphas()[2] )
    
    assert Mini == pytest.approx( MF.GetMtot() )


def test_ComputeIMFFromToday():
    """checks that the initial mass function is computed correctly"""
    
    #create the initial values of a realistic GC:
    Mini = 1e6
    ZH = -2.0
    
    #orbital parameters
    Rapo = 8.0
    e = 0.5
    
    Rperi = Rapo * ( 1.0 - e ) / ( 1.0 + e )
    
    #the age of the GC in Myr
    t = 12.0
    
    #compute the present-day mass for that GC
    IMFGen = cIMFGenerator( ZH, 20.0 )
    M = IMFGen.ComputeCurrentMass( Mini, Rapo, Rperi, t )
    
    #do the calculations
    IMF = IMFGen.ComputeIMFFromToday( M, t, Rapo, Rperi )
    
    assert Mini == pytest.approx( IMF.GetMtot() )
    
    #2nd case
    Mini2 = 1e7
    ZH2 = -1.5
    
    #compute the present-day mass for that GC
    IMFGen2 = cIMFGenerator( ZH2, 20.0 )
    M2 = IMFGen2.ComputeCurrentMass( Mini2, Rapo, Rperi, t )
    
    #do the calculations
    IMF2 = IMFGen2.ComputeIMFFromToday( M2, t, Rapo, Rperi )
    
    assert Mini2 == pytest.approx( IMF2.GetMtot() )
    
    ## test NGC 2298
    Mini3 = 4.98e4
    ZH3 = -1.62
    Rapo3 = 16.49
    Rperi3 = 0.97
    
    #compute the present-day mass for that GC
    IMFGen2 = cIMFGenerator( ZH3, 30.0 )
    M3 = IMFGen2.ComputeCurrentMass( Mini3, Rapo3, Rperi3, t )
    
    #do the calculations
    IMF3 = IMFGen2.ComputeIMFFromToday( M3, t, Rapo3, Rperi3 )
    
    assert Mini3 == pytest.approx( IMF3.GetMtot() )
