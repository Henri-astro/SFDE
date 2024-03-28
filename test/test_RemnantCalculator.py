# general modules
import pytest
import numpy as np

#own modules
from src.RemnantCalculator import cRemnantCalculator
from src.massfunction import cMassFunction


def Setup():
    """instantiates a standart version of the class for testing"""
    
    data = { "mass[Msun]": [1.0,2.0], "t_-1.0": [6.0,3.0], "t_0.0": [6.5,3.2], "Mfin_-1.0": [0.8,1.2], "Mfin_0.0":[0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc
   

def Setup2():
    """instantiates a standart version of the class for testing"""
    
    data = { "mass[Msun]": [0.1,0.4,0.5,0.8,1.0,2.0], "t_-1.0": [6.0,3.0,2.7,2.5,2.2,2.1], "t_0.0": [6.5,3.2,3.0,2.9,2.8,2.7], "Mfin_-1.0": [0.08,0.3,0.4,0.6,0.8,1.2], "Mfin_0.0":[0.08,0.25,0.4,0.7,0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc
    

def test_GetTimeFromMass():
    """tests the function to return the time corresponding to a given mass and vice versa"""
    
    RemCalc = Setup()
    
    assert pow( 10.0, 0.4 ) == pytest.approx( RemCalc.GetTimeFromMass( 1.0 ) )
    assert pow( 10.0, -2.75 ) == pytest.approx( RemCalc.GetTimeFromMass( 10.0 ) )
    
    assert 10.0 == pytest.approx( RemCalc.GetMassFromTime( pow( 10.0, -2.75 ) ) )
    assert 8.5 == pytest.approx( RemCalc.GetMassFromTime( pow( 10.0, 0.47233038399997789002 - 3.0 ) ) )
    
    masses = [0.1, 0.5, 1.0, 2.0, 4.0, 8.5, 20.0, 35.5, 40.7, 79.8, 121.0, 134.5 ]
    
    for mass in masses:
        assert mass == pytest.approx( RemCalc.GetMassFromTime( RemCalc.GetTimeFromMass( mass ) ))
        
    RemCalc2 = Setup2()
    
    for mass in masses:
        assert mass == pytest.approx( RemCalc2.GetMassFromTime( RemCalc2.GetTimeFromMass( mass ) ))


def test_GetMfinFromMass():
    """tests that the remnant masses are determined correctly"""
    
    RemCalc = Setup()
    
    assert pow( 10.0, 0.35 ) == pytest.approx( RemCalc.GetMfinFromMass( 1.0 ) )
    assert pow( 10.0, 0.85 ) == pytest.approx( RemCalc.GetMfinFromMass( 10.0 ) )
    assert pow( 10.0, 1.1 ) == pytest.approx( RemCalc.GetMfinFromMass( pow( 10.0, 1.5 ) ) )
    assert pow( 10.0, 1.35 ) == pytest.approx( RemCalc.GetMfinFromMass( 100.0 ) )
    assert pow( 10.0, 1.85 ) == pytest.approx( RemCalc.GetMfinFromMass( 1000.0 ) )


def test_GetMfinFromMassFunct():
    """tests that the left-over mass is computed correctly"""
    
    data = { "mass[Msun]": [0.5,1.0,2.0], "t_-1.0": [8.0,6.0,3.0], "t_0.0": [8.0,6.5,3.2], "Mfin_-1.0": [0.5,0.8,1.2], "Mfin_0.0":[0.4,0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    IMF = cMassFunction( 2.0, [0.08,0.5,1.0], [1.3,2.3] )
    
    assert IMF.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF, 1.0 ))
    assert IMF.GetMass( IMF.Getbounds()[0], RemCalc.GetMassFromTime( 10.0 ) ) + IMF.GetMass( RemCalc.GetMassFromTime( 10.0 ), IMF.Getbounds()[-1] ) / IMF.GetMtot() * RemCalc.GetMfinFromLogMass( 0.5 * ( np.log10( RemCalc.GetMassFromTime( 10.0 )) + np.log10( IMF.Getbounds()[-1] ))) == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF, 10.0 ))
    
    IMF2 = cMassFunction( 100.0, [0.08,0.5,1.0,15.0], [1.3,2.3,2.2] )
    
    assert IMF2.GetMass( IMF2.Getbounds()[0], RemCalc.GetMassFromTime( 1.0 ) ) + ( IMF2.GetMass( RemCalc.GetMassFromTime( 1.0 ), pow( 10.0, 0.75 ) ) * RemCalc.GetMfinFromLogMass( 0.5 * ( np.log10( RemCalc.GetMassFromTime( 1.0 )) + 0.75 )) + IMF2.GetMass( pow( 10.0, 0.75 ), 15.0 ) * RemCalc.GetMfinFromLogMass( 0.5 * ( np.log10( 15.0 ) + 0.75  )) ) / IMF2.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF2, 1.0 ))
    assert IMF2.GetMass( IMF2.Getbounds()[0], RemCalc.GetMassFromTime( 1.1e-3 ) ) + ( IMF2.GetMass( RemCalc.GetMassFromTime( 1.1e-3 ), 15.0 ) * RemCalc.GetMfinFromLogMass( 0.5 * ( np.log10( RemCalc.GetMassFromTime( 1.1e-3 )) + np.log10( 15.0 ) )) ) / IMF2.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF2, 1.1e-3 ))
