# general modules
import pytest
import numpy as np

#own modules
from src.RemnantCalculator import cRemnantCalculator
from src.massfunction import cMassFunction
from src.DataReader import cDataReader


def Setup0():
    """instantiates a standart version of the class for testing"""
    
    data = { "mass[Msun]": [1.0,2.0], "t_-1.0": [6.0,3.0], "t_0.0": [6.5,3.2], "Mfin_-1.0": [0.8,1.2], "Mfin_0.0":[0.9,1.5] }
    ZH = -1.0
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc


def Setup1():
    """instantiates a standart version of the class for testing"""
    
    data = { "mass[Msun]": [1.0,8.0], "t_-1.0": [6.0,3.0], "t_0.0": [6.5,3.2], "Mfin_-1.0": [0.9,5.5], "Mfin_0.0":[0.8,4.0] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc
   

def Setup2():
    """instantiates a standart version of the class for testing"""
    
    data = { "mass[Msun]": [0.1,0.4,0.5,0.8,1.0,2.0], "t_-1.0": [6.0,3.0,2.7,2.5,2.2,2.1], "t_0.0": [6.5,3.2,3.0,2.9,2.8,2.7], "Mfin_-1.0": [0.08,0.3,0.4,0.6,0.8,1.2], "Mfin_0.0":[0.08,0.25,0.4,0.7,0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc


def SetupRealData( ZH ):
    """instantiates a class for testing using actual data from Yan et al."""
    
    dataReader = cDataReader( "test/mockdata/RemnantDataComplete.dat", ["mass[Msun]"] )
    data = dataReader.GetData()
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc
    

def test_GetTimeFromMass():
    """tests the function to return the time corresponding to a given mass and vice versa"""
    
    RemCalc0 = Setup0()
    
    assert pow( 10.0, -3.0 ) == pytest.approx( RemCalc0.GetTimeFromMass( 1.0 ) )
    
    RemCalc1 = Setup1()
    
    assert pow( 10.0, 6.1201265 - 9.0 ) == pytest.approx( RemCalc1.GetTimeFromMass( 1.0 ) )
    assert pow( 10.0, 2.1703146 - 9.0 ) == pytest.approx( RemCalc1.GetTimeFromMass( 10.0 ) )
    
    assert 10.0 == pytest.approx( RemCalc1.GetMassFromTime( pow( 10.0, 2.1703146 - 9.0 ) ) )
    assert 8.5 == pytest.approx( RemCalc1.GetMassFromTime( pow( 10.0, 2.8286166 - 9.0 ) ) )
    
    masses = [ 0.1, 0.5, 1.0, 2.0, 4.0, 8.5, 20.0, 35.5, 40.7, 79.8, 121.0, 134.5 ]
    
    for mass in masses:
        assert mass == pytest.approx( RemCalc1.GetMassFromTime( RemCalc1.GetTimeFromMass( mass ) ))
        
    RemCalc2 = Setup2()
    
    for mass in masses:
        assert mass == pytest.approx( RemCalc2.GetMassFromTime( RemCalc2.GetTimeFromMass( mass ) ))
        
        
def test_GetTimeFromMassRealData():
    """read in the actual data file to thst the TimeFromMass function"""
    
    RemCalc1 = SetupRealData( -1.67 )
    
    assert pow( 10.0, 14.781965892799773 - 9.0 ) == pytest.approx( RemCalc1.GetTimeFromMass( 0.07999999999999999 ) )
    assert pow( 10.0, 6.548735818688513 - 9.0 ) == pytest.approx( RemCalc1.GetTimeFromMass( 149.9999999999999 ) )
    
    RemCalc2 = SetupRealData( -1.25 )
    
    assert 0.0034701733382704793 == pytest.approx( RemCalc2.GetTimeFromMass( 135.32400682968088 ) )


def test_GetMfinFromMass():
    """tests that the remnant masses are determined correctly"""
    
    RemCalc = Setup1()
    
    assert 0.8759747 == pytest.approx( RemCalc.GetMfinFromMass( 1.0 ) )
    assert 6.3578049 == pytest.approx( RemCalc.GetMfinFromMass( 10.0 ) )
    assert 19.5280704 == pytest.approx( RemCalc.GetMfinFromMass( pow( 10.0, 1.5 ) ) )
    assert 61.1761067 == pytest.approx( RemCalc.GetMfinFromMass( 100.0 ) )
    assert 609.3591249 == pytest.approx( RemCalc.GetMfinFromMass( 1000.0 ) )


def test_GetMfinFromMassFunct():
    """tests that the left-over mass is computed correctly"""
    
    data = { "mass[Msun]": [0.5,1.0,2.0], "t_-1.0": [11.0,9.5,8.5], "t_0.0": [10.0,9.0,8.0], "Mfin_-1.0": [0.5,0.8,1.2], "Mfin_0.0":[0.4,0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    IMF = cMassFunction( 2.0, [0.08,0.5,1.0], [1.3,2.3] )
    
    assert IMF.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF, 1.0 ))
    assert IMF.GetMass( IMF.Getbounds()[0], RemCalc.GetMassFromTime( 10.0 ) ) + IMF.GetMass( RemCalc.GetMassFromTime( 10.0 ), IMF.Getbounds()[-1] ) / IMF.GetMtot() * RemCalc.GetMfinFromMass( 0.5 * ( RemCalc.GetMassFromTime( 10.0 ) +  IMF.Getbounds()[-1] )) == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF, 10.0 ))
    
    IMF2 = cMassFunction( 100.0, [0.08,0.5,1.0,15.0], [1.3,2.3,2.2] )
    
    print( RemCalc.GetMassFromTime( 1.0 ) )
    
    assert IMF2.GetMass( IMF2.Getbounds()[0], RemCalc.GetMassFromTime( 1.0 ) ) + IMF2.GetMass( RemCalc.GetMassFromTime( 1.0 ), 15.0 ) * RemCalc.GetMfinFromMass( 0.5 * ( 15.0 + RemCalc.GetMassFromTime( 1.0 )  )) / IMF2.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF2, 1.0 ))
    assert IMF2.GetMass( IMF2.Getbounds()[0], RemCalc.GetMassFromTime( 1.1e-3 ) ) + ( IMF2.GetMass( RemCalc.GetMassFromTime( 1.1e-3 ), 15.0 ) * RemCalc.GetMfinFromMass( 0.5 * ( RemCalc.GetMassFromTime( 1.1e-3 ) + 15.0 )) ) / IMF2.GetMtot() == pytest.approx( RemCalc.GetMfinFromMassFunct( IMF2, 1.1e-3 ))
