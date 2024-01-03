# general modules
import pytest

#own modules
from src.RemnantCalculator import cRemnantCalculator


def Setup():
    """instantiates a standart version of the class for testing"""
    
    data = { "Mstar": [1.0,2.0], "t_-1.0": [6.0,3.0], "t_0.0": [6.5,3.2], "Mfin_-1.0": [0.8,1.2], "Mfin_0.0":[0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc


def test_FindValues():
    """tests the find values function from the remnant calculator"""
    
    data = { "Mstar": [1.0,2.0], "t_-1.0": [6.0,3.0], "t_0.0": [6.5,3.2], "Mfin_-1.0": [0.8,1.2], "Mfin_0.0":[0.9,1.5] }
    ZH = -0.5
    
    RemCalc = cRemnantCalculator( data, ZH )
    
    Res = RemCalc.FindValues( data, ZH, "t_" )
    
    assert 6.25 == pytest.approx( Res[0] )
    assert 3.1 == pytest.approx( Res[1] )
    
    Res2 = RemCalc.FindValues( data, -1.0, "Mfin_" )
    
    assert 0.8 == pytest.approx( Res2[0] )
    assert 1.2 == pytest.approx( Res2[1] )
    
    
def test_GetValFromList():
    """tests if values are returned """
    
    RemCalc = Setup()
    
    ListA = [ -1.0,1.0,3.0,6.0 ]
    ListB = [ 2.0,2.5,3.1,4.0 ]
    
    assert 1.75 == pytest.approx( RemCalc.GetValFromList( -2.0, ListA, ListB ))
    assert 2.0 == pytest.approx( RemCalc.GetValFromList( -1.0, ListA, ListB ))
    assert 2.25 == pytest.approx( RemCalc.GetValFromList( 0.0, ListA, ListB ))
    assert 2.95 == pytest.approx( RemCalc.GetValFromList( 2.5, ListA, ListB ))
    assert 3.1 == pytest.approx( RemCalc.GetValFromList( 3.0, ListA, ListB ))
    assert 4.3 == pytest.approx( RemCalc.GetValFromList( 7.0, ListA, ListB ))


def test_GetTimeFromMass():
    """tests the function to return the time corresponding to a given mass and vice versa"""
    
    RemCalc = Setup()
    
    assert pow( 10.0, 3.4 ) == pytest.approx( RemCalc.GetTimeFromMass( 1.0 ) )
    assert pow( 10.0, 0.25 ) == pytest.approx( RemCalc.GetTimeFromMass( 10.0 ) )
    
    assert 10.0 == pytest.approx( RemCalc.GetMassFromTime( pow( 10.0, 0.25 ) ) )
    assert 8.5 == pytest.approx( RemCalc.GetMassFromTime( pow( 10.0, 0.47233038399997789002 ) ) )
    
    masses = [8.5, 20.0, 35.5, 40.7, 79.8, 121.0, 134.5 ]
    
    for mass in masses:
        assert mass == pytest.approx( RemCalc.GetMassFromTime( RemCalc.GetTimeFromMass( mass ) ))


def test_GetMfinFromMass():
    """tests that the remnant masses are determined correctly"""
    
    RemCalc = Setup()
    
    assert pow( 10.0, 0.35 ) == pytest.approx( RemCalc.GetMfinFromMass( 1.0 ) )
    assert pow( 10.0, 0.85 ) == pytest.approx( RemCalc.GetMfinFromMass( 10.0 ) )
    assert pow( 10.0, 1.1 ) == pytest.approx( RemCalc.GetMfinFromMass( pow( 10.0, 1.5 ) ) )
    assert pow( 10.0, 1.35 ) == pytest.approx( RemCalc.GetMfinFromMass( 100.0 ) )
    assert pow( 10.0, 1.85 ) == pytest.approx( RemCalc.GetMfinFromMass( 1000.0 ) )
