import pytest
import numpy as np

from src.massfunction import cMassFunction
    

def test_GetterInit():
    """checks that all the getters of the MF are initialized correctly"""
    
    Mini = 1.0
    bounds = [0.08, 1.0]
    alphas = [1.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert MF.GetMtot() == Mini
    assert MF.Getbounds() == bounds
    assert MF.Getalphas() == alphas
    assert MF.GetAlpha( 0 ) == alphas[0]
    
    Mini = 3.0
    bounds = [0.08, 1.0, 2.0]
    alphas = [1.0, 2.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert MF.GetMtot() == Mini
    assert MF.Getbounds() == bounds
    assert MF.Getalphas() == alphas
    
    for nElem in range( len( alphas )):
        assert MF.GetAlpha( nElem ) == alphas[nElem]

    
def test_GetterValueErrors():
    """tests that the getters through value errors when expected"""
    
    MF = cMassFunction( 3.0, [0.08, 1.0, 2.0], [1.0, 2.0] )
    
    with pytest.raises( ValueError, match = "massfunction: GetFunctionValue: mass out of bounds!" ):
        MF.GetFunctionValue( 0.05 )
        
    with pytest.raises( ValueError, match = "massfunction: GetFunctionValue: mass out of bounds!" ):
        MF.GetFunctionValue( 3.0 )
        
    with pytest.raises( ValueError, match = "massfunction: GetMassDensity: mass out of bounds!" ):
        MF.GetMassDensity( 0.05 )
        
    with pytest.raises( ValueError, match = "massfunction: GetMassDensity: mass out of bounds!" ):
        MF.GetMassDensity( 3.0 )
        
    with pytest.raises( ValueError, match = "massfunction: GetAlpha: index out of bounds!" ):
        MF.GetAlpha( -1 )
        
    with pytest.raises( ValueError, match = "massfunction: GetAlpha: index out of bounds!" ):
        MF.GetAlpha( 3 )
        
    with pytest.raises( ValueError, match = "massfunction: GetK: index out of bounds!" ):
        MF.GetK( -1 )
        
    with pytest.raises( ValueError, match = "massfunction: GetK: index out of bounds!" ):
        MF.GetK( 3 )


def test_ComputeKs():
    """checks that the k-values for a mass function are computed correctly"""
    
    MF = cMassFunction( 1.0, [0.08, 1.0], [1.0] )
    
    assert len( MF.Getks() ) == 1
    assert MF.Getks()[0] == 0.0
    
    MF = cMassFunction( 2.0, [0.08, 1.0], [1.0] )
    
    assert len( MF.Getks() ) == 1
    assert pytest.approx( MF.Getks()[0] ) == 1.0 / 0.92
    assert pytest.approx( MF.GetK( 0 ) ) == 1.0 / 0.92
    
    MF = cMassFunction( 3.0, [0.08, 1.0, 2.0], [1.0, 2.0] )
    
    assert len( MF.Getks() ) == 2
    assert pytest.approx( MF.Getks()[1] ) == MF.Getks()[0]
    assert pytest.approx( MF.Getks()[0] ) == 1.0 / ( 0.92 + np.log( 2.0 ) )
    assert pytest.approx( MF.GetK( 0 ) ) == MF.Getks()[0]
    assert pytest.approx( MF.GetK( 1 ) ) == MF.Getks()[1]
    
    MF = cMassFunction( 2.0, [0.08, 0.5, 1.0], [1.3, 2.3] )
    
    assert pytest.approx( MF.Getks()[0] ) == 1.0 / 1.02081201235762066279
    assert pytest.approx( MF.Getks()[1] ) == 0.5 / 1.02081201235762066279
    assert pytest.approx( MF.GetK( 0 ) ) == MF.Getks()[0]
    assert pytest.approx( MF.GetK( 1 ) ) == MF.Getks()[1]
    

def test_GetFunctionValue():
    """tests if the function value is returned correctly"""
    
    MF = cMassFunction( 2.0, [0.08, 1.0], [1.0] )
    
    assert pytest.approx( MF.GetFunctionValue( 0.08 ) ) == 1.0 / ( 0.92 * 0.08 )
    assert pytest.approx( MF.GetFunctionValue( 0.5 ) ) == 1.0 / ( 0.92 * 0.5 )
    assert pytest.approx( MF.GetFunctionValue( 1.0 ) ) == 1.0 / 0.92
    
    MF = cMassFunction( 2.0, [0.08, 0.5, 1.0], [1.3, 2.3] )
    
    assert pytest.approx( MF.GetFunctionValue( 0.08 ) ) == pow( 0.08, -1.3 ) * MF.Getks()[0]
    assert pytest.approx( MF.GetFunctionValue( 0.5 ) ) == pow( 0.5, -1.3 ) * MF.Getks()[0]
    assert pytest.approx( MF.GetFunctionValue( 0.5 ) ) == pow( 0.5, -2.3 ) * MF.Getks()[1]
    assert pytest.approx( MF.GetFunctionValue( 1.0 ) ) == 1.0 * MF.Getks()[1]


def test_ComputeIntegral():
    """tests if the integrals are computed correctly, also contains some tests for GetNumbers"""
    
    Mini = 1.0
    bounds = [0.08, 1.0]
    alphas = [1.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert pytest.approx( MF.ComputeIntegral( 0.08, 1.0 ) ) == 0.0
    assert pytest.approx( MF.ComputeIntegral( 0.08, 2.0 ) ) == 0.0
    assert pytest.approx( MF.GetNumbers( 0.08, 1.0 ) ) == 1.0
    assert pytest.approx( MF.GetNumbers( 0.08, 2.0 ) ) == 1.0
    assert pytest.approx( MF.GetTotNumbers() ) == 1.0
    
    Mini = 2.0
    bounds = [0.08, 1.0]
    alphas = [1.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert pytest.approx( MF.ComputeIntegral( 0.08, 1.0, 1.0 ) ) == 1.0
    assert pytest.approx( MF.ComputeIntegral( 0.08, 2.0 ) ) == 1.0 / 0.92 * ( np.log( 1.0 ) - np.log( 0.08 ))
    
    Mini = 3.0
    bounds = [0.08, 1.0, 2.0]
    alphas = [1.0, 2.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert pytest.approx( MF.ComputeIntegral( 0.08, 2.0, 1.0 ) ) == 1.0
    assert pytest.approx( MF.ComputeIntegral( 0.08, 1.0, 1.0 ) ) == 1.0 / ( 0.92 + np.log( 2.0 ) ) * 0.92
    
    
def test_GetMass():
    """tests if the mass and mass fraction in a given mass interval is computed correctly"""
    
    Mini = 1.0
    bounds = [0.08, 1.0]
    alphas = [1.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert pytest.approx( MF.GetMass( 0.08, 1.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.08, 2.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.5, 2.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.5, 1.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.5, 0.6 ) ) == 0.0
    
    assert pytest.approx( MF.GetMassPortion( 0.08, 1.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.08, 2.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.5, 2.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.5, 1.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.5, 0.6 ) ) == 0.0
    
    Mini = 3.0
    bounds = [0.08, 1.0, 2.0]
    alphas = [1.0, 2.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    assert pytest.approx( MF.GetMass( 0.08, 2.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.08, 3.0 ) ) == Mini
    assert pytest.approx( MF.GetMass( 0.08, 1.0 ) ) == 1.0 / ( 0.92 + np.log( 2.0 ) ) * 0.92
    
    assert pytest.approx( MF.GetMassPortion( 0.08, 2.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.08, 3.0 ) ) == 1.0
    assert pytest.approx( MF.GetMassPortion( 0.08, 1.0 ) ) == 1.0 / ( 0.92 + np.log( 2.0 ) ) * 0.92 / Mini


def test_GetMassStarMinX():
    """checks that the mass of a star X below the given one is computed correctly"""
    
    Mini = 2.0
    bounds = [0.08, 1.0]
    alphas = [1.0]
    
    MF = cMassFunction( Mini, bounds, alphas )
    
    Masses = [ 0.1 * i for i in range( 1,10 ) ]
    
    for mass1 in Masses:
        for mass2 in Masses:
            if not mass2 > mass1:
                break
    
            Num = MF.GetNumbers( mass1, mass2 )
            
            assert mass1 == pytest.approx( MF.GetMassStarMinX( mass2, Num ) )
    
    alphas2 = [2.0]
    
    MF2 = cMassFunction( Mini, bounds, alphas )
    
    for mass1 in Masses:
        for mass2 in Masses:
            if not mass2 > mass1:
                break
    
            Num = MF2.GetNumbers( mass1, mass2 )
            
            assert mass1 == pytest.approx( MF2.GetMassStarMinX( mass2, Num ) )
            
    Mini3 = 1000.0
    bounds3 = [0.08,0.5,100.0]
    alphas3 = [1.3,2.3]
    
    MF3 = cMassFunction( Mini3, bounds3, alphas3 )
    
    Masses3 = [ 1.1 * i for i in range(1,90) ]
    
    for mass1 in Masses3:
        for mass2 in Masses3:
            if not mass2 > mass1:
                break
    
            Num = MF3.GetNumbers( mass1, mass2 )
            
            assert mass1 == pytest.approx( MF3.GetMassStarMinX( mass2, Num ) )
    
    #that case that didn't want to work
    Mini4 = 1000.0
    bounds4 = [0.3,0.5,0.7]
    alphas4 = [1.3,2.3]
    
    MF4 = cMassFunction( Mini4, bounds4, alphas4 )
    
    assert 0.301058871413287 == pytest.approx( MF4.GetMassStarMinX( 0.7, 2188.800573307484 ))
