# general libs
import pytest

# own libs
from src.massfunction import cMassFunction
from src.StarExtractor import cStarExtractor


def Setup():
    """sets up a StarExtractor with a given mass function"""
    
    Mini = 1000
    alphas = [1.3, 2.3]
    bounds = [0.08, 0.5, 100.0]
    
    massfunction = cMassFunction( Mini, bounds, alphas )
    
    return cStarExtractor( massfunction )


def test_GetNextMostMassiveStar():
    """checks that the most massive star is returned correctly"""

    StarExtractor = Setup()
    
    assert 100.0 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 37.25701942381608 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 24.516667055863117 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 18.72454481455272 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 15.342693919909495 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 13.099324107107401 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 11.490204391104232 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 10.273250422565809 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 9.316940190015194 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 8.543315163945714 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 7.903065746515376 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 7.363391521696144 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 6.901577171843005 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 6.501364477254261 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 6.150792506827555 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 5.840854504041673 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 5.564631384355116 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    assert 5.316714887089405 == pytest.approx( StarExtractor.GetNextMostMassiveStar() )
    
    # more automated test
    Mini = 10000
    alphas = [1.3, 2.3, 1.8]
    bounds = [0.08, 0.5, 1.0, 140.0]
    
    massfunction = cMassFunction( Mini, bounds, alphas )
    
    StarExtractor2 = cStarExtractor( massfunction )
    
    for nStar in range( 1000 ):
        nextMass = StarExtractor2.GetNextMostMassiveStar()
        
        assert nStar + 1 == pytest.approx( massfunction.GetNumbers( nextMass, bounds[-1] ) )
