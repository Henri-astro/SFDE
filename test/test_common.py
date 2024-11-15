#general modules
import pytest

#own modules
from src.common import *


def test_Isnumber():
    """tests that numbers are recognised correctly"""
    
    assert True == IsNumber( "0" )
    assert True == IsNumber( "5.0" )
    assert True == IsNumber( "-6.7843" )
    assert False == IsNumber( "bla" )
    assert False == IsNumber( "73s" )
    assert False == IsNumber( "1k4i2je" )
    assert False == IsNumber( "" )


def test_LinInterExtrapolate():
    """tests that the inter- and extrapolation is done correctly"""
    
    assert 1.0 == pytest.approx( LinInterExtrapolate( [0.5,1.0], [0.8,9.2], 0.5 ))
    assert 9.2 == pytest.approx( LinInterExtrapolate( [0.5,1.0], [0.8,9.2], 0.8 ))
    assert 2.0 == pytest.approx( LinInterExtrapolate( [0.5,1.0], [0.8,1.6], 1.0 ))
    assert -0.2 == pytest.approx( LinInterExtrapolate( [0.5,-1.0], [0.8,1.4], 0.6 ))
    assert 3.0 == pytest.approx( LinInterExtrapolate( [0.5,-1.0], [0.8,1.4], 1.0 ))
