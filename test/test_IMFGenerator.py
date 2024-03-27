# general libraries
import pytest
import numpy as np

# own libraries
from src.IMFGenerator import cIMFGenerator
from src.massfunction import cMassFunction
from src.BaseSCCalc import ComputeDens
from src.RemnantCalculator import cRemnantCalculator


def SetupRemCalc( ZH ):
    """sets up a remnand calculator with the desired metallicity
    ZH: the metallicity to use"""
    
    #initialize the remnant calculator
    data = { "mass[Msun]": [0.0,1.0,2.0], "t_-1.0": [9.0,6.0,3.0], "t_0.0": [9.0,6.5,3.2], "Mfin_-1.0": [0.5,0.8,1.2], "Mfin_0.0":[0.4,0.9,1.5] }
    RemCalc = cRemnantCalculator( data, ZH )
    
    return RemCalc


def test_ComputeAlpha():
    """tests the compute alpha function"""
    
    IMFGen = cIMFGenerator( SetupRemCalc( 0.0 ) )
    
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
    
    
    IMFGen = cIMFGenerator( SetupRemCalc( -2.0 ) )
    
    alpha = IMFGen.ComputeAlpha( 1e5 )
    
    assert 3 == len( alpha )
    assert 1.3 - 0.885654 == pytest.approx( alpha[0] )
    assert 2.3 - 0.885654 == pytest.approx( alpha[1] )
    assert 2.178333 == pytest.approx( alpha[2] )
    
    
def test_CheckUpperEnd():
    """tests the CheckUpperEnd function"""
    
    MF = cMassFunction( 100.0, [0.08,150.0], [1.3] )
    IMFGen = cIMFGenerator( SetupRemCalc( 0.0 ) )
    
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
    
    IMFGen = cIMFGenerator( SetupRemCalc( ZH ) )
    
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
    
    #initialize the remnant calculator
    data = { "mass[Msun]": [0.0,1.0,2.0], "t_-1.0": [9.0,6.0,3.0], "t_0.0": [9.0,6.5,3.2], "Mfin_-1.0": [0.5,0.8,1.2], "Mfin_0.0":[0.4,0.9,1.5] }
    RemCalc = cRemnantCalculator( data, ZH )
    
    #compute the present-day mass for that GC
    IMFGen = cIMFGenerator( RemCalc )
    M = IMFGen.ComputeCurrentMass( Mini, Rapo, Rperi, t )
    
    assert M < Mini
    
    #do the calculations
    IMF = IMFGen.ComputeIMFFromToday( M, t, Rapo, Rperi )
    
    assert Mini == pytest.approx( IMF.GetMtot() )
    
    #2nd case
    Mini2 = 1e7
    ZH2 = -1.5
    
    RemCalc2 = cRemnantCalculator( data, ZH2 )
    
    #compute the present-day mass for that GC
    IMFGen2 = cIMFGenerator( RemCalc2 )
    M2 = IMFGen2.ComputeCurrentMass( Mini2, Rapo, Rperi, t )
    
    #do the calculations
    IMF2 = IMFGen2.ComputeIMFFromToday( M2, t, Rapo, Rperi )
    
    assert Mini2 == pytest.approx( IMF2.GetMtot() )
    
    # test NGC 2298 (not anymore)
    Mini3 = 1e6
    ZH3 = -1.62
    Rapo3 = 16.49
    Rperi3 = 0.97
    
    data2 = { "mass[Msun]": [-1.0969100130080565,1.2811318113306203,1.4905601716686168,1.7004783792473719,1.9095126511112555,2.109930732634646,2.176091259055681],\
        "t_-1.67": [14.781965892799771,7.0116178624137095,6.782656630008358,6.6271309462829295,6.5458504858751345,6.536877443455831,6.548735818688512],\
        "t_-0.67": [14.500836618524938,7.036950102574134,6.808820219192039,6.649878535398668,6.560969046845507,6.540720164040784,6.548004881970323],\
        "t_-0.37": [14.278246979261908,7.0436651762878135,6.816705090939301,6.655845470724989,6.56195154148037,6.534042677503096,6.538179560916891],\
        "t_-0.2": [14.082641523009046,7.031120566944764,6.810206123382868,6.65376944320497,6.5626514317370495,6.535878176086212,6.540061982173659],\
        "Mfin_-1.67": [-1.0683387603155188,0.30531051804439324,0.38214857243620165,0.45199662615017105,0.514416249102295,0.5675803560346638,0.5836940180895561],\
        "Mfin_-0.67": [-1.199145508496439,0.4917341867301212,0.6157796654378671,0.7360768395473158,0.8518482214315842,0.9590824335965844,0.9936720522957804],\
        "Mfin_-0.37": [-1.2547690154718592,0.6007366361903759,0.749345119848031,0.8958958127226768,1.0394371222929097,1.1748198121819418,1.2190293852418737],\
        "Mfin_-0.2": [-1.2964507017617695,0.6508813453057455,0.8132025672636466,0.9744135854483941,1.1334635389646466,1.2845685581693704,1.3341516505881326] }
    
    RemCalc3 = cRemnantCalculator( data2, ZH3 )
    
    #compute the present-day mass for that GC
    IMFGen3 = cIMFGenerator( RemCalc3 )
    M3 = IMFGen3.ComputeCurrentMass( Mini3, Rapo3, Rperi3, t )
    
    #do the calculations
    IMF3 = IMFGen3.ComputeIMFFromToday( M3, t, Rapo3, Rperi3 )
    
    assert Mini3 == pytest.approx( IMF3.GetMtot() )
