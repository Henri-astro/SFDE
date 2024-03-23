#general libs
import pytest

#own libs
from src.common import *
from src.DataProcessor import cDataProcessor
from src.IMFGenerator import cIMFGenerator
from src.RemnantCalculator import cRemnantCalculator
from src.DataReader import cDataReader


class cMockData():
    """a mock class to hold all the required data"""
    
    def __init__( self ):
        self.__GCData = {}
        
    
    def AccessGCData( self, name ):
        """get an element of the GC data"""
        return self.__GCData[name]
    
    
    def AddGCData( self, name, values ):
        """add an element to the GC data"""
        self.__GCData[name] = values
    
    
    def SetRemnantData( self, data ):
        """sets the remnant data"""
        
        self.__RemnantData = data
    
    
    def GetRemnantData( self ):
        """returns a copy of the complete remnant data"""
        
        return self.__RemnantData.copy()
    
    
    def SNExplodes( self, mass ):
        return True
    
    
    def Ejecta( self, nextMass ):
        return 0.074


def test_ComputeIMF():
    """tests that the IMF is generated correctly"""
    
    ###############################
    # setup some mock data        #
    ###############################
    
    data = cMockData()
    
    compareData = { "Mini": (5e5, 1e6),\
                    "Age": (12.0, 11.0),\
                    "R_a": (8.0, 8.0),\
                    "R_p": (6.0, 10.0),\
                    "Fe-H": (-2.0, -0.5),\
                    "FeSpread": (0.05, 0.1),\
                    "SFE": (0.3, 0.3)}
    
    RemnantReader = cDataReader( "test/mockdata/RemnantData.dat", ["Mstar"] )
    data.SetRemnantData( RemnantReader.GetData() )
    
    compareData["Mass"] = []
    compareData["alphas"] = []
    compareData["ProducedIron"] = []
    compareData["NSN"] = []
    compareData["mlast"] = []
    compareData["SFD"] = []
    
    for nGC in range( len( compareData["Mini"] ) ):
        RemCalc = cRemnantCalculator( RemnantReader.GetData(), compareData["Fe-H"][nGC] + 0.3 )
        IMFGen = cIMFGenerator( RemCalc )
        MF = IMFGen.ComputeMF( compareData["Mini"][nGC] )
        compareData["Mass"].append( IMFGen.ComputeCurrentMass( compareData["Mini"][nGC], compareData["R_a"][nGC], compareData["R_p"][nGC], compareData["Age"][nGC] ) )
        compareData["alphas"].append( IMFGen.ComputeAlpha( compareData["Mini"][nGC] ))
        
        compareData["ProducedIron"].append( pIronSun * ( pow( 10.0, compareData["Fe-H"][nGC] + compareData["FeSpread"][nGC] ) - pow( 10.0, compareData["Fe-H"][nGC] - compareData["FeSpread"][nGC] )) * compareData["Mini"][nGC] * ( 1.0 - compareData["SFE"][nGC] ) / compareData["SFE"][nGC] )
        
        compareData["NSN"].append( int( compareData["ProducedIron"][nGC] / data.Ejecta( None ) + 0.999999 ) )
        compareData["mlast"].append( MF.GetMassStarMinX( MF.Getbounds()[-1], compareData["NSN"][nGC] - 1 ) )
        compareData["SFD"].append( RemCalc.GetTimeFromMass( compareData["mlast"][nGC] ) )
    
    compareData["Mass"] = tuple( compareData["Mass"] )
    
    for Elem in compareData:
        if Elem in ["alphas", "ProducedIron", "NSN", "mlast"]:
            continue
        
        data.AddGCData( Elem, compareData[Elem] )
    
    ###############################
    # do the calculations         #
    ###############################
    
    Proc = cDataProcessor()
    
    Proc.ProcessData( data )
    
    ###############################
    # check the results           #
    ###############################
    
    for nCompElem in range( len( compareData["Mini"] ) ):
        Minis = data.AccessGCData( "Mini" )
        IMFs = data.AccessGCData( "IMF" )
        
        assert compareData["Mini"][nCompElem] == pytest.approx( Minis[nCompElem] )
        assert compareData["Mini"][nCompElem] == pytest.approx( IMFs[nCompElem].GetMtot() )
        assert 0.08 == pytest.approx( IMFs[nCompElem].Getbounds()[0] )
        assert 0.5 == pytest.approx( IMFs[nCompElem].Getbounds()[1] )
        assert 1.0 == pytest.approx( IMFs[nCompElem].Getbounds()[2] )
        
        for nAlpha in range( len( compareData["alphas"][nCompElem] )):
            assert compareData["alphas"][nCompElem][nAlpha] == pytest.approx( IMFs[nCompElem].Getalphas()[nAlpha] )
        
        assert compareData["ProducedIron"][nCompElem] == pytest.approx( data.AccessGCData( "ProducedIron" )[nCompElem] )
        assert compareData["NSN"][nCompElem] == pytest.approx( data.AccessGCData( "NSN" )[nCompElem] )
        assert compareData["mlast"][nCompElem] == pytest.approx( data.AccessGCData( "mlast" )[nCompElem] )
        assert compareData["SFD"][nCompElem] == pytest.approx( data.AccessGCData( "SFD" )[nCompElem] )
