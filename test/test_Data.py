import pytest

from src.Data import cData


def SetupTest():
    """a function to setup the tests for the cData class"""
    
    return cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/Ejecta.dat", "test/mockdata/RemnantData.dat" )


def test_ConstructorExceptions():
    """tests that all exceptions in the constructor are thrown correctly"""
    
    with pytest.raises(ValueError, match=r"cData: Empty SN data. Please check your SN file."):
        cData( "test/mockdata/GCData.dat", "test/mockdata/0SNe.dat", "test/mockdata/Ejecta.dat", "test/mockdata/RemnantData.dat" )
        
    with pytest.raises(ValueError, match=r"cData: Empty Ejecta data. Please check your Ejecta file."):
        cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/0Ejecta.dat", "test/mockdata/RemnantData.dat" )
        
    with pytest.raises(ValueError, match=r"cData: Empty Remnant data. Please check your Remnant file."):
        cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/Ejecta.dat", "test/mockdata/0RemnantData.dat" )
        
    with pytest.raises(ValueError, match=r"cData: time information in Remnant data missing. Please check your Remnant file."):
        cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/Ejecta.dat", "test/mockdata/RemnantData_brokent.dat" )
        
    with pytest.raises(ValueError, match=r"cData: final masses in Remnant data missing. Please check your Remnant file."):
        cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/Ejecta.dat", "test/mockdata/RemnantData_brokenMfin.dat" )


def test_AccessGCData():
    """tests if the GCData is accessed apropriately"""
    
    data = SetupTest()
    
    assert 2 == len( data.AccessGCData( "Name" ) )
    assert "47_Tuc" == data.AccessGCData( "Name" )[0]
    assert "NGC_288" == data.AccessGCData( "Name" )[1]
    
    assert 2 == len( data.AccessGCData( "Mass" ) )
    assert 807000.0 == pytest.approx( data.AccessGCData( "Mass" )[0] )
    assert 98400.0 == pytest.approx( data.AccessGCData( "Mass" )[1] )
    
    assert 2 == len( data.AccessGCData( "R_a" ) )
    assert 7.44 == pytest.approx( data.AccessGCData( "R_a" )[0] )
    assert 12.26 == pytest.approx( data.AccessGCData( "R_a" )[1] )
    
    assert 2 == len( data.AccessGCData( "R_p" ) )
    assert 5.47 == pytest.approx( data.AccessGCData( "R_p" )[0] )
    assert 2.01 == pytest.approx( data.AccessGCData( "R_p" )[1] )
    
    assert 2 == len( data.AccessGCData( "SFE" ) )
    assert 0.3 == pytest.approx( data.AccessGCData( "SFE" )[0] )
    assert 0.3 == pytest.approx( data.AccessGCData( "SFE" )[1] )
    
    assert 2 == len( data.AccessGCData( "Fe-H" ) )
    assert -0.747 == pytest.approx( data.AccessGCData( "Fe-H" )[0] )
    assert -1.226 == pytest.approx( data.AccessGCData( "Fe-H" )[1] )
    
    assert 2 == len( data.AccessGCData( "FeSpread" ) )
    assert 0.033 == pytest.approx( data.AccessGCData( "FeSpread" )[0] )
    assert 0.037 == pytest.approx( data.AccessGCData( "FeSpread" )[1] )
    
    assert 2 == len( data.AccessGCData( "Age" ) )
    assert 12.0 == pytest.approx( data.AccessGCData( "Age" )[0] )
    assert 12.0 == pytest.approx( data.AccessGCData( "Age" )[1] )


def test_AccessGCDataPrinteable():
    """checks that the GC data is displayed correctly in printeable form"""
    
    data = SetupTest()
    
    ResultString = "Name       Mass        R_a      R_p     SFE    Fe-H      FeSpread                Age     \n"\
                 + "47_Tuc     807000.0    7.44     5.47    0.3    -0.747    0.033                   12.0    \n"\
                 + "NGC_288    98400.0     12.26    2.01    0.3    -1.226    0.037000000000000005    12.0    \n"
    
    assert ResultString == data.AccessGCDataPrinteable()


def test_AddGCData():
    """tests adding data to the GC dataset"""
    
    data = SetupTest()
    
    data.AddGCData( "NewColumn", [3,4] )
    
    assert 2 == len( data.AccessGCData( "NewColumn" ) )
    assert 3 == pytest.approx( data.AccessGCData( "NewColumn" )[0] )
    assert 4 == pytest.approx( data.AccessGCData( "NewColumn" )[1] )
    
    #if data of the wrong size shall be added an error should be raised
    with pytest.raises(ValueError, match=r"cData: Attempt to add Column '.*.' to GCData failed. Length of data does not match!"):
        data.AddGCData( "FaultyData", [1,2,3] )
        
        
def test_SNExplodes():
    """tests that the masses for which SN explode are determined correctly"""
    
    data = SetupTest()
    
    assert True == data.SNExplodes( 10.0 )
    assert True == data.SNExplodes( 12.3 )
    assert False == data.SNExplodes( 15.03 )
    assert False == data.SNExplodes( 15.11 )
    assert True == data.SNExplodes( 25.21 )
    assert True == data.SNExplodes( 140.0 )
    
    
def test_Ejecta():
    """tests that the iron ejecta are computed correctly"""
    
    #test with multiple values in ejecta and SN files
    data = SetupTest()
    
    assert 0.0 == pytest.approx( data.Ejecta( 7.0 ) )
    assert 0.002078521939953809 == pytest.approx( data.Ejecta( 9.180887372013653 ) )
    assert 0.006235565819861427 == pytest.approx( data.Ejecta( 9.399317406143345 ) )
    assert 0.041916859122401853 == pytest.approx( data.Ejecta( 12.034129692832765 ) )
    assert 0.060508083140877605 == pytest.approx( data.Ejecta( 13.358361774744028 ) )
    assert 0.0 == pytest.approx( data.Ejecta( 21.604095563139936 ) )
    assert 0.12222222222222223 == pytest.approx( data.Ejecta( 25.194444444444443 ) )
    assert 0.04791666666666668 == pytest.approx( data.Ejecta( 120.0 ) )
    
    #test with only one ejecta value given
    
    data1ejecta = cData( "test/mockdata/GCData.dat", "test/mockdata/SNe.dat", "test/mockdata/1Ejecta.dat", "test/mockdata/RemnantData.dat" )
    
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 7.0 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 9.180887372013653 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 9.399317406143345 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 12.034129692832765 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 13.358361774744028 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 21.604095563139936 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 25.194444444444443 ) )
    assert 0.074 == pytest.approx( data1ejecta.Ejecta( 120.0 ) )


def test_GetRemnantData():
    """tests whether the remnant data is read in correctly"""
    
    data = SetupTest()
    
    RemnantData = data.GetRemnantData()
    
    assert -1.0969100130080565 == pytest.approx( RemnantData["mass[Msun]"][0] )
    assert -0.8804142250382161 == pytest.approx( RemnantData["mass[Msun]"][1] )
    
    assert 14.781965892799771 == pytest.approx( RemnantData["t_-1.67"][0] )
    assert 13.68137744911132 == pytest.approx( RemnantData["t_-1.67"][1] )
    
    assert 14.500836618524938 == pytest.approx( RemnantData["t_-0.67"][0] )
    assert 13.450735279577588 == pytest.approx( RemnantData["t_-0.67"][1] )
    
    assert 14.278246979261908 == pytest.approx( RemnantData["t_-0.37"][0] )
    assert 13.26544876054881 == pytest.approx( RemnantData["t_-0.37"][1] )
    
    assert 14.082641523009046 == pytest.approx( RemnantData["t_-0.2"][0] )
    assert 13.09520930625935 == pytest.approx( RemnantData["t_-0.2"][1] )
    
    assert -1.0683387603155188 == pytest.approx( RemnantData["Mfin_-1.67"][0] )
    assert -0.9051796196452 == pytest.approx( RemnantData["Mfin_-1.67"][1] )
    
    assert -1.199145508496439 == pytest.approx( RemnantData["Mfin_-0.67"][0] )
    assert -1.023741750742955 == pytest.approx( RemnantData["Mfin_-0.67"][1] )
    
    assert -1.2547690154718592 == pytest.approx( RemnantData["Mfin_-0.37"][0] )
    assert -1.073091898294564 == pytest.approx( RemnantData["Mfin_-0.37"][1] )
    
    assert -1.2964507017617695 == pytest.approx( RemnantData["Mfin_-0.2"][0] )
    assert -1.111202932543319 == pytest.approx( RemnantData["Mfin_-0.2"][1] )


def test_ConstructorSameFile():
    """tests whether or not the data is read in correctly of the same file is used for all values"""
    
    CombinedData = cData( "test/mockdata/GCData.dat", "test/mockdata/CombinedData.dat", "test/mockdata/CombinedData.dat", "test/mockdata/CombinedData.dat" )
    IndividualData = cData( "test/mockdata/GCData.dat", "test/mockdata/2SNe.dat", "test/mockdata/2Ejecta.dat", "test/mockdata/RemnantData.dat" )
    
    
    #test that the ejecta data and SNe data are treated correctly
    for mass in [0.08, 0.1, 0.5, 1.0]:
        assert CombinedData.Ejecta( mass ) == IndividualData.Ejecta( mass )
        assert CombinedData.SNExplodes( mass ) == IndividualData.SNExplodes( mass )
    
    #test that the remnant data is treated correctly
    CombinedRemData = CombinedData.GetRemnantData()
    IndRemnantData = IndividualData.GetRemnantData()
    
    for Elem in CombinedRemData:
        if not Elem == "mass[Msun]" or Elem[:2] == "t_" or Elem[:5] == "Mfin_":
            continue
        
        for nElem in range( len( CombinedRemData["mass[Msun]"] ) ):
            assert CombinedRemData[Elem][nElem] == IndRemnantData[Elem][nElem]
