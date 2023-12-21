import pytest
import pandas as pd

from src.DataReader import cDataReader


def CompareData( data ):
    """Helper function for test_data. Compares the given data to the expected mock data."""
    
    assert data["Column1"][0] == 1.9
    assert data["Column1"][1] == 9.4
    assert data["Column2"][0] == 3.1
    assert data["Column2"][1] == 1.3
    assert data["Column3"][0] == "bla"
    assert data["Column3"][1] == "keks"


def test_data():
    """tests if the data is read in correctly based on the mock file"""
    
    data = cDataReader( "test/mockdata/TestData.dat", ["Column1","Column2","Column3"] )
    
    CompareData( data.GetPandasSheet() )
    CompareData( data.GetData() )


def test_Exceptions():
    """tests if the exceptions in this class are raised appropriately"""
    
    with pytest.raises(NameError, match=r"cDataReader: Error importing data. Column '.*.' missing!"):
        data = cDataReader( "test/mockdata/TestData.dat", ["Column1","Column2","Column3","Column4"] )
