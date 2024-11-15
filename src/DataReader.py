#general libs
import pandas as pd


class cDataReader:
    """A class to read in the given datafiles"""
    
    def __init__( self, FileName, ExpectedColumns ):
        """Constructor
        FileName: The name of the file that contains the GC information
        ExpectedColumns: column headers that need to be in the input file. If one column header is missing a NameError is raised
        The data is expected to be in coloums seperated by spaces."""
        
        self.ReadFile( FileName, ExpectedColumns )
    
    
    def ReadFile( self, FileName, ExpectedColumns ):
        """reads in the data from the given file
        FileName: the file that contains the data"""
        
        #read in the data
        self.__datasheet = pd.read_table( FileName, sep='\\s+', comment='#' )
        
        #check that all columns are present        
        header = self.__datasheet.columns.values.tolist()
        
        for Name in ExpectedColumns:
            if not Name in header:
                raise NameError( "cDataReader: Error importing data. Column '" + Name + "' missing!" )
        
    
    def GetPandasSheet( self ):
        """returns the pandas datasheet"""
        
        return self.__datasheet
    
    
    def GetData( self ):
        """returns the data in form of a plain python list"""
        
        return { header: tuple( self.__datasheet[header] ) for header in self.__datasheet.columns.values.tolist() }
