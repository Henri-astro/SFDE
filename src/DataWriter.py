class cDataWriter():
    """a class to write all the output data into dedicated output files"""
    
    def __init__( self, OutputFolder ):
        """OutputFolder: the folder the data shall be written to"""
        
        self.__OutputFolder = OutputFolder
    
    
    def WriteGCData( self, data ):
        """writes out all of the GCData using the data class
        data: the object in which all of the data is stored"""
        
        GCFile = open( self.__OutputFolder + "/GCData.dat" )
        
        GCFile.write( data.AccessGCDataPrinteable() )
        
        GCFile.close()
