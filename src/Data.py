from src.DataReader import cDataReader


class cData:
    """A class to hold all the data used by the program"""
    
    def __init__( self, GCFile, SNFile, EjectaFile ):
        """Constructor
        GCFile: a file that holds general information about the GCs
        SNFile: a file that holds information about which stars explode in SNe
        EjectaFile: a file that contains the amount of iron produced by exploding stars depending on their mass"""
        
        #import the GC data
        GCData = cDataReader( GCFile, ["Name", "Mass", "R_a", "R_p", "SFE", "Fe-H", "FeSpread"] )
        self.__GCData = GCData.GetData()
        
        #import the SN data
        SNData = cDataReader( SNFile, ["mass[Msun]", "SN"] )
        self.__SNData = SNData.GetData()
        
        #make sure the data is not empty
        if 0 == len( self.__SNData["mass[Msun]"] ):
            raise ValueError( "cData: Empty SN data. Please check your SN file." )
        
        #import the ejecta data
        EjectaData = cDataReader( EjectaFile, ["mass[Msun]", "Fe[Msun]"] )
        self.__EjectaData = EjectaData.GetData()
        
        #make sure the data is not empty
        if 0 == len( self.__EjectaData["mass[Msun]"] ):
            raise ValueError( "cData: Empty Ejecta data. Please check your Ejecta file." )
        
        
    def AccessGCData( self, ColumnName ):
        """return a single column from the GC data
        ColumnName: the name of the column from which the data shall be returned"""
        
        return self.__GCData[ ColumnName ]
    
    
    def AddGCData( self, ColumnName, data ):
        """adds a single column to the GCdata
        ColumnName: the name of the column to be added
        data: a list of the data to be added, the length of the list needs to match the length of the other columns"""
        
        if not len( data ) == len( self.__GCData["Name"] ):
            raise ValueError( "cData: Attempt to add Column '" + ColumnName + "' to GCData failed. Length of data does not match!" )
        
        self.__GCData[ColumnName] = data
    
    
    def SNExplodes( self, mass ):
        """returns whether or not a star of a certain mass is going to explode in a SNe
        mass: the mass of the star
        returns True if the star explodes and False if it doesn't"""
        
        #find the SN value corresponding to the given mass in the table
        if mass <= self.__SNData["mass[Msun]"][0]:
            return bool( self.__SNData["SN"][0] )
        
        if mass > self.__SNData["mass[Msun]"][-1]:
            return bool( self.__SNData["SN"][-1] )

        for nIndex in range( len( self.__SNData["mass[Msun]"] ) - 1 ):
            if mass <= self.__SNData["mass[Msun]"][nIndex + 1]:
                if self.__SNData["mass[Msun]"][nIndex + 1] - mass < mass - self.__SNData["mass[Msun]"][nIndex]:
                    return bool( self.__SNData["SN"][nIndex + 1] )
                else:
                    return bool( self.__SNData["SN"][nIndex] )
                    
    
    def Ejecta( self, mass ):
        """returns the amount of iron produced by a star of the given mass. This function assumes all stars explode.
        mass: the mass of the star"""
        
        #if only one value is known return this one value
        if 1 == len( self.__EjectaData["mass[Msun]"] ):
            Ejecta = self.__EjectaData["Fe[Msun]"][0]
        
        #check the ends of the ejecta table and extrapolate if neccesary
        elif mass == self.__EjectaData["mass[Msun]"][0]:
            Ejecta = self.__EjectaData["Fe[Msun]"][0]
        
        elif mass < self.__EjectaData["mass[Msun]"][0]:
            Ejecta = ( self.__EjectaData["Fe[Msun]"][0] - self.__EjectaData["Fe[Msun]"][1] ) / ( self.__EjectaData["mass[Msun]"][0] - self.__EjectaData["mass[Msun]"][1] ) * ( mass - self.__EjectaData["mass[Msun]"][0] ) + self.__EjectaData["Fe[Msun]"][0]
        
        elif mass > self.__EjectaData["mass[Msun]"][-1]:
            Ejecta = ( self.__EjectaData["Fe[Msun]"][-1] - self.__EjectaData["Fe[Msun]"][-2] ) / ( self.__EjectaData["mass[Msun]"][-1] - self.__EjectaData["mass[Msun]"][-2] ) * ( mass - self.__EjectaData["mass[Msun]"][-1] ) + self.__EjectaData["Fe[Msun]"][-1]
        
        else:
            for nIndex in range( len( self.__EjectaData["mass[Msun]"] ) - 1 ):
                if mass == self.__EjectaData["mass[Msun]"][nIndex + 1]:
                    Ejecta = self.__EjectaData["Fe[Msun]"][nIndex + 1]
                    
                    break
                
                if mass < self.__EjectaData["mass[Msun]"][nIndex + 1]:
                    Ejecta = ( self.__EjectaData["Fe[Msun]"][nIndex] - self.__EjectaData["Fe[Msun]"][nIndex + 1] ) / ( self.__EjectaData["mass[Msun]"][nIndex] - self.__EjectaData["mass[Msun]"][nIndex + 1] ) * ( mass - self.__EjectaData["mass[Msun]"][nIndex] ) + self.__EjectaData["Fe[Msun]"][nIndex]
                    
                    break
                
        #the amount of iron ejected cannot be negative
        if Ejecta < 0.0:
            Ejecta = 0.0
            
        return Ejecta
