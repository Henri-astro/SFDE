#general modules
import numpy as np

#own module
from src.common import IsNumber
from src.DataReader import cDataReader


class cData:
    """A class to hold all the data used by the program"""
    
    def __init__( self, GCFile, SNFile, EjectaFile, RemnantFile ):
        """Constructor
        GCFile: a file that holds general information about the GCs
        SNFile: a file that holds information about which stars explode in SNe
        EjectaFile: a file that contains the amount of iron produced by exploding stars depending on their mass"""
        
        #import the GC data
        GCData = cDataReader( GCFile, ["Name", "Mass", "R_a", "R_p", "SFE", "Fe-H", "FeSpread", "Age"] )
        self.__GCData = GCData.GetData()
        
        #if the last three files are the same I only need to read in the data once
        if SNFile == EjectaFile == RemnantFile:
            StarData = cDataReader( SNFile, ["mass[Msun]", "SN", "Fe[Msun]"] )

            self.__SNData = StarData.GetData()
            self.__EjectaData = self.__SNData
            
        else:
            #import the SN data
            SNData = cDataReader( SNFile, ["mass[Msun]", "SN"] )
            self.__SNData = SNData.GetData()
        
            #import the ejecta data
            EjectaData = cDataReader( EjectaFile, ["mass[Msun]", "Fe[Msun]"] )
            self.__EjectaData = EjectaData.GetData()
        
        #import the remnant data
        RemnantData = cDataReader( RemnantFile, ["mass[Msun]"] )
        self.__RemnantData = RemnantData.GetData()
        
        #make sure the SN data is not empty
        if 0 == len( self.__SNData["mass[Msun]"] ):
            raise ValueError( "cData: Empty SN data. Please check your SN file." )
        
        #make sure the ejecta data is not empty
        if 0 == len( self.__EjectaData["mass[Msun]"] ):
            raise ValueError( "cData: Empty Ejecta data. Please check your Ejecta file." )
        
        #check that at least one column with t and Mfin is present in the remnant data
        t_present = False
        Mfin_present = False
        
        for Heading in self.__RemnantData:
            if Heading[:2] == "t_":
                if IsNumber( Heading[2:] ):
                    t_present = True
                
            if Heading[:5] == "Mfin_":
                if IsNumber( Heading[5:] ):
                    Mfin_present = True
                        
        if not t_present:
            raise ValueError( "cData: time information in Remnant data missing. Please check your Remnant file." )
        
        if not Mfin_present:
            raise ValueError( "cData: final masses in Remnant data missing. Please check your Remnant file." )
        
        #make sure the data is not empty
        if 0 == len( self.__RemnantData["mass[Msun]"] ):
            raise ValueError( "cData: Empty Remnant data. Please check your Remnant file." )
        
        #change all the initial and remnant masses into their logarithms
        self.__RemnantData["mass[Msun]"] = tuple( np.log10( logmass ) for logmass in self.__RemnantData["mass[Msun]"] )
        
        for Elem in self.__RemnantData:
            if Elem[:5] == "Mfin_":
                self.__RemnantData[Elem] = tuple( np.log10( logmass ) for logmass in self.__RemnantData[Elem] )
        
        
    def AccessGCData( self, ColumnName ):
        """return a single column from the GC data
        ColumnName: the name of the column from which the data shall be returned"""
        
        return self.__GCData[ ColumnName ]
    
    
    def AccessGCDataPrinteable( self ):
        """returns a string of the GC data for printing"""
        
        Header = ""
        lines = [ "" for Elem in self.__GCData["Name"] ]
        
        for Name in self.__GCData:
            if not type(self.__GCData[ Name ][0]  ) in [ str, int, float, np.float64 ]:
                continue
            
            #find maximum length
            length = len( str( Name ) )
            for Elem in self.__GCData[ Name ]:
                newlen = len( str( Elem ) )
                
                if newlen > length:
                    length = newlen
            
            length += 4
            
            Header += Name.ljust(length)
            
            for nElem in range( len( self.__GCData[ Name ] ) ):
                lines[nElem] += str( self.__GCData[ Name ][nElem] ).ljust( length )
                
        ReturnString = Header
        
        for line in lines:
            ReturnString += "\n"
            ReturnString += line
            
        return ReturnString
    
    
    def AddGCData( self, ColumnName, data ):
        """adds a single column to the GCdata
        ColumnName: the name of the column to be added
        data: a tuple of the data to be added, the length of the tuple needs to match the length of the other columns"""
        
        if not len( data ) == len( self.__GCData["Name"] ):
            raise ValueError( "cData: Attempt to add Column '" + ColumnName + "' to GCData failed. Length of data does not match!" )
        
        self.__GCData[ColumnName] = tuple( data )
    
    
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
        mass: the mass of the star
        returns the amount of iron produced [Msun]"""
        
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
    
    
    def GetRemnantData( self ):
        """returns a copy of the complete remnant data"""
        
        return self.__RemnantData.copy()
