#general modules
import numpy as np
import copy

#own modules
from src.common import IsNumber, LinInterExtrapolate


class cRemnantCalculator():
    """a class that computes the remnant mass and life time of a given star"""
    
    def __init__( self, data, ZH ):
        """Constructor
        data: a lib containing initial masses, development times and remnant masses for stars, the masses are expected to be given in accending order
        ZH: the metallicity for which the remnant cRemnantCalculator shall be used"""
        
        #copy the initial masses of the stars
        self.__Mstar = data["Mstar"]
        self.__t = self.FindValues( data, ZH, "t_" )
        self.__Mfin = self.FindValues( data, ZH, "Mfin_" )
    
    
    def FindValues( self, data, ZH, Name ):
        """finds the closests columns for a Z-value and turns them into one list of values
        data: a lib containing initial masses, development times and remnant masses for stars
        ZH: the metallicity for which the remnant cRemnantCalculator shall be used
        Name: the beginning of the column names before the metallicity data
        returns a list of interpolated values"""
        
        #variable to store the closest columns in
        Below = { "ColumnName": None, "metal": None }
        Below2 = { "ColumnName": None, "metal": None }
        Above = { "ColumnName": None, "metal": None }
        Above2 = { "ColumnName": None, "metal": None }
        
        #find those closest columns
        lenName = len( Name )
        
        for Heading in data:
            if not ( Heading[:lenName] == Name ):
                continue
            if not IsNumber( Heading[lenName:] ):
                continue
            
            ColumnMetal = float( Heading[lenName:] )
            
            if ColumnMetal == ZH:
                return copy.deepcopy( data[ Heading ] )   #if it's the same, no interpolation needed
            elif ColumnMetal < ZH:
                if None == Below["ColumnName"] or ColumnMetal > Below["metal"]:
                    Below2 = Below.copy()
                    Below["ColumnName"] = Heading
                    Below["metal"] = ColumnMetal
            elif ColumnMetal > ZH:
                if None == Above["ColumnName"] or ColumnMetal < Above["metal"]:
                    Above2 = Above.copy()
                    Above["ColumnName"] = Heading
                    Above["metal"] = ColumnMetal
        
        #cover special cases (no metal above/below known) e.t.c.
        if None == Below["ColumnName"]:
            if None == Above2["ColumnName"]:
                return copy.deepcopy( data[ Above["ColumnName"] ] )
            else:
                Below = Above2.copy()   #interpolation and extrapolation can be done in the same way
                
        if None == Above["ColumnName"]:
            if None == Below2["ColumnName"]:
                return copy.deepcopy( data[ Below["ColumnName"] ] )
            else:
                Above = Below2.copy()   #interpolation and extrapolation can be done in the same way
                
        #make the final list:
        return [ LinInterExtrapolate( ( Below["metal"], data[ Below["ColumnName"] ][nElem] ), ( Above["metal"], data[ Above["ColumnName"] ][nElem] ), ZH ) for nElem in range( len( self.__Mstar ) ) ]


    def GetValFromList( self, val, valList, dataList ):
        """searches for the value corresponding to a given value in a given valueList and returns the corresponding value from the dataList (treating it as a function data( val ))
        val: the mass of a star for which a property shall be found [log(Msun)]
        dataList: the list to retrieve the value from
        returns the value corresponding to the given mass, the mass cannot be found exactly it is interpolated/extrapolated"""
        
        #make sure to handle masses smaller than the smallest element correctly
        if val < valList[0]:
            return LinInterExtrapolate( ( valList[0], dataList[0] ),( valList[1], dataList[1] ), val )
        
        #make sure to handle masses larger than the largest element correctly
        if val > valList[-1]:
            return LinInterExtrapolate( ( valList[-1], dataList[-1] ),( valList[-2], dataList[-2] ), val )
        
        for nElem in range( len( valList ) ):
            if val == valList[nElem]:
                return dataList[nElem]
            elif val < valList[nElem]:
                return LinInterExtrapolate( ( valList[nElem-1], dataList[nElem-1] ), ( valList[nElem], dataList[nElem] ), val )
            
            
    def GetTimeFromMass( self, mass ):
        """returns the life-expectancy of a star of a given mass
        mass: the mass of the star [Msun]
        returns the life expectancy of the star [Gyr]"""
        
        logMass = np.log10( mass )
        
        logTimeYr = self.GetValFromList( logMass, self.__Mstar, self.__t )
        
        return pow( 10.0, logTimeYr - 9.0 )
    
    
    def GetMassFromTime( self, t ):
        """returns the mass of a star corresponding to the given life-expectancy
        time: the life-expectancy of the star [Gyr]
        returns the mass of the star [Msun]"""
        
        logTimeYr = np.log10( t ) + 9.0
        
        logMass = self.GetValFromList( logTimeYr, self.__t, self.__Mstar )
        
        return pow( 10.0, logMass )
    
    
    def GetMfinFromMass( self, mass ):
        """computes the mass of stellar remnant from its initial mass
        mass: the mass of a star [Msun]
        returns the mass of the stellar remnant [Msun]"""
        
        logMass = np.log10( mass )
        
        logMfin = self.GetValFromList( logMass, self.__Mstar, self.__Mfin )
        
        return pow( 10.0, logMfin )


    def GetMfinFromMassFunct( self, MF, time, N = 100 ):
        """computes the mass left in a cluster of a given mass function after a given time
        MF: the mass function describing the cluster [cMassFunction]
        time: the time for which the left-over mass shall be computed [Gyr]
        N: number of sampling points to be used"""
        
        #compute the lowest mass of a star to have died
        minMass = GetMassFromTime( time )
        
        logMinMass = np.log10( minMass )
        logMaxMass = np.log10( MF.Getbounds[-1] )
        
        #compute the current SC mass:
        CurSCMass = MF.GetMass( MF.Getbounds[0], minMass )  #the stars below minMass contribute with their entire mass
        
        Step = ( logMaxMass - logMinMass ) / N
        
        for sample in range( N ):
            MidMass = pow( 10.0, logMinMass + ( sample + 0.5 ) * Step )
            CurSCMass += MF.GetMass( pow( 10.0, logMinMass + sample * Step ), pow( 10.0, logMinMass + ( sample + 1 ) * Step )  ) * GetMfinFromMass( MidMass ) / MidMass     #the mass enclosed in the step on the MF * the mass portion left after the star dies
            
        return CurSCMass
