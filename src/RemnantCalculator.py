#general modules
import numpy as np

#own modules
from src.common import IsNumber, LinInterExtrapolate


class cRemnantCalculator():
    """a class that computes the remnant mass and life time of a given star"""
    
    def __init__( self, data, ZH ):
        """Constructor
        data: a lib containing initial masses, development times and remnant masses for stars, the masses are expected to be given in accending order
        ZH: the metallicity for which the remnant cRemnantCalculator shall be used"""
        
        self.__ZH = ZH
        
        #copy the initial masses of the stars
        self.__Mstar = data["mass[Msun]"]
        self.__t = self.__FindValues( data, ZH, "t_" )
        self.__Mfin = self.__FindValues( data, ZH, "Mfin_" )
    
    
    def __FindValues( self, data, ZH, Name ):
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
                return data[ Heading ]   #if it's the same, no interpolation needed
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
                return data[ Above["ColumnName"] ]
            else:
                Below = Above2.copy()   #interpolation and extrapolation can be done in the same way
                
        if None == Above["ColumnName"]:
            if None == Below2["ColumnName"]:
                return data[ Below["ColumnName"] ]
            else:
                Above = Below2.copy()   #interpolation and extrapolation can be done in the same way
                
        #make the final list:
        return [ LinInterExtrapolate( ( pow( 10.0, Below["metal"] ), data[ Below["ColumnName"] ][nElem] ), ( pow( 10.0, Above["metal"] ), data[ Above["ColumnName"] ][nElem] ), pow( 10.0, ZH ) ) for nElem in range( len( self.__Mstar ) ) ]

        
    def GetZH( self ):
        return self.__ZH
    

    def __GetValFromList( self, val, valList, dataList ):
        """searches for the value corresponding to a given value in a given valueList and returns the corresponding value from the dataList (treating it as a function data( val ))
        val: the mass of a star for which a property shall be found [log(Msun)]
        valList: a list of possible values sorted from smallest to largest
        dataList: the list to retrieve the value from
        returns the value corresponding to the given mass, the mass cannot be found exactly it is interpolated/extrapolated"""
        
        #make sure to handle masses smaller than the smallest element correctly
        if val < valList[1]:
            return LinInterExtrapolate( ( valList[0], dataList[0] ),( valList[1], dataList[1] ), val )
        
        #make sure to handle masses larger than the largest element correctly
        if val > valList[-2]:
            return LinInterExtrapolate( ( valList[-1], dataList[-1] ),( valList[-2], dataList[-2] ), val )
        
        lenValList = len( valList )
        pos = int( lenValList / 2 )
        
        for i in range( int( np.log2( pos ) ) + 1 ):
            if val == valList[pos]:
                return dataList[pos]
            elif valList[pos] < val:
                if valList[pos + 1] > val:
                    break
                else:
                    pos += int( lenValList / pow( 2, i + 2 ))
            else:
                if valList[pos - 1] < val:
                    pos = pos - 1
                    break
                else:
                    pos -= int( lenValList / pow( 2, i + 2 ))
        
        return LinInterExtrapolate( ( valList[pos], dataList[pos] ), ( valList[pos + 1], dataList[pos + 1] ), val )
    
    
    def __GetValFromReversedList( self, val, valList, dataList ):
        """searches for the value corresponding to a given value in a given valueList and returns the corresponding value from the dataList (treating it as a function data( val ))
        val: the mass of a star for which a property shall be found [log(Msun)]
        valList: a list of possible values sorted from largest to smallest
        dataList: the list to retrieve the value from
        returns the value corresponding to the given mass, the mass cannot be found exactly it is interpolated/extrapolated"""
        
        #make sure to handle masses larger than the largest element correctly
        if val > valList[1]:
            return LinInterExtrapolate( ( valList[0], dataList[0] ),( valList[1], dataList[1] ), val )
        
        #make sure to handle masses smaller than the smallest element correctly
        if val < valList[-2]:
            return LinInterExtrapolate( ( valList[-1], dataList[-1] ),( valList[-2], dataList[-2] ), val )
        
        lenValList = len( valList )
        pos = int( lenValList / 2 )
        
        for i in range( int( np.log2( pos ) ) + 1 ):
            if val == valList[pos]:
                return dataList[pos]
            elif valList[pos] < val:
                if valList[pos - 1] > val:
                    pos = pos - 1
                    break
                else:
                    pos -= int( lenValList / pow( 2, i + 2 ))
            else:
                if valList[pos + 1] < val:
                    break
                else:
                    pos += int( lenValList / pow( 2, i + 2 ))
        
        return LinInterExtrapolate( ( valList[pos], dataList[pos] ), ( valList[pos + 1], dataList[pos + 1] ), val )
    
            
    def GetTimeFromMass( self, mass ):
        """returns the life-expectancy of a star of a given mass
        mass: the mass of the star [Msun]
        returns the life expectancy of the star [Gyr]"""
        
        logTimeYr = self.__GetValFromList( mass, self.__Mstar, self.__t )
        
        return pow( 10.0, logTimeYr - 9.0 )
    
    
    def GetMassFromTime( self, t ):
        """returns the mass of a star corresponding to the given life-expectancy
        time: the life-expectancy of the star [Gyr]
        returns the mass of the star [Msun]"""
        
        logTimeYr = np.log10( t ) + 9.0
        
        return self.__GetValFromReversedList( logTimeYr, self.__t, self.__Mstar )
    
    
    def GetMfinFromLogMass( self, logMass ):
        """computes the mass of stellar remnant from its initial mass
        mass: the mass of a star log[Msun]
        returns the mass of the stellar remnant [Msun]"""
        
        logMfin = self.__GetValFromList( logMass, self.__Mstar, self.__Mfin )
        
        return pow( 10.0, logMfin )
    
    
    def GetMfinFromMass( self, mass ):
        """computes the mass of stellar remnant from its initial mass
        mass: the mass of a star [Msun]
        returns the mass of the stellar remnant [Msun]"""
        
        return self.__GetValFromList( mass, self.__Mstar, self.__Mfin )


    def GetMfinFromMassFunct( self, MF, time ):
        """computes the mass left in a cluster of a given mass function after a given time
        MF: the mass function describing the cluster [cMassFunction]
        time: the time for which the left-over mass shall be computed [Gyr]
        N: number of sampling points to be used"""
        
        #compute the lowest mass of a star to have died
        minMass = self.GetMassFromTime( time )
        
        if minMass > MF.Getbounds()[-1]:
            return MF.GetMtot()
        
        maxMass = MF.Getbounds()[-1]
        
        #compute the current SC mass:
        CurSCMass = MF.GetMass( MF.Getbounds()[0], minMass )  #the stars below minMass contribute with their entire mass
        
        #for the other stars only the remnant masses contribute
        MassPart = lambda LowBound, HighBound: self.GetMfinFromMass( 0.5 * ( LowBound + HighBound ) ) * MF.GetMass( LowBound, HighBound )
        
        #boundary cases
        if minMass > self.__Mstar[-1] or maxMass < self.__Mstar[0]:
            CurSCMass += MassPart( minMass, maxMass ) / MF.GetMtot()
            return CurSCMass
        
        #normal cases
        Start = 0
        Finnish = -1
        
        for nElem in range( len( self.__Mstar ) ):
            if minMass < self.__Mstar[ nElem ]:
                Start = nElem
                break
            
        for nElem in range( len( self.__Mstar ) - 1, -1, -1 ):
            if maxMass > self.__Mstar[ nElem ]:
                Finnish = nElem
                break
        
        if Start == Finnish or Finnish < Start:
            CurSCMass += MassPart( minMass, maxMass ) / MF.GetMtot()
            return CurSCMass
        
        CompRemnantMass = 0.0       #add the remnant mass separately first (so I only have to do one division at the end)
        
        CompRemnantMass += MassPart( minMass, 0.5 * ( self.__Mstar[ Start ] + self.__Mstar[ Start + 1 ] ) )
        
        for nElem in range( Start + 1, Finnish ):
            CompRemnantMass += MassPart( 0.5 * ( self.__Mstar[ nElem - 1 ] + self.__Mstar[ nElem ] ), 0.5 * ( self.__Mstar[ nElem ] + self.__Mstar[ nElem + 1 ] ) )
            
        CompRemnantMass += MassPart( 0.5 * ( self.__Mstar[ Finnish - 1 ] + self.__Mstar[ Finnish ] ), maxMass )
        
        CurSCMass += CompRemnantMass / MF.GetMtot()
            
        return CurSCMass
