class cStarExtractor():
    """a class to extract the masses of individual stars from a mass function"""
    
    def __init__( self, massfunction ):
        """Constructor
        massfunction: the mass function this extractor shall work on"""
        
        self.__massfunction = massfunction
        self.__lastmass = None
    
    
    def GetNextMostMassiveStar( self ):
        """returns the next most massive star"""
        
        if self.__lastmass == None:
            self.__lastmass = self.__massfunction.Getbounds()[-1]
        else:
            self.__lastmass = self.__massfunction.GetMassStarMinX( self.__lastmass, 1.0 )
            
        return self.__lastmass
        
