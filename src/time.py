import numpy as np

from src.timedata import Mini, T_SN, Mfin


def FindClosestInList( val, List ):
    """finds the index of the element closest to val in the list
    val: the value
    List: the list
    returns the index of the element closest to val"""
    
    diff = abs( val - List[0] )
    index = 0
    
    for nElem in range( 1, len( List ) ):
        newDiff = abs( val - List[nElem] )
        if newDiff < diff:
            index = nElem
            diff = newDiff
    
    return index


def MassToTime( mass ):
    """gives the time at which a star of the given initial mass dies
    mass: stellar mass [Msun]
    returns the time [Myr]"""
    
    logmass = np.log10( mass )
    
    for nElem in range( len( Mini ) ):
        if logmass < Mini[nElem]:
            if nElem == 0 or ( logmass - Mini[nElem - 1] > Mini[nElem] - logmass ):
                return pow( 10.0, T_SN[nElem] - 6.0 )
            else:
                return pow( 10.0, T_SN[nElem - 1] - 6.0 )
        
    return pow( 10.0, T_SN[-1] - 6.0 )


def TimeToMass( time ):
    """gives mass of a star that dies at a given time
    time: the time [Myr]
    returns the stellar mass [Msun]
    WARNING: since time(Mini) is not injective, this function may yield wrong results"""
    
    logtime = np.log10( time ) + 6 #the time in the timedata is log10( yr )
    
    for nElem in range( len( T_SN ) ):
        if logtime > T_SN[nElem]:
            if nElem == 0 or ( logtime - T_SN[nElem] < T_SN[nElem-1] - logtime ):
                return pow( 10.0, Mini[nElem] )
            else:
                return pow( 10.0, Mini[nElem-1] )
        
    return pow( 10.0, Mini[-1] )


def MassToMfin( mass ):
    """gives the remnant mass for a star of the given initial mass
    mass: stellar mass [Msun]
    returns the time [Myr]"""
    
    logmass = np.log10( mass )
    
    for nElem in range( len( Mini ) ):
        if logmass < Mini[nElem]:
            if nElem == 0 or ( logmass - Mini[nElem - 1] > Mini[nElem] - logmass ):
                return pow( 10.0, Mfin[nElem] )
            else:
                return pow( 10.0, Mfin[nElem - 1] )
        
    return pow( 10.0, Mfin[-1] )
