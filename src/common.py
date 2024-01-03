def IsNumber( str ):
    """checks, whether a string is a floating point number
    str: the string to check
    returns True if the string is a number and False otherwise"""
    
    try:
        float( str )
        return True
    except:
        return False


def LinInterExtrapolate( X1Y1, X2Y2, X3 ):
    """inter- or extrapolates linearily using to points of a linear function y(x)
    X1Y1: the first point in the form (X1,Y1)
    X2Y2: the second point in the form (X2Y2)
    X3: the x-value for which y shall be found
    returns Y3"""
    
    if X3 == X1Y1[0]:
        return X1Y1[1]
    if X3 == X2Y2[0]:
        return X2Y2[1]
    
    return ( X2Y2[1] - X1Y1[1] ) * ( X3 - X1Y1[0] ) / ( X2Y2[0] - X1Y1[0] ) + X1Y1[1]
