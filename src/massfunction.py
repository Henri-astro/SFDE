import sys
import numpy as np


class cMassFunction:

    def __init__( self, Mtot, bounds, alphas ):
        """instantiation of the class
        self: the "this" equivalent
        Mtot: the total mass of the cluster in question
        bounds: a list of m_min, boundaries between segments and m_max, it is assumed that m_max is the mass of the most massive star in the SC, therefore M_SC = int limits_{m_min}^{m_max} m xi(m) dm + m_max
        alphas: a list of alphas (there always needs to be one alpha more then boundaries)"""
        
        #default for faulty inputs
        if len( bounds ) != len( alphas ) + 1:
            sys.stderr.write( "IMF __init__: Invalid number of alphas, must be one less then boundaries!\nDefault to Kroupa 2001 IMF." )
            self.__bounds = [0.08,0.5,150.0]
            self.__alphas = [1.3,2.3]
        else:
            self.__bounds = bounds.copy()
            self.__alphas = alphas.copy()
        
        #transfer the parameters
        self.__Mtot = Mtot
        
        #compute the k-values for the mass function
        self.__ComputeKs()
        
        
    def __ComputeKs( self ):
        """compute the k values from the given masses and alphas"""
        
        #compute preliminary k's
        self.__ks = [1.0]
        
        for nAlpha in range(len(self.__alphas) - 1):
            self.__ks.append( self.__ks[-1] * pow( self.__bounds[nAlpha + 1], -self.__alphas[nAlpha] )/pow( self.__bounds[nAlpha + 1], -self.__alphas[nAlpha + 1] ))
        
        #correct the ks
        Masses = []
        
        for nAlpha in range(len(self.__alphas)):
            if 2.0 == self.__alphas[nAlpha]:
                Masses.append(self.__ks[nAlpha] * ( np.log( self.__bounds[nAlpha + 1] ) - np.log( self.__bounds[nAlpha] )))
                continue
            
            Masses.append(self.__ks[nAlpha] * ( pow( self.__bounds[nAlpha + 1], 2.0 - self.__alphas[nAlpha] ) - pow( self.__bounds[nAlpha], 2.0 - self.__alphas[nAlpha] )) / ( 2.0 - self.__alphas[nAlpha] ) )
        
        MassTot = sum( Masses )
        
        Factor = ( self.__Mtot - self.__bounds[-1] ) / MassTot            # Mtot = \int\limits_{m_min}^{m_max} m \xi( m ) dm + m_max
        self.__ks = [ Factor * Elem for Elem in self.__ks ]
    
    
    def GetMtot( self ):
        return self.__Mtot
    
    
    def Getbounds( self ):
        return self.__bounds.copy()
    
    
    def Getalphas( self ):
        return self.__alphas.copy()
    
    
    def Getks( self ):
        return self.__ks.copy()


    def GetFunctionValue( self, mass ):
        """gets the function value for a given mass
        self: the "this" equivalent
        mass: the mass at which the function value shall be returned
        throws a ValueError if the mass is out of bounds"""
        
        if mass < self.__bounds[0] or mass > self.__bounds[-1]:
            raise ValueError( "massfunction: GetFunctionValue: mass out of bounds!" )
        
        for nElem in range( len( self.__bounds ) - 1 ):
            if mass <= self.__bounds[nElem+1]:
                return self.__ks[nElem] * pow( mass, -self.__alphas[nElem] )
            
    
    def GetMassDensity( self, mass ):
        """gets the mass density (function * mass) value for a given mass
        self: the "this" equivalent
        mass: the mass at which the function value shall be returned
        throws a ValueError if the mass is out of bounds"""
        
        if mass < self.__bounds[0] or mass > self.__bounds[-1]:
            raise ValueError( "massfunction: GetMassDensity: mass out of bounds!" )
        
        for nElem in range( len( self.__bounds ) - 1 ):
            if mass <= self.__bounds[nElem+1]:
                return mass * self.__ks[nElem] * pow( mass, -self.__alphas[nElem] )
    
    
    def GetAlpha( self, Index ):
        """gets the stored alpha for the given index
        self: the "this" equivalent
        Index: the index at which alpha shall be found
        throws a ValueError if the index is out of bounds"""
        
        if Index >= len( self.__alphas ) or 0 > Index:
            raise ValueError( "massfunction: GetAlpha: index out of bounds!" )
        
        return self.__alphas[Index]
    
    
    def GetK( self, Index ):
        """gets the stored k for the given index
        self: the "this" equivalent
        Index: the index at which k shall be found
        throws a ValueError if the index is out of bounds"""
        
        if Index >= len( self.__ks ) or 0 > Index:
            raise ValueError( "massfunction: GetK: index out of bounds!" )
        
        return self.__ks[Index]
    
    
    def ComputeIntegral( self, mass1, mass2, summand = 0 ):
        """computes the integral of the IMF adding summand to -alpha int limits_mass1^mass2 k_i m^{summand-alpha_i}
        mass1: the low-mass end of the mass interval to be investigated (the interval is truncated when outside of the boundaries of the MF)
        mass2: the high-mass end of the mass interval to be investigated (the interval is truncated when outside of the boundaries of the MF)
        summand: the summand added to the exponent -alpha (1 to compute the mass within an interval, 0 for the number of stars)"""
            
        IntRan = lambda m1, m2, alpha, k: k * (np.log(m2) - np.log(m1)) if summand + 1 == alpha else k / (1.0 + summand - alpha) * ( pow(m2, 1.0 + summand - alpha) - pow(m1, 1.0 + summand - alpha) )
        
        Int = 0
        InRange = False
        
        if mass1 < self.__bounds[0]:
            mass1 = self.__bounds[0]
            
        if mass2 > self.__bounds[-1]:
            mass2 = self.__bounds[-1]
    
        for nBound in range( len( self.__bounds ) - 1 ):
            if self.__bounds[nBound] <= mass1 and mass1 < self.__bounds[nBound + 1]:
                if self.__bounds[nBound] < mass2 and mass2 <= self.__bounds[nBound + 1]:
                    Int = IntRan( mass1, mass2, self.__alphas[nBound], self.__ks[nBound] )
                    break
                
                Int += IntRan( mass1, self.__bounds[nBound + 1], self.__alphas[nBound], self.__ks[nBound] )
                InRange = True
                continue
            
            if InRange == True:
                if self.__bounds[nBound] < mass2 and mass2 <= self.__bounds[nBound + 1]:
                    Int += IntRan( self.__bounds[nBound], mass2, self.__alphas[nBound], self.__ks[nBound] )
                    InRange = False
                    break
                else:
                    Int += IntRan( self.__bounds[nBound], self.__bounds[nBound + 1], self.__alphas[nBound], self.__ks[nBound] )
                    
        return Int
    
    
    def GetNumbers( self, mass1, mass2 ):
        """computes what number of stars is between mass1 and mass2
        mass1: the low-mass end of the mass interval to be investigated
        mass2: the high-mass end of the mass interval to be investigated
        a closed interval [mass1,mass2] is expected"""
        
        maxMass = False
        
        if mass2 >= self.__bounds[-1]:
            maxMass = True
    
        Number = self.ComputeIntegral( mass1, mass2 )
        
        if maxMass:
            return Number + 1
                    
        return Number
    
    
    def GetTotNumbers( self ):
        """returns the total number of stars in the cluster"""
        
        return self.GetNumbers( self.__bounds[0], self.__bounds[-1] )
    
    
    def GetMass( self, mass1, mass2 ):
        """computes what mass is between mass1 and mass2
        mass1: the low-mass end of the mass interval to be investigated
        mass2: the high-mass end of the mass interval to be investigated
        a closed interval [mass1,mass2] is expected"""
        
        maxMass = False
        
        if mass2 >= self.__bounds[-1]:
            maxMass = True
    
        Mass = self.ComputeIntegral( mass1, mass2, 1 )
        
        if maxMass:
            return Mass + self.__bounds[-1]
                    
        return Mass
    
    
    def GetMassPortion( self, mass1, mass2 ):
        """computes what portion of the overall GC mass is between mass1 and mass2
        mass1: the low-mass end of the mass interval to be investigated
        mass2: the high-mass end of the mass interval to be investigated"""
        
        return self.GetMass( mass1, mass2 ) / self.__Mtot
    
    
    def GetMassStarMinX( self, mass, NumStars ):
        """returns the mass of teh Xth less massive star than a star of a given mass assuming optimal sampling
        mass: the mass of the given star
        NumStars: the number of stars to the next one, can be float"""
        
        #check the outer boundaries
        if mass > self.__bounds[-1] or mass < self.__bounds[0]:
            raise ValueError( "massfunction: GetMassStarMinX: mass out of bounds!" )
        
        #check how much above the last boundary the star is
        AlphaIndex = 0          #the index of the alpha used
        
        for nBound in range( len( self.__bounds ) - 1 ):
            if mass <= self.__bounds[ nBound + 1 ]:
                AlphaIndex = nBound
                break
        
        NumToBound = self.GetNumbers( self.__bounds[AlphaIndex], mass )
        
        if NumToBound < NumStars:
            return self.GetMassStarMinX( self.__bounds[AlphaIndex], NumStars - NumToBound )
        
        if self.__alphas[ AlphaIndex ] == 1.0:
            return np.exp( np.log( mass ) - NumStars / self.__ks[ AlphaIndex ] )
        else:
            return pow( pow( mass, 1.0 - self.__alphas[ AlphaIndex ] ) - NumStars * ( 1.0 - self.__alphas[ AlphaIndex ] ) / self.__ks[ AlphaIndex ], 1.0 / ( 1.0 - self.__alphas[ AlphaIndex ] ) )
