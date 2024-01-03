#general libs
import warnings
import numpy as np

#own packages
from src.massfunction import cMassFunction
from src.BaseSCCalc import ComputeDens


METAL_IN_SUN = 0.0142

class cIMFGenerator():
    """a class to generate initial mass functions based on a GC's current chemical and orbital properties"""
    
    def __init__( self, ZH ):
        """initialises the generator
        ZH: the metallicity of the GC in dex"""
        
        self.__ZH = ZH
    
    
    def ComputeAlpha( self, Mini ):
        """computes the alpha for a GC
        Mini: initial mass Msun
        returns a list of alphas"""
        
        #compute the alphas acording to Yan et al. (2021)
        alpha = [0.0,0.0,0.0]
        
        Dalpha = 63
        
        alpha[0] = 1.3 + Dalpha * ( pow( 10.0, self.__ZH ) - 1 ) * METAL_IN_SUN
        alpha[1] = 2.3 + Dalpha * ( pow( 10.0, self.__ZH ) - 1 ) * METAL_IN_SUN
        
        y = -0.14 * self.__ZH + 0.99 * np.log10( ComputeDens( Mini ) * pow( 10.0, -6.0 ) )
        
        if y < -0.87:
            alpha[2] = 2.3
        else:
            alpha[2] = -0.41 * y + 1.94
            
        return alpha
    
    
    def CheckUpperEnd( self, MF ):
        """checks that the integral between m_max and 150 Msun of the MF is 1
        MF: the mass function to check
        returns int_m_max^150Msun MF - 1.0
        This function assumes that there are no breakpoints in the IMF above m_max"""
        
        m_max = MF.Getbounds()[ -1 ]
        
        if MF.Getalphas()[-1] == 1.0:
            return MF.Getks()[-1] * ( np.log( 150.0 ) - np.log( m_max ) ) - 1.0
        
        return MF.Getks()[-1] / ( 1.0 - MF.Getalphas()[-1] ) * ( pow( 150.0, 1.0 - MF.Getalphas()[-1] ) - pow( m_max, 1.0 - MF.Getalphas()[-1] ) ) - 1.0
    
    
    def ComputeMF( self, Mini ):
        """computes the mass function for a given Mini
        Mini: the initial mass of the cluster [Msun]"""
        
        alpha = self.ComputeAlpha( Mini )
        
        logMini = np.log10( Mini )
        
        #PflammAltenburg (2007) approximation
        m_max = pow( 10.0, 2.56 * logMini * pow( pow( 3.82, 9.17 ) + pow( logMini, 9.17 ), -1.0/9.17 ) - 0.38 )
        
        #improve using Newton-Raphson
        epsilon = 1e-6      #allowed error
        h = 1.0             #step width
        
        for i in range( 100 ):
            MF = cMassFunction( Mini, [0.08,0.5,1.0,m_max], alpha )
            
            Delta = self.CheckUpperEnd( MF )
            
            if abs( Delta ) < epsilon:
                return MF
            
            Derivative = 0.5 * ( self.CheckUpperEnd( cMassFunction( Mini, [0.08,0.5,1.0,m_max + h], alpha ) ) - self.CheckUpperEnd( cMassFunction( Mini, [0.08,0.5,1.0,m_max - h], alpha ) )) / h
            
            m_max -= Delta / Derivative
        
        warnings.warn( "cIMFGenerator: ComputeMF: Newton-Raphson did not converge!" )
        
        return MF
    
    
    def ComputeMFFromToday( M, Age, Rperi, Rapo ):
        """computes the initial mass function based on todays data
        M: current mass of the GC
        Age: the current age of the GC [Gyr]
        Rperi: the pericentre of the GC [kpc]
        Rapo: the apocentre of the GC [kpc]"""
        
        e = ( Rapo - Rperi ) / ( Rapo + Rperi )
        
        #TODO
