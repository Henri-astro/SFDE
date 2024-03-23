# general libs
import numpy as np

# own libs
from src.common import *
from src.IMFGenerator import cIMFGenerator
from src.RemnantCalculator import cRemnantCalculator
from src.StarExtractor import cStarExtractor


class cDataProcessor():
    """class responsible for the main data processing"""
    
    def __init__( self ):
        pass
    
    
    def ProcessData( self, data ):
        """performs the data processing on the data
        data: the data to work on"""
        
        self.ComputeIMF( data )
        self.ComputeIron( data )
        self.ComputeSNe( data )
    
    
    def ComputeIMF( self, data ):
        """computes all the IMF's and initial masses
        data: the data to work on
        This function adds the columns 'IMF' and 'Mini' to the data"""
        
        #extract the required columns
        FHs = data.AccessGCData( "Fe-H" )
        ZHs = [ Elem + 0.3 for Elem in FHs ]
        
        Ms = data.AccessGCData( "Mass" )
        Ages = data.AccessGCData( "Age" )
        Rapos = data.AccessGCData( "R_a" )
        Rperis = data.AccessGCData( "R_p" )
        
        #do the computation
        RemCalcs = [ cRemnantCalculator( data.GetRemnantData(), ZH ) for ZH in ZHs ]
        IMFGenerators = [ cIMFGenerator( RemCalc ) for RemCalc in RemCalcs ]
        IMFs = [ IMFGenerators[ nElem ].ComputeIMFFromToday( Ms[ nElem ], Ages[ nElem ], Rapos[ nElem ], Rperis[ nElem ] ) for nElem in range( len( Ms ) ) ]
        
        #write the results into the data
        data.AddGCData( "IMF", IMFs )
        data.AddGCData( "Mini", [ IMF.GetMtot() for IMF in IMFs ] )
    
    
    def ComputeIron( self, data ):
        """computes the amount of iron to be produced
        data: the data to work on
        adds the column 'ProducedIron' to the data"""
        
        #extract the required columns
        Minis = data.AccessGCData( "Mini" )
        SFEs = data.AccessGCData( "SFE" )
        FHs = data.AccessGCData( "Fe-H" )
        FeSpreads = data.AccessGCData( "FeSpread" )
        
        #do the calculations
        ProducedIron = [ pIronSun * ( pow( 10.0, FHs[nSC] + FeSpreads[nSC] ) - pow( 10.0, FHs[nSC] - FeSpreads[nSC] ) ) * Minis[nSC] * ( 1.0 / SFEs[nSC] - 1.0 ) for nSC in range( len( Minis ) ) ]
        
        #write results back into the data
        data.AddGCData( "ProducedIron", ProducedIron )
        
        
    def ComputeSNe( self, data ):
        """computes the number of Sne and the time star formation lasts
        data: the data to work on
        adds the columns 'NSN','mlast' and 'SFD'"""
        
        #extract the required columns
        IMFs = data.AccessGCData( "IMF" )
        ProducedIrons = data.AccessGCData( "ProducedIron" )
        FHs = data.AccessGCData( "Fe-H" )
        ZHs = [ Elem + 0.3 for Elem in FHs ]
        
        #do the calculations
        NSNe = []       #number of SNe exploded
        mlasts = []     #the mass of the last star to contribute to SF
        SFDs = []        #the star formation duration
        
        for nIMF in range( len( IMFs )):
            StarExtractor = cStarExtractor( IMFs[nIMF] )
            ProducedIron = ProducedIrons[nIMF]
            RemCalc = cRemnantCalculator( data.GetRemnantData(), ZHs[nIMF] )
            
            NSN = 0
            mlast = IMFs[nIMF].Getbounds()[-1]
            
            while ProducedIron > 0.0:
                nextMass = StarExtractor.GetNextMostMassiveStar()
                
                if nextMass < 8.0:
                    mlast = float( "NaN" )
                    NSN = float( "NaN" )
                    break
                
                if data.SNExplodes( nextMass ):
                    ProducedIron -= data.Ejecta( nextMass )
                    NSN += 1
                    mlast = nextMass
                    
            NSNe.append( NSN )
            mlasts.append( mlast )
            
            if np.isnan( mlast ):
                SFDs.append( float( "NaN" ) )
            else:
                SFDs.append( RemCalc.GetTimeFromMass( mlast ) )
            
        # add columns to data
        data.AddGCData( "NSN", NSNe )
        data.AddGCData( "mlast", mlasts )
        data.AddGCData( "SFD", SFDs )
