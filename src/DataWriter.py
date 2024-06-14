# general libs
import os
import sys
import matplotlib.pyplot as plt
import numpy as np


class cDataWriter():
    """a class to write all the output data into dedicated output files"""
    
    def __init__( self, OutputFolder ):
        """OutputFolder: the folder the data shall be written to"""
        
        self.__OutputFolder = OutputFolder
        
        try:
            os.mkdir( OutputFolder )
        except FileExistsError:
            print( "The choosen output folder already exists. Please choose another name." )
            sys.exit()
    
    
    def WriteAllData( self, data ):
        """generates all possible outputdata
        data: the data to generate the output from"""
        
        self.WriteGCData( data )
        self.PlotBarcodes( data )
    
    
    def WriteGCData( self, data ):
        """writes out all of the GCData using the data class
        data: the object in which all of the data is stored"""
        
        GCFile = open( self.__OutputFolder + "/GCData.dat", "w" )
        
        GCFile.write( data.AccessGCDataPrinteable() )
        
        GCFile.close()
    
    
    def PlotBarcodes( self, data ):
        """plots the barcodes for the GCs and puts them into dedicated output folders
        data: the data to work with"""
        
        #make dedicated output folder
        OutFolder = self.__OutputFolder + "/Barcodes"
        os.mkdir( OutFolder )
        
        #get required columns
        Names = data.AccessGCData( "Name" )
        IMFs = data.AccessGCData( "IMF" )
        mlasts = data.AccessGCData( "mlast" )
        SFDs = data.AccessGCData( "SFD" )
        
        # create the colour scheme
        nSamples = 550
        mmin = 8.0
        
        fig = plt.figure( num = 0, figsize = [3.5,1.5], dpi = 200)
        plt.rcParams.update({'font.size': 8})
        
        for nSC in range( len( Names )):
            if np.isnan( mlasts[nSC] ):
                continue
            
            mmax = IMFs[ nSC ].Getbounds()[-1]
            
            logDist = ( np.log10( mmax ) - np.log10( mmin )) / nSamples
            
            Masses = [ pow( 10.0, np.log10( mmin ) + i * logDist ) for i in range( nSamples + 1 ) ]
            
            Colours = []
            
            for mass in Masses:
                if not data.SNExplodes( mass ):
                    Colours.append( "black" )
                elif mass < mlasts[nSC]:
                    Colours.append( "grey" )
                else:
                    Colours.append( "red" )
                    
            plt.bar( Masses, 1, color = Colours, lw = 10 )
            plt.xscale( "log" )
            
            plt.xlabel( "$m [M_\\odot]$" )
            plt.text( mlasts[nSC], 1.05, " $t_{SF} = " + "{:.1f}".format( SFDs[nSC] * 1000.0 ) + "$ Myr ", ha = "left" if mlasts[nSC] < 80 else "right" )
            plt.arrow( mlasts[nSC] - 0.1, 1.15, 0.0, -0.11, color = "black", lw = 1, clip_on=False, head_width = 0.03 * mlasts[nSC], head_length = 0.03, head_starts_at_zero = False )
            
            plt.gca().axes.get_yaxis().set_visible(False)
            
            plt.xlim( xmin = mmin, xmax = mmax )
            plt.ylim( ymin = 0, ymax = 1 )
            
            plt.subplots_adjust(left=0.02, right=0.98, top=0.88, bottom=0.3)
            
            plt.savefig( OutFolder + "/" + Names[nSC] + ".png" )
            plt.clf()
