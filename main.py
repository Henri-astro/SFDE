# general libs
import sys
import os

# own files
from src.Data import cData


#check that all input files are there
if 6 > len( sys.argv ):
    print( "Missing parameter.\nUsage: python " + sys.argv[0] + " <GC_property_file> <SN_file> <Ejecta_file> <Remnant_file> <output_folder>" )
    sys.exit()


#prepare the output folder
try:
    os.mkdir( sys.argv[4] )
except OSError as error:
    print( "The chosen output directory already exists! Please choose a new directory name." )


#read in data and create data struct
data = cData( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )
