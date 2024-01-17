# general libs
import sys

# own files
from src.Data import cData
from src.DataProcessor import cDataProcessor
from src.DataWriter import cDataWriter


#check that all input files are there
if 6 > len( sys.argv ):
    print( "Missing parameter.\nUsage: python " + sys.argv[0] + " <GC_property_file> <SN_file> <Ejecta_file> <Remnant_file> <output_folder>" )
    sys.exit()

#read in data and create data struct
data = cData( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )

#process the data
Processor = cDataProcessor()
Processor.ProcessData( data )

#write the output data
DataWriter = cDataWriter( sys.argv[5] )
DataWriter.WriteAllData( data )
