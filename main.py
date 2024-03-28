# general libs
import sys

# own files
from src.Data import cData
from src.DataProcessor import cDataProcessor
from src.DataWriter import cDataWriter


#check that all input files are there and read in data and create data struct
if 4 == len( sys.argv ):
    data = cData( sys.argv[1], sys.argv[2], sys.argv[2], sys.argv[2] )
    DataWriter = cDataWriter( sys.argv[3] )
elif 6 <= len( sys.argv ):
    data = cData( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] )
    DataWriter = cDataWriter( sys.argv[5] )
else:
    print( "Missing parameter.\nUsage: python " + sys.argv[0] + " <GC_property_file> <SN_file> <Ejecta_file> <Remnant_file> <output_folder>\nor: python " + sys.argv[0] + " <GC_property_file> <Remnant_file> <output_folder>\nMore information in the documentation." )
    sys.exit()

#process the data
Processor = cDataProcessor()
Processor.ProcessData( data )

#write the output data
DataWriter.WriteAllData( data )
