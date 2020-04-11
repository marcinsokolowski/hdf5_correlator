from __future__ import print_function

import pdb
import h5py
import numpy
import sys
import os
import errno
from optparse import OptionParser,OptionGroup

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise


# function to return a list of paths to each dataset
def getdatasets(key,archive):

  if key[-1] != '/': key += '/'

  out = []

  for name in archive[key]:

    path = key + name

    if isinstance(archive[path], h5py.Dataset):
      out += [path]
    else:
       out += getdatasets(path,archive)

  return out

def parse_options(idx):
   usage="Usage: %prog [options]\n"
   usage+='\tpython get_hdf5_part.py IN_FILE OUTFILE\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-n','--n_int','--n_integrations',dest="out_n_int",default=30, help="Number of output integrations [default %]",type="int")
   parser.add_option('-o','--outdir','--merged_dir','--out_dir',dest="outdir",default="merged/", help="Name of output directory [default %]")
   parser.add_option('-f','--outlistfile','--output_hdf5_list',dest="outlistfile",default="merged/merged_hdf5_list.txt", help="Name of output file with list of merged HDF5 files [default %]")
   parser.add_option("-F", "--force", '--overwrite',action="store_true", dest="force", default=False, help="Force to overwrite otherwise already merged files are skipped [default %default]")

   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)



def print_usage() :
   print("get_hdf5_part.py IN_FILE OUTFILE")
   print("Creates a new hdf5 file which is subset of the original one") 

# Script entry point
if __name__ == "__main__":

   infile = None
   if len(sys.argv) > 1:
      infile = sys.argv[1]

   outfile = None
   if len(sys.argv) > 2:
      outfile = sys.argv[2]
      
   (options, args) = parse_options( 2 )       
   outdir = options.outdir

   print("######################################################################")
   print("PARAMETERS :")
   print("######################################################################")
   print("infile       = %s" % (infile))
   print("outfile      = %s" % (outfile))
   print("outdir       = %s" % (outdir))
   print("Out list file = %s" % (options.outlistfile))
   print("######################################################################")

   # create output directory if does not exist already :
   mkdir_p( outdir )

   print("Reading HDF5 file %s" % (infile))
   tile0_file = h5py.File( infile , mode="r" )
   datasets = getdatasets('/',tile0_file)
   
   new_data = h5py.File( outfile , 'w' )
   # get the group-names from the lists of datasets
   groups = list(set([i[::-1].split('/',1)[1][::-1] for i in datasets]))
   groups = [i for i in groups if len(i)>0]

   # sort groups based on depth
   idx    = numpy.argsort(numpy.array([len(i.split('/')) for i in groups]))
   groups = [groups[i] for i in idx]

   new_data.create_group('/correlation_matrix')
   new_data.create_group('/sample_timestamps')      
   tile0_file.copy('observation_info',new_data)
   tile0_file.copy('root',new_data)      
   new_data['correlation_matrix'].clear()
   new_data['sample_timestamps'].clear()
   n_integrations = tile0_file['correlation_matrix']['data'].shape[0]
   n_blocks = new_data['root'].attrs['n_blocks']
   new_data['root'].attrs['n_blocks'] = n_integrations
   print("Changing number of blocks in the %s file from %d -> %d" % (outfile,n_blocks,n_integrations))

#   f['correlation_matrix']['data'][378:478,:,:,:,]
   alldata = tile0_file['correlation_matrix']['data'][(n_integrations-options.out_n_int):(n_integrations),:,:,:,]
   times   = tile0_file['sample_timestamps']['data'][(n_integrations-options.out_n_int):(n_integrations)]
   print("Shapes : %s , %s" % (alldata.shape,times.shape))
   shape = tile0_file['correlation_matrix']['data'].shape
   new_shape = tile0_file['correlation_matrix']['data'].shape
   t_shape = tile0_file['sample_timestamps']['data'].shape
#   t_shape[0] = options.out_n_int
   data_type = alldata[0,0,0,0].dtype
   data_type_t = times[0].dtype
  
   new_data['correlation_matrix'].create_dataset("data", (options.out_n_int,shape[1],shape[2],shape[3]), dtype=data_type, data=alldata )
   new_data['sample_timestamps'].create_dataset("data", numpy.array([options.out_n_int]), dtype=data_type_t, data=times )
   print("INFO : saved merged files %s" % (outfile))
   
#   if options.outlistfile is not None :
#      out_list_f = open( options.outlistfile , "a+" )
#      line = "%s\n" % (merged_file_base)
#      out_list_f.write( line )
#      out_list_f.close()
#      
#      print("Saved merged filename (%s) to output list file %s" % (merged_file,options.outlistfile))

   new_data.close()

   # out_hdf_file['/chan_']['data'] = numpy.hstack( (tile0_file['/chan_']['data'],tile1_file['/chan_']['data'],tile2_file['/chan_']['data']) )
   # out_hdf_file['/chan_']['data'] = numpy.zeros( (tile0_file['/chan_']['data'].shape[0],tile0_file['/chan_']['data'].shape[1]) )
