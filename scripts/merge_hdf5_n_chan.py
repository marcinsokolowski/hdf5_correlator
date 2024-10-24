from __future__ import print_function
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
   usage+='\tpython merge_n_hdf5_files.py channel_cont_0_20190531_43409_0.hdf5 N_TILES[default = 3]\n'
   parser = OptionParser(usage=usage,version=1.00)
   # merged/merged_hdf5_list.txt
   parser.add_option('-o','--outdir','--merged_dir','--out_dir',dest="outdir",default="merged_channel0/", help="Name of output directory [default %]")
   parser.add_option('-f','--outlistfile','--output_hdf5_list',dest="outlistfile",default=None, help="Name of output file with list of merged HDF5 files [default %]")
   parser.add_option("-F", "--force", '--overwrite',action="store_true", dest="force", default=False, help="Force to overwrite otherwise already merged files are skipped [default %default]")
   parser.add_option('-c','--chan','--channel','--out_channel',dest="out_channel",default=0, help="Select one of the channels to be merged [default %]",type="int")
   parser.add_option('-n','--n_chan','--n_channel','--n_channels',dest="n_channels",default=2, help="Total number of channels in the dataset [default %]",type="int")

   (options, args) = parser.parse_args(sys.argv[idx:])
   
   if options.outlistfile is None :
      options.outlistfile = options.outdir + "/merged_hdf5_list.txt"

   return (options, args)



def print_usage() :
   print("merge_n_hdf5_files.py channel_cont_0_20190531_43409_0.hdf5 N_TILES[default = 3]")
   print("N_TILES specifies data from how many tiles is to be merged") 

# Script entry point
if __name__ == "__main__":

# file_list=["channel_cont_0_20190531_43409_0.hdf5","channel_cont_1_20190531_43409_0.hdf5","channel_cont_2_20190531_43409_0.hdf5"]
# file_list=["channel_cont_0_20190531_65848_0.hdf5","channel_cont_1_20190531_65848_0.hdf5","channel_cont_2_20190531_65848_0.hdf5"]
# file_list=["channel_cont_0_20190531_66632_0.hdf5","channel_cont_2_20190531_66632_0.hdf5","channel_cont_1_20190531_66632_0.hdf5"]
# file_list=["channel_cont_0_20190531_69143_0.hdf5","channel_cont_1_20190531_69143_0.hdf5","channel_cont_2_20190531_69143_0.hdf5"]

# 20190607/08 Galaxy transit files :
   file_list=[]
   file0 = None
   if len(sys.argv) > 1:
      file0 = sys.argv[1]
#      file1 = file0.replace("channel_cont_0_","channel_cont_1_")
#      file2 = file0.replace("channel_cont_0_","channel_cont_2_")
#      file_list = [ file0 , file1 , file2 ]
   if file0.find(".hdf5") >= 0 :
       file0 = file0.replace(".hdf5","")

   n_hdf5_files = 3 
   if len(sys.argv) > 2:
      n_hdf5_files = int( sys.argv[2] )
      
#   merged_file = "new.hdf5"
   merged_file = file0.replace("channel_cont_0_","channel_cont_") + ".hdf5"
   if len(sys.argv) > 3 and sys.argv[3] != "-" :   
       merged_file = sys.argv[3]

   (options, args) = parse_options( 3 )       
   outdir = options.outdir
#   outdir="merged/"
#   if len(sys.argv) > 4:    
#       outdir = sys.argv[4]

   merged_file_base = merged_file       
   merged_file = outdir + "/" + merged_file
   
   print("######################################################################")
   print("PARAMETERS :")
   print("######################################################################")
   print("file 0       = %s" % (file0))
   print("N hdf5 files = %d" % (n_hdf5_files))
   print("merged_file  = %s (base = %s)" % (merged_file,merged_file_base))
   print("outdir       = %s" % (outdir))
   print("Out list file = %s" % (options.outlistfile))
   print("######################################################################")

   # create output directory if does not exist already :
   mkdir_p( outdir )

   tile_file_list = []   
   if file0 is None :
      print("ERROR : name of the first file must be specified")
      print_usage()
      sys.exit(0)      
   else :
      for f in range(0,n_hdf5_files) :
          file_prefix = ( "channel_cont_%d_" % f )
          file = file0.replace("channel_cont_0_",file_prefix)
          file += ".hdf5"
          file_list.append( file )
          
          print("Opening file %s" % (file))         
          tile_file = h5py.File( file , mode="r" )
          tile_file_list.append( tile_file )
#          print "\tFile %s added to list : tile%d_file.shape = %d x %d" % (tile_file,f,tile_file['/chan_']['data'].shape[0],tile_file['/chan_']['data'].shape[1])
          print("File %s added to a list (%s)" % (file,tile_file))

   tile0_file = tile_file_list[0]
   print("tile0_file = %s" % (tile0_file))
   print("Merging %d hdf5 files : %s -> into one (%s)" % (len(tile_file_list),file_list,merged_file))   
   for f in range(0,len(tile_file_list)) :
      tile_file = tile_file_list[f]
#      print "\tfile = %s" % (tile_file)
      print("\t%s : tile%d_file.shape = %d x %d" % (tile_file,f,tile_file['/chan_']['data'].shape[0],tile_file['/chan_']['data'].shape[1]))
   
   # initialise merged file 
   datasets = getdatasets('/',tile0_file)
   new_data = h5py.File( merged_file , 'w' )
   # get the group-names from the lists of datasets
   groups = list(set([i[::-1].split('/',1)[1][::-1] for i in datasets]))
   groups = [i for i in groups if len(i)>0]

   # sort groups based on depth
   idx    = numpy.argsort(numpy.array([len(i.split('/')) for i in groups]))
   groups = [groups[i] for i in idx]

   # create all groups that contain dataset that will be copied
   for group in groups:
     new_data.create_group(group)

   # copy datasets
   for path in datasets:

      # - get group name
      group = path[::-1].split('/',1)[1][::-1]
      print("Copying group = %s" % (group))

      # - minimum group name
      if len(group) == 0: group = '/'

      # - copy data
      tile0_file.copy(path, new_data[group])


   # new_data['/chan_']['data'] = numpy.zeros( (10,10) )
   # new_data['/chan_']['data'].resize( (7408000,96) )
   # new_data['/chan_']['data'].clear()
   new_data['/chan_'].clear()
   # alldata = numpy.hstack( (tile0_file['/chan_']['data'],tile1_file['/chan_']['data'],tile2_file['/chan_']['data']) )
   alldata = None

   if len(tile_file_list) == 16 :
       tile0_file = tile_file_list[0]
       tile1_file = tile_file_list[1]
       tile2_file = tile_file_list[2]
       tile3_file = tile_file_list[3]
       tile4_file = tile_file_list[4]
       tile5_file = tile_file_list[5]
       tile6_file = tile_file_list[6]
       tile7_file = tile_file_list[7]

       tile8_file = tile_file_list[8]
       tile9_file = tile_file_list[9]
       tile10_file = tile_file_list[10]
       tile11_file = tile_file_list[11]
       tile12_file = tile_file_list[12]
       tile13_file = tile_file_list[13]
       tile14_file = tile_file_list[14]
       tile15_file = tile_file_list[15]
              
       # alldata = numpy.hstack( (tile0_file['/chan_']['data'] , tile1_file['/chan_']['data'] , tile2_file['/chan_']['data'] , tile3_file['/chan_']['data'] , tile4_file['/chan_']['data'] , tile5_file['/chan_']['data'] , tile6_file['/chan_']['data'] , tile7_file['/chan_']['data'] , tile8_file['/chan_']['data'] , tile9_file['/chan_']['data'] , tile10_file['/chan_']['data'] , tile11_file['/chan_']['data'] , tile12_file['/chan_']['data'] , tile13_file['/chan_']['data'] , tile14_file['/chan_']['data'] , tile15_file['/chan_']['data']   ) )
       alldata = numpy.hstack( (tile0_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile1_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile2_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile3_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile4_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile5_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile6_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile7_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile8_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile9_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile10_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile11_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile12_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile13_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile14_file['/chan_']['data'][:,options.out_channel::options.n_channels] , tile15_file['/chan_']['data'][:,options.out_channel::options.n_channels]   ) )
       print("INFO : 16 tiles is supported OK !")

   else :
       print("ERROR : only option to merge 16 tiles is implemented")
   
   shape0 = tile0_file['/chan_']['data'].shape[0]
   shape1 = tile0_file['/chan_']['data'].shape[1]
   data_type = tile0_file['/chan_']['data'][0,0].dtype
   
   merged_shape0 = shape0
   merged_shape1 = int(shape1*(1/options.n_channels)) * len(tile_file_list) # only saving out of of options.n_channels channels
   
   new_data['/chan_'].create_dataset("data", (merged_shape0,merged_shape1), dtype=data_type, data=alldata )
   print("INFO : saved merged files %s with shape (%d,%d)" % (merged_file,merged_shape0,merged_shape1))
   
   if options.outlistfile is not None :
      out_list_f = open( options.outlistfile , "a+" )
      line = "%s\n" % (merged_file_base)
      out_list_f.write( line )
      out_list_f.close()
      
      print("Saved merged filename (%s) to output list file %s" % (merged_file,options.outlistfile))


   # out_hdf_file['/chan_']['data'] = numpy.hstack( (tile0_file['/chan_']['data'],tile1_file['/chan_']['data'],tile2_file['/chan_']['data']) )
   # out_hdf_file['/chan_']['data'] = numpy.zeros( (tile0_file['/chan_']['data'].shape[0],tile0_file['/chan_']['data'].shape[1]) )
