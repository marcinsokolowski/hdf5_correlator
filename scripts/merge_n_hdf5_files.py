import h5py
import numpy
import sys
from optparse import OptionParser,OptionGroup


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
   parser.add_option('-o','--outdir','--merged_dir','--out_dir',dest="outdir",default="merged/", help="Name of output directory [default %]")
   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)



def print_usage() :
   print "merge_n_hdf5_files.py channel_cont_0_20190531_43409_0.hdf5 N_TILES[default = 3]"
   print "N_TILES specifies data from how many tiles is to be merged" 

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
       
   merged_file = outdir + "/" + merged_file

   print "######################################################################"
   print "PARAMETERS :"
   print "######################################################################"
   print "file 0       = %s" % (file0)
   print "N hdf5 files = %d" % (n_hdf5_files)
   print "merged_file  = %s" % (merged_file)
   print "outdir       = %s" % (outdir)
   print "######################################################################"

   tile_file_list = []   
   if file0 is None :
      print "ERROR : name of the first file must be specified"
      print_usage()
      sys.exit(0)      
   else :
      for f in range(0,n_hdf5_files) :
          file_prefix = ( "channel_cont_%d_" % f )
          file = file0.replace("channel_cont_0_",file_prefix)
          file += ".hdf5"
          file_list.append( file )
          
          print "Opening file %s" % (file)         
          tile_file = h5py.File( file , mode="r" )
          tile_file_list.append( tile_file )
#          print "\tFile %s added to list : tile%d_file.shape = %d x %d" % (tile_file,f,tile_file['/chan_']['data'].shape[0],tile_file['/chan_']['data'].shape[1])
          print "File %s added to a list (%s)" % (file,tile_file)

   tile0_file = tile_file_list[0]
   print "tile0_file = %s" % (tile0_file)
   print "Merging %d hdf5 files : %s -> into one (%s)" % (len(tile_file_list),file_list,merged_file)   
   for f in range(0,len(tile_file_list)) :
      tile_file = tile_file_list[f]
#      print "\tfile = %s" % (tile_file)
      print "\t%s : tile%d_file.shape = %d x %d" % (tile_file,f,tile_file['/chan_']['data'].shape[0],tile_file['/chan_']['data'].shape[1])
   
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
      print "Copying group = %s" % (group)

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
              
       alldata = numpy.hstack( (tile0_file['/chan_']['data'] , tile1_file['/chan_']['data'] , tile2_file['/chan_']['data'] , tile3_file['/chan_']['data'] , tile4_file['/chan_']['data'] , tile5_file['/chan_']['data'] , tile6_file['/chan_']['data'] , tile7_file['/chan_']['data'] , tile8_file['/chan_']['data'] , tile9_file['/chan_']['data'] , tile10_file['/chan_']['data'] , tile11_file['/chan_']['data'] , tile12_file['/chan_']['data'] , tile13_file['/chan_']['data'] , tile14_file['/chan_']['data'] , tile15_file['/chan_']['data']   ) )
       print "INFO : 16 tiles is supported OK !"

   elif len(tile_file_list) == 8 :
       tile0_file = tile_file_list[0]
       tile1_file = tile_file_list[1]
       tile2_file = tile_file_list[2]
       tile3_file = tile_file_list[3]
       tile4_file = tile_file_list[4]
       tile5_file = tile_file_list[5]
       tile6_file = tile_file_list[6]
       tile7_file = tile_file_list[7]
       
       alldata = numpy.hstack( (tile0_file['/chan_']['data'] , tile1_file['/chan_']['data'] , tile2_file['/chan_']['data'] , tile3_file['/chan_']['data'] , tile4_file['/chan_']['data'] , tile5_file['/chan_']['data'] , tile6_file['/chan_']['data'] , tile7_file['/chan_']['data']) )
       print "INFO : 8 tiles is supported OK !"
   elif len(tile_file_list) == 3 :
       tile0_file = tile_file_list[0]
       tile1_file = tile_file_list[1]
       tile2_file = tile_file_list[2]

       alldata = numpy.hstack( (tile0_file['/chan_']['data'] , tile1_file['/chan_']['data'] , tile2_file['/chan_']['data']) )
       print "INFO : 3 tiles is supported OK !"
   elif len(tile_file_list) == 4 :
       tile0_file = tile_file_list[0]
       tile1_file = tile_file_list[1]
       tile2_file = tile_file_list[2]
       tile3_file = tile_file_list[3]

       alldata = numpy.hstack( (tile0_file['/chan_']['data'] , tile1_file['/chan_']['data'] , tile2_file['/chan_']['data'], tile3_file['/chan_']['data']) )
       print "INFO : 4 tiles is supported OK !"
   else :
       print "ERROR : only option to merge 3, 4 or 8 tiles is implemented"
   
   shape0 = tile0_file['/chan_']['data'].shape[0]
   shape1 = tile0_file['/chan_']['data'].shape[1]
   data_type = tile0_file['/chan_']['data'][0,0].dtype
   
   merged_shape0 = shape0
   merged_shape1 = shape1 * len(tile_file_list)
   
   new_data['/chan_'].create_dataset("data", (merged_shape0,merged_shape1), dtype=data_type, data=alldata )
   print "INFO : saved merged files %s with shape (%d,%d)" % (merged_file,merged_shape0,merged_shape1)


   # out_hdf_file['/chan_']['data'] = numpy.hstack( (tile0_file['/chan_']['data'],tile1_file['/chan_']['data'],tile2_file['/chan_']['data']) )
   # out_hdf_file['/chan_']['data'] = numpy.zeros( (tile0_file['/chan_']['data'].shape[0],tile0_file['/chan_']['data'].shape[1]) )
