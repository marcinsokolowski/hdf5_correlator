from __future__ import print_function
import h5py 
import sys
from optparse import OptionParser,OptionGroup


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tGets basic information from the HDF5 file\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-s','--sample_time',dest="sample_time",default=1.08, help="Sample time in usec [default %default usec]",type="float")
   parser.add_option('-n','--n_chan','--n_channels','--fft_samples',dest="n_fft_samples",default=32, help="Number of FFT samples/channels [default %default]",type="int")
   parser.add_option('-i','--inttime',dest="inttime",default=0.2, help="Requested integration time to calculate number of chunks per uvfits file [default %default sec]",type="float")
   parser.add_option('-t','--type','--file_type','--filetype',dest="file_type",default="auto", help="Type of file if not specified -> automatic detection from file name [default %default]")
   
   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)

def hdf5_file_info( hdf5file , options, station_beam_file=False ) :
   sample_time   = options.sample_time # usec 
   n_fft_samples = options.n_fft_samples
   inttime       = options.inttime
            
   f = h5py.File( hdf5file , "r" )
   
   
   # is it correlated file or normal channelised file 
   is_corr_file = ( "correlation_matrix" in f.keys())
   inttime_hdf5 = f['root'].attrs['tsamp']
   inttime      = inttime_hdf5
   channel_id   = f['root'].attrs['channel_id']
   
   print("Information about HDF5 file %s" % (hdf5file))
   print("keys              = %s" % (f.keys()))
   print("Correlation file ? = %d" % (is_corr_file))
   print("Integration time   = %.8f [seconds]" % (inttime_hdf5))
   print("Frequency channel  = %d" % (channel_id))
   
   data_keyword="chan_"
   if is_corr_file :
      data_keyword = "correlation_matrix"
   if station_beam_file : 
      data_keyword = ("polarization_%d" % (0))
         
   print("f[data_keyword].keys() = %s" % (f[data_keyword].keys()))
   
   l = len( f[data_keyword]['data'].shape )
   print("len( f[%s]['data'].shape ) = %d" % (data_keyword,l))
   
   shape_str=""
   for a in range(0,l) :
      if a < (l-1) :
         shape_str += ("%d x " % f[data_keyword]['data'].shape[a])
      else :
         shape_str += ("%d" % f[data_keyword]['data'].shape[a])
    
   print("f['%s']['data'].shape = %s" % (data_keyword,shape_str))

   t0 = float( f['sample_timestamps']['data'][0] )
   print("%s : first timestamp = %.2f (unix time)" % (hdf5file,t0))      
   
   if l >= 2 :
      n_samples = int( f[data_keyword]['data'].shape[0] )
      n_inputs  = int( f[data_keyword]['data'].shape[1] )
      print("N_samples                 = %d" % (n_samples))
      print("N_inputs                  = %d" % (n_inputs))
      
      total_time_usec = n_samples * sample_time
      total_time_sec  = (total_time_usec/1000000.0)      
      print("Total time = %d [usec] = %.2f seconds" % (total_time_usec,total_time_sec))
      
      n_integrations = (n_samples / n_fft_samples)
      print("Number of  %d samples/channels spectra  = %d" % (n_fft_samples,n_integrations))
      
      n_integrations_per_uvfits = ( total_time_sec / inttime )
      print("n_integrations_per_uvfits (-n option of Lfile2uvfits_eda.sh ) = %d ( to get required inttime = %.4f [sec] from total_time_sec = %.4f [sec] )" % (n_integrations_per_uvfits,inttime,total_time_sec))
   

   timestamps = f['sample_timestamps']['data']
   last = timestamps.shape[0]
   t_start = timestamps[0]
   t_end   = timestamps[last-1]
   t_center = (t_start + t_end)/2.00
   print("Unixtime start = %.4f , end = %.4f -> center = %.4f" % (t_start,t_end,t_center))


def main() :   
   hdf5file="data.hdf5"
   if len(sys.argv) > 0:
      hdf5file = sys.argv[1]
      
   (options, args) = parse_options(1)
   
   is_station_beam = False
   if hdf5file.find("stationbeam_") >= 0 or hdf5file.find("station_beam") >= 0 or options.file_type.find("station") >= 0 :
      is_station_beam = True      
      print("is_station_beam set to True")
   print("Station beam = %s" % (is_station_beam))
 
   hdf5_file_info( hdf5file , options, station_beam_file=is_station_beam )     

#   if is_station_beam or options.file_type.find("station")>=0 :
#      station_file_info( hdf5file , options ) 
#   else "
      # if not is_station_beam or options.file_type.find("corr")>=0 :
#      corr_file_info( hdf5file , options )   
      

if __name__ == "__main__":
   main()
  
 