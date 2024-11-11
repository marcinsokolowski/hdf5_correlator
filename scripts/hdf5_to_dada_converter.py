from __future__ import print_function
# see my notes in /home/msok/Desktop/AAVS1/logbook/20190520_feeding_newmwacorr_with_tpm_data_stream.odt

# debug :
# import pdb

import sys
import os
import math
import numpy
import time, datetime
from optparse import OptionParser,OptionGroup
import h5py
import logging

g_file_number = 0

SKA_low_channel_separation = (400.00/512.00)
SKA_low_oversampling_ratio = (32.00/27.00)
SKA_sampling_time_usec     = 1.08 # micro-seconds


def get_dada_filename( uxtime = 0, tile_id=0, station_id=0 ) :
    utc_string = datetime.datetime.utcfromtimestamp( uxtime ).strftime( "%Y%m%d_%H%M%S" )
    uxtime_usec = ( uxtime - int(uxtime) )*1000000.00

    if tile_id is not None :    
        dada_filename = ( "psrdada_%s_%06d_station%04d_tile%04d.dada" % (utc_string,int(uxtime_usec),station_id,tile_id) )
    else :
        dada_filename = ( "psrdada_%s_%06d_station%04d.dada" % (utc_string,int(uxtime_usec),station_id) )
    
    return dada_filename


# out_header += ("%d\n") % ()
def generate_dada_header( start_uxtime=0, 
                          obsid = 0, 
                          nbit=32, 
                          npol=2, 
                          ndim=2, # 1-real, 2-complex
                          ntimesamples=10240000, 
                          ninputs=256,
                          ninputs_xgpu=256,
                          inttime_msec=1000,
                          proj_id="G0008",
                          exptime_sec=8,
                          file_size  = 5269094400,
                          file_number=0,
                          n_fine_channels=40,
                          bandwidth_hz=((400.00/512.00)*(32.00/27.00))*1e6, # MWA : 1280000,
                          frequency_mhz=(204*(400.00/512.00)),
                          telescope="LFAASP",
                          source="B0950+08",
                          ra=0.00,
                          dec=0.00
                        ) :
   # 
   if start_uxtime <= 0 :
       start_uxtime = time.time()

   # default values from Ian's example file ( /home/msok/askap/craft/data/J0953/incoherent_beam_1kHz_1ms.dada )
   header_size = 4096
   obs_offset  = 0     
   inttime_usec = inttime_msec*1000.00 
   

   # see  : http://dspsr.sourceforge.net/manuals/dspsr/dada.shtml
   out_header = ("HDR_VERSION 1.0\n")
   out_header += ("HDR_SIZE %d\n") % (header_size)
   out_header += ("BW %.4f\n") % (bandwidth_hz/1e6)
   out_header += ("FREQ %.4f\n") % (frequency_mhz)
   out_header += ("TELESCOPE %s\n") % (telescope)
   out_header += ("RECEIVER %s\n") % (telescope)
   out_header += ("INSTRUMENT %s\n") % (telescope)
   out_header += ("SOURCE %s\n") % (source)
   out_header += ("MODE PSR\n")
   out_header += ("NBIT %d\n") % (nbit)
   out_header += ("NPOL %d\n") % (npol)
   out_header += ("NCHAN %d\n") % (n_fine_channels)
   out_header += ("NDIM %d\n") % (ndim)
   out_header += ("OBS_OFFSET %d\n") % (obs_offset)
   out_header += ("TSAMP %.4f\n") % (inttime_usec)
   # utc_string = time.strftime("%Y-%m-%d-%H:%M:%S", str(start_uxtime))
   utc_string = datetime.datetime.utcfromtimestamp( start_uxtime ).strftime( "%Y-%m-%d-%H:%M:%S" )
   out_header += ("UTC_START %s\n" % (utc_string))
#   out_header += ("RA  %s\n" % ("08:35:20.61149"))
#   out_header += ("DEC %s\n" % ("-45:10:34.8751"))
      

   # non-crucial :         
   out_header += "POPULATED 1\n"
   out_header += ("OBS_ID %d\n") % (obsid)
   out_header += ("SUBOBS_ID %d\n") % (obsid)
   out_header += "COMMAND CAPTURE\n"
   
   # UTC_START 2018-10-11-05:26:14
   out_header += ("NTIMESAMPLES %d\n") % (ntimesamples)
   out_header += ("NINPUTS %d\n") % (ninputs)
   out_header += ("NINPUTS_XGPU %d\n") % (ninputs_xgpu)
   out_header += ("METADATA_BEAMS %d\n") % (2)
   out_header += ("APPLY_PATH_WEIGHTS %d\n") % (1)
   out_header += ("APPLY_PATH_DELAYS %d\n") % (2)
   out_header += ("INT_TIME_MSEC %d\n") % (inttime_msec)
   out_header += ("FSCRUNCH_FACTOR %d\n") % ((ntimesamples/2)/n_fine_channels)
   out_header += ("TRANSFER_SIZE %d\n") % (81920000)
   out_header += ("PROJ_ID %s\n") % (proj_id)
   out_header += ("EXPOSURE_SECS %d\n") % (exptime_sec)
   out_header += ("COARSE_CHANNEL %d\n") % (76)
   out_header += ("CORR_COARSE_CHANNEL %d\n") % (2)
   out_header += ("SECS_PER_SUBOBS %d\n") % (8)
   out_header += ("UNIXTIME %d\n") % (int(start_uxtime))
   uxtime_msec = ( start_uxtime - int(start_uxtime) )*1000.00
   out_header += ("UNIXTIME_MSEC %d\n") % (int(uxtime_msec))
   out_header += ("FINE_CHAN_WIDTH_HZ %d\n") % (bandwidth_hz/n_fine_channels)
   out_header += ("NFINE_CHAN %d\n") % (n_fine_channels)
   out_header += ("BANDWIDTH_HZ %d\n") % (bandwidth_hz)
   out_header += ("SAMPLE_RATE %d\n") % (bandwidth_hz)
   out_header += ("MC_IP 0.0.0.0\n")
   out_header += ("MC_PORT 0\n")
   out_header += ("MC_SRC_IP 0.0.0.0\n")
   out_header += ("FILE_SIZE %d\n") % (file_size)
   out_header += ("FILE_NUMBER %d\n") % (file_number)
#   out_header += (" %d\n") % ()
#   out_header += (" %d\n") % ()

   missing = ( header_size - len(out_header) )
   for i in range(0,missing) :
       out_header += '\0'


   return out_header
   
   
def save_psrdada_file( data_filename, data=None,
                       start_uxtime=0, 
                       obsid = 0, 
                       nbit=32, 
                       npol=2, 
                       ndim=2,
                       ntimesamples=10240000, 
                       ninputs=256,
                       ninputs_xgpu=256,
                       inttime_msec=1000,
                       proj_id="G0008",
                       exptime_sec=8,
                       file_size  = 5269094400,
                       file_number=0,
                       n_fine_channels=40,
                       bandwidth_hz=((400.00/512.00)*(32.00/27.00))*1e6, # MWA : 1280000,
                     ) :
    header = generate_dada_header( start_uxtime=start_uxtime, obsid=obsid, nbit=nbit, npol=npol, ndim=ndim, ntimesamples=ntimesamples, ninputs=ninputs, ninputs_xgpu=ninputs_xgpu, 
                                   inttime_msec=inttime_msec, proj_id=proj_id, exptime_sec=exptime_sec, file_size=file_size, file_number=file_number , n_fine_channels=n_fine_channels, bandwidth_hz=bandwidth_hz
                                 )
    
    dada_file = open(  data_filename, "wb" )
    
    # write header :
    dada_file.write( header )
    
    # write data :
    # dada_file.write( data )    
    # numpy.save( dada_file, data, allow_pickle=False )
#    data.astype('float32').tofile( dada_file )
#    data.astype('complex64').tofile( dada_file )
#    data.astype('complex16').tofile( dada_file )
#    dtype=[('real', 'i1'), ('imag', 'i1')]
    dtype=data.dtype
    data.astype(dtype).tofile( dada_file )
    
    dada_file.close()


    return True


def read_and_convert_dat2psrdada_file( data_filename, dada_filename, 
                       start_uxtime=0, 
                       obsid = 0, 
                       nbit=32, 
                       npol=2, 
                       ndim=2,
                       ntimesamples=10240000, 
                       ninputs=256,
                       ninputs_xgpu=256,
                       inttime_msec=1000,
                       proj_id="G0008",
                       exptime_sec=8,
                       file_size  = 5269094400,
                       file_number=0,
                       n_fine_channels=40,
                       bandwidth_hz=((400.00/512.00)*(32.00/27.00))*1e6, # MWA : 1280000,
                       buffer_size=1000000,
                       complex_mult=None,
                       frequency_mhz=(204*(400.00/512.00)),
                       conjugate=False
                     ) :
    header = generate_dada_header( start_uxtime=start_uxtime, obsid=obsid, nbit=nbit, npol=npol, ndim=ndim, ntimesamples=ntimesamples, ninputs=ninputs, ninputs_xgpu=ninputs_xgpu, 
                                   inttime_msec=inttime_msec, proj_id=proj_id, exptime_sec=exptime_sec, file_size=file_size, file_number=file_number , n_fine_channels=n_fine_channels, bandwidth_hz=bandwidth_hz,
                                   frequency_mhz=frequency_mhz
                                 )
    
    dada_file = open(  dada_filename, "wb" )
    # write header :
    dada_file.write( header )

    complex_8t = numpy.dtype([('real', numpy.int8), ('imag', numpy.int8)])

    file_pos = 0
    f = open( data_filename, 'rb')
    
    conj_f = None
    if conjugate :
       # for debugging and test purposes :
       conj_f = open( "conjugate.dat" , "wb" )
    
    read_ok = True
    while read_ok :
       f.seek( file_pos )
       data = numpy.fromfile(f, dtype=numpy.int8, count=buffer_size)
       
       if len(data) > 0 :
          file_pos += len(data)
#          if mult is not None :
#             data = data * complex_mult

          if conjugate :          
             for i in range(1,len(data),2) :
#            for i in range(0,len(data),2) :
                data[i] = -data[i]

          # swap RE/IM
#          for i in range(1,len(data)) :
#             if ( i % 2 ) == 1 : 
#                tmp = data[i-1]
#                data[i-1] = data[i]
#                data[i] = tmp
                          
          data.astype(numpy.int8).tofile( dada_file )
          
          if conj_f is not None :
             data.astype(numpy.int8).tofile( conj_f )
          
          read_ok = True
          print("DEBUG : saved %d bytes in total" % (file_pos))
       else :
          print("DEBUG : no data after reading %d bytes" % (file_pos))
          read_ok = False
    

#    read_ok = True
#    while read_ok :
#       data = f.read( buffer_size )
#       
#       if len(data) > 0 :
#          file_pos += len(data)
##          data.astype(numpy.int8).tofile( dada_file )
#          read_ok = True
#          print("DEBUG : saved %d bytes in total" % (file_pos))
#       else :
#          print("DEBUG : no data after reading %d bytes" % (file_pos))
#          read_ok = False
       
    
    f.close()
    dada_file.close()
    
    if conj_f is not None :
       conj_f.close()

    return True


def create_test_dada_file() :
   test_dada_file = "test.dada"      
   
   # create test data :
   nbit=32 
   npol=2 
   ntimesamples=10240000 
   ninputs=256
   ninputs_xgpu=256
   inttime_msec=1000
   proj_id="G0008"
   exptime_sec=8
   file_size  = 5269094400
   file_number=0

   n_channels  = 1280
   n_timesteps = 1000 # same as in file /home/msok/askap/craft/data/J0953/incoherent_beam_1kHz_1ms.dada   
   n_blocks    = 4 
   n_reim      = 2    # RE/IM
   n_seconds   = 8
   

   # generate template spectrum to fill the test data :
   template_spectrum = numpy.zeros( n_channels * n_reim, dtype=numpy.float32 )
   for ch in range(0,n_channels) :
       idx = ch*n_reim 
       value = n_reim*ch

       template_spectrum[idx]   = value # fill only RE part with ch value
       template_spectrum[idx+1] = value + 1
       
       print("%d %.4f" % (idx,value))
       print("%d %.4f" % (idx+1,value+1))

   # sizes in bytes :
   one_spectrum_size = n_channels*n_reim
   one_block_size    = n_timesteps*one_spectrum_size 
   one_second_size   = one_block_size*n_blocks
      
   data = numpy.zeros( one_second_size*n_seconds, dtype=numpy.float32 )   
   print("Allocated array of %d bytes" % (one_second_size*n_seconds*4))

   for s in range(0,n_seconds) :
       skip_seconds_size = s*(one_second_size)
       for b in range(0,n_blocks) :
          skip_blocks_size = b*one_block_size + s*(one_second_size)
          for t in range(0,n_timesteps) :
              for ch in range(0,n_channels) :
                  # idx = n_reim*ch + t*one_spectrum_size + b*one_block_size + s*(one_second_size)
                  idx = n_reim*ch + t*one_spectrum_size + skip_blocks_size
              
                  data[idx]   = template_spectrum[n_reim*ch]
                  data[idx+1] = template_spectrum[n_reim*ch + 1]
                  
                  if t==0 and b==0 and s==0 and True :
                      print("DEBUG1 : %d %.4f" % (2*ch,data[idx]))
                      print("DEBUG1 : %d %.4f" % (2*ch+1,data[idx+1]))

   print("Filled data array")
   for ch in range(0,n_channels) :
      idx = ch*n_reim 
      print("DEBUG2 : %d %.4f" % (idx,data[idx]))
      print("DEBUG2 : %d %.4f" % (idx+1,data[idx+1]))
             
   save_psrdada_file( test_dada_file, data=data )

def convert_hdf5data_to_dada( data_ptr, timestamp=0, dada_filename=None, debug_level=0, do_reshape=False, 
                              n_blocks   = 200, # for 7408000 samples, 25 for 926000 samples  
                                                # Ian's number 1209600 samples x 32 = 38707200 total samples
                              n_fine_channels = 40,
                              bandwidth_hz    = 926000, # really it is rather ~ ((800.00/2.00)/512.00)*(32.00/27.00) # ( 400 MHz / 512 )*(32/27) ~= 926000
                              n_antennas      = 16,
                              n_pols          =  2
                            ): 
   global g_file_number
   
   if timestamp <= 0 :
       timestamp = time.time()   
   
   if dada_filename is None :
      dada_filename = get_dada_filename( timestamp )

   n_total_samples =  data_ptr.shape[0]
   if len(data_ptr.shape) > 1 :
       n_total_samples = 1
       for i in range(0,len(data_ptr.shape)) :
           n_total_samples = n_total_samples * data_ptr.shape[i]

   print("convert_hdf5data_to_dada : n_total_samples = %d , n_blocks = %d, n_fine_channels = %d, bandwidth_hz = %d, n_antennas = %d" % (n_total_samples,n_blocks,n_fine_channels,bandwidth_hz,n_antennas))
   
   # n_samp     = 11575
   n_inputs = (n_antennas*n_pols)
   n_samples_per_input = n_total_samples / n_inputs
   
   if do_reshape :       
       data_ptr = numpy.reshape(data_ptr, (n_samples_per_input, n_inputs))
       logging.info('\tpost re-shape = %d x %d' % (data_ptr.shape[0],data_ptr.shape[1]))
   else :
       logging.info('\tre-shape is not required')

   # creating n_blocks buffers :
   n_samples_per_block = n_samples_per_input / n_blocks
   n_tot_samples_per_block = n_total_samples / n_blocks
   logging.info('\tn_samples_per_input / n_blocks = %d / %d = %d , n_tot_samples_per_block = %d / %d = %d' % (n_samples_per_input,n_blocks,n_samples_per_block,n_total_samples,n_blocks,n_tot_samples_per_block))

   data_buffer = numpy.zeros((n_blocks+1,n_tot_samples_per_block),dtype=data_ptr.dtype)# was numpy.complex64
   for b in range(0,n_blocks) :
       for inp in range(0,n_inputs) :
           # print "DEBUG : first sample in block %d / input %d = %s" % ( b, inp, data_ptr[0,inp] )
           # b+1 - do leave the 0-block empty (ZEROS/METADATA) :
           data_buffer[b+1,inp*n_samples_per_block:(inp+1)*n_samples_per_block] = data_ptr[b*n_samples_per_block:(b+1)*n_samples_per_block,inp].copy()
           if debug_level > 0 :
               print("DEBUG : first sample in block %d / input %d = %s vs. %s" % ( b, inp, data_buffer[b+1,inp*n_samples_per_block], data_ptr[b*n_samples_per_block,inp] ))

   gpstime = timestamp - 315964783
   save_psrdada_file( dada_filename, 
                      data_buffer, 
                      start_uxtime = timestamp, 
                      obsid = gpstime, 
                      nbit  = data_buffer[0,0]['real'].nbytes*8, # nbits 
                      npol  = 2,
                      ntimesamples = n_tot_samples_per_block / n_inputs,
                      ninputs = n_inputs,
                      ninputs_xgpu = n_inputs,
                      inttime_msec = 0,
                      proj_id      = "EDA1TPM20",
                      exptime_sec  = ((n_samples_per_input*1.08)/1000000.00),                           
                      file_size    = n_tot_samples_per_block*4*2, # 4 = sizeof(float)*2 (real/imag)
                      file_number  = g_file_number,
                      n_fine_channels = n_fine_channels,
                      bandwidth_hz = bandwidth_hz
                    )

   g_file_number = g_file_number + 1
   
   return g_file_number
   
def convert_stationbeam_hdf5data_to_dada( data_ptr, timestamp=0, dada_filename=None, debug_level=0, do_reshape=False, 
                                          n_blocks   = 200, # for 7408000 samples, 25 for 926000 samples  
                                                # Ian's number 1209600 samples x 32 = 38707200 total samples
                                          n_fine_channels = 40,
                                          bandwidth_hz    = 926000, # really it is rather ~ ((800.00/2.00)/512.00)*(32.00/27.00) # ( 400 MHz / 512 )*(32/27) ~= 926000
                                          n_antennas      = 16,
                                          n_pols          =  2
                                        ): 
   global g_file_number
   
   if timestamp <= 0 :
       timestamp = time.time()   
   
   if dada_filename is None :
      dada_filename = get_dada_filename( timestamp )

   n_total_samples =  data_ptr.shape[0]
   if len(data_ptr.shape) > 1 :
       n_total_samples = 1
       for i in range(0,len(data_ptr.shape)) :
           n_total_samples = n_total_samples * data_ptr.shape[i]

   print("convert_stationbeam_hdf5data_to_dada : n_total_samples = %d , n_blocks = %d, n_fine_channels = %d, bandwidth_hz = %d, n_antennas = %d" % (n_total_samples,n_blocks,n_fine_channels,bandwidth_hz,n_antennas))
   
   gpstime = timestamp - 315964783
   save_psrdada_file( dada_filename, 
                      data_ptr, 
                      start_uxtime = timestamp, 
                      obsid = gpstime, 
                      nbit  = data_ptr[0].nbytes*8, # nbits 
                      npol  = 1,
                      ntimesamples = n_total_samples,
                      ninputs = 1,
                      ninputs_xgpu = 1,
                      inttime_msec = 0,
                      proj_id      = "EDA2_STATION_BEAM",
                      exptime_sec  = ((n_total_samples*1.08)/1000000.00),                           
                      file_size    = n_total_samples*data_ptr[0].nbytes, # 4 = sizeof(float)*2 (real/imag)
                      file_number  = g_file_number,
                      n_fine_channels = 1,
                      bandwidth_hz = bandwidth_hz
                    )

   g_file_number = g_file_number + 1
   
   return g_file_number
                                        
   

def convert_hdf5_to_dada( hdf5file, dadafile, timestamp=0, debug_level=0, n_antennas=16 ) :
   print("Converting HDF5 file %s to dada file %s" % (hdf5file,dadafile))
   f = h5py.File( hdf5file )
   data_ptr = f['/chan_']['data']
   print("Data shape = %d x %d" % (data_ptr.shape[0],data_ptr.shape[1]))

   convert_hdf5data_to_dada( data_ptr, timestamp=timestamp, dada_filename=dadafile, n_antennas=n_antennas )

def convert_stationbeam_hdf5_to_dada( hdf5file, dadafile, timestamp=0, debug_level=0, skip_n_samples=0, polarisation=0 ) :
   print("Converting HDF5 file %s to dada file %s" % (hdf5file,dadafile))
   f = h5py.File( hdf5file )
   
   data_keyword=( 'polarization_%d' % (polarisation) )
   
   data_ptr = None
   if skip_n_samples > 0 :
      data_ptr = f[data_keyword]['data'][skip_n_samples:,0]
   else :
      data_ptr = f[data_keyword]['data'][:,0]
      
   print("Data shape = %d , first and last values = %.4f , %.4f" % (data_ptr.shape[0],data_ptr[0],data_ptr[data_ptr.shape[0]-1]))
   convert_stationbeam_hdf5data_to_dada( data_ptr, timestamp=timestamp, dada_filename=dadafile, n_antennas=1 )


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tConvert hdf5 file with continous channel data to PSRDADA format for the new correlator\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-d','--debug','--debug_level',dest="debug_level",default=0, help="Debug level [default %]",type="int")
   parser.add_option('-a','--nants','--n_antennas','--n_ants',dest="n_antennas",default=256, help="Number of antennas [default %]",type="int")
   parser.add_option('-P','--npol','--n_pol',dest="npol",default=2, help="Number of polarisations [default %]",type="int")
   parser.add_option('-p','--pol','--polarisation',dest="polarisation",default=0, help="Polarisation 0 is X and 1 is Y [default %]",type="int")
   parser.add_option('--skip_n_samples','--skip','--skip_n',dest="skip_n_samples",default=0,help="Skip N samples [default %]",type="int")
   parser.add_option('-c','--hdf52dada',action="store_true",dest="hdf52dada",default=False, help="Convert hdf5 file to .dada file [default %]")
   parser.add_option('-s','--station_beam',action="store_true",dest="station_beam",default=False, help="Treat HDF5 file as station beam file [default %]")
   parser.add_option('-u','--uxtime','--unix_time','--unixtime','--start_uxtime','--start_unix_time',dest="start_unix_time",default=0, help="Unixtime of the first sample [default %]",type="float")
   parser.add_option('--hdr','--dadahdr','--psrdadahdr',dest="generete_dada_header",action="store_true",default=False, help="Generate PSRDADA header [default %]")
   parser.add_option('-o','--outfile','--out_file','--output_file',dest="output_file",default=None, help="Output file name to save DADA header [default %]")
   parser.add_option('--dat2dada',action="store_true",dest="dat2dada",default=False, help="Convert dat file to .dada file [default %]")
   parser.add_option('--dat2dada_old',action="store_true",dest="dat2dada_old",default=False, help="OBSOLATE VERSION : Convert dat file to .dada file [default %]")
   parser.add_option('--byte2float',action="store_true",dest="byte2float",default=False, help="Byte to float [default %]")
   parser.add_option('--ra','--right_ascension','--ra_deg',dest="ra_deg",default=0.00, help="RA [deg] [default %]",type="float")
   parser.add_option('--dec','--declination','--dec_deg',dest="dec_deg",default=0.00, help="DEC [deg] [default %]",type="float")
   parser.add_option('--ch','--freq_ch',dest="freq_ch",default=204, help="Frequency channel [default %]",type="int")
   parser.add_option('--source','--object',dest="source",default="B0950+08", help="Observed source [default %]")
   parser.add_option('--multiplier','--complex_multiplier','--mult_complex',dest="complex_multiplier",default=None, help="Complex multiplier [default %]")
   
   parser.add_option('--conjugate',dest="conjugate",action="store_true",default=False, help="Conjugate [default %default]")

   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)

    
if __name__ == '__main__':
    hdf5file="channel_cont_0_20190524_06638_0.hdf5"
    if len(sys.argv) > 1:
       hdf5file = sys.argv[1]   


    (options, args) = parse_options(1) 

    if options.hdf52dada :
        dadafile=hdf5file.replace('.hdf5', '.dada')

        if options.station_beam :
           convert_stationbeam_hdf5_to_dada( hdf5file, dadafile, skip_n_samples=options.skip_n_samples, 
                                             debug_level=options.debug_level, timestamp=options.start_unix_time,
                                             polarisation=options.polarisation )
        else :    
           convert_hdf5_to_dada( hdf5file, dadafile,
                                 n_antennas = options.n_antennas,
                                 debug_level=options.debug_level )
    elif options.dat2dada_old :
        dadafile=hdf5file.replace('.dat', '.dada')
        if options.output_file is not None :
           dadafile = options.output_file
           
        print("Converting file %s to %s" % (hdf5file,dadafile))
        
        # Custom numpy type for creating complex signed 8-bit data
        complex_8t = numpy.dtype([('real', numpy.int8), ('imag', numpy.int8)])

        # Convert to complex
        data = numpy.fromfile(hdf5file, dtype=complex_8t) # count=all ?
        inttime_msec = SKA_sampling_time_usec / 1000.00 # 1.08 usec -> miliseconds
        
        npol = 2 # X and Y 
        ndim = 2 # re/im
        file_size = os.stat( hdf5file ).st_size # was 1073741824
        data_complex = data['real'].astype(numpy.float32) + 1j * data['imag'].astype(numpy.float32)
        n_timestamps = data_complex.shape[0] / options.npol
        bandwidth_hz  = SKA_low_channel_separation*SKA_low_oversampling_ratio*1e6
        inttime_msec = (SKA_sampling_time_usec / 1000.00)
                            
#        header = generate_dada_header( start_uxtime=options.start_unix_time, obsid=0, nbit=16, npol=2, ntimesamples=	
        data_file = save_psrdada_file( dadafile, data=data_complex, start_uxtime=options.start_unix_time, obsid=0, nbit=8, npol=npol, ndim=ndim,	
                                       ntimesamples=n_timestamps, ninputs=2, ninputs_xgpu=2, inttime_msec=inttime_msec, proj_id = "LFAASP", 
                                       exptime_sec = (n_timestamps*inttime_msec/1000.00), file_size=data.shape[0], n_fine_channels=1, bandwidth_hz=bandwidth_hz
                                     )

    elif options.dat2dada :
        dadafile=hdf5file.replace('.dat', '.dada')
        if options.output_file is not None :
           dadafile = options.output_file
           
        print("Converting file %s to %s" % (hdf5file,dadafile))
        
        inttime_msec = SKA_sampling_time_usec / 1000.00 # 1.08 usec -> miliseconds
        
        npol = 2 # X and Y 
        ndim = 2 # re/im
        file_size = os.stat( hdf5file ).st_size # was 1073741824
        n_timestamps = file_size / options.npol
        bandwidth_hz  = SKA_low_channel_separation*SKA_low_oversampling_ratio*1e6
        inttime_msec = (SKA_sampling_time_usec / 1000.00)
        frequency_mhz = (options.freq_ch*SKA_low_channel_separation)
                            
#        header = generate_dada_header( start_uxtime=options.start_unix_time, obsid=0, nbit=16, npol=2, ntimesamples=	
        data_file = read_and_convert_dat2psrdada_file( hdf5file , dadafile, start_uxtime=options.start_unix_time, obsid=0, nbit=8, npol=npol, ndim=ndim,	
                                       ntimesamples=n_timestamps, ninputs=2, ninputs_xgpu=2, inttime_msec=inttime_msec, proj_id = "LFAASP", 
                                       exptime_sec = (n_timestamps*inttime_msec/1000.00), file_size=file_size, n_fine_channels=1, bandwidth_hz=bandwidth_hz,
                                       frequency_mhz=frequency_mhz, conjugate=options.conjugate
                                     )


    elif options.generete_dada_header : 
        hdrfile=hdf5file.replace('.dat', '.hdr')                        
        if options.output_file is not None : 
           hdrfile = options.output_file

        npol = 2 # X and Y 
        ndim = 2 # re/im
        file_size = os.stat( hdf5file ).st_size # was 1073741824
        ntimesamples = file_size / (npol*ndim) # was 268435456
        inttime_msec=(SKA_sampling_time_usec / 1000.00)
        frequency_mhz = (options.freq_ch*SKA_low_channel_separation)
        bandwidth_hz  = SKA_low_channel_separation*SKA_low_oversampling_ratio*1e6
        exptime_sec = (ntimesamples*inttime_msec/1000.00)
        
        header = generate_dada_header( start_uxtime=options.start_unix_time, obsid = 0, nbit=8, npol=npol, ndim=ndim, ntimesamples=ntimesamples, ninputs=2, ninputs_xgpu=2, inttime_msec=inttime_msec, 
                                       proj_id = "SKA1", exptime_sec = exptime_sec, file_size=file_size, n_fine_channels=1, bandwidth_hz=bandwidth_hz,
                                       frequency_mhz=frequency_mhz, telescope="LFAASP", source=options.source, ra=options.ra_deg, dec=options.dec_deg )
        out_f = open(  hdrfile , "w" )
        out_f.write( header )
        out_f.close()    
        