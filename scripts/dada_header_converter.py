from __future__ import print_function
import sys
from optparse import OptionParser,OptionGroup


def parse_options(idx=0):
   usage="Usage: %prog [options]\n"
   usage+='\tGets basic information from the HDF5 file\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('-f','--frequency',dest="frequency",default=None, help="Final frequency [default %default MHz]",type="float")
   parser.add_option('-c','--channel',dest="channel",default=None, help="Final frequency channel [default %default]",type="int")
#   parser.add_option('-n','--n_chan','--n_channels','--fft_samples',dest="n_fft_samples",default=32, help="Number of FFT samples/channels [default %default]",type="int")
#   parser.add_option('-i','--inttime',dest="inttime",default=0.2, help="Requested integration time to calculate number of chunks per uvfits file [default %default sec]",type="float")
#   parser.add_option('-t','--type','--file_type','--filetype',dest="file_type",default="auto", help="Type of file if not specified -> automatic detection from file name [default %default]")
   
   (options, args) = parser.parse_args(sys.argv[idx:])

   return (options, args)

def convert_header_file( dada_header, out_file, options ) :

   file=open(dada_header,'r')
   data=file.readlines()
   out_f=open(out_file,"w")
   
   output_string=""

   for line in data :    
      words = line.split(' ')
      
      if words[0] in ["FREQ" , "BW" , "NCHAN" , "BANDWIDTH_HZ" , "NINPUTS" , "NINPUTS_XGPU" ] :
         if words[0] == "FREQ" :
            if options.frequency is not None :
               output_string += ("FREQ %.4f\n" % (options.frequency))
            else :
               output_string += line

         if words[0] == "BW" :
            output_string += "BW 0.9259\n"
            
         if words[0] == "NCHAN" :
            output_string += "NCHAN 1\n"         
         
         if words[0] == "BANDWIDTH_HZ" :
            output_string += "BANDWIDTH_HZ 925925\n"
            
         if words[0] == "NINPUTS" :
            output_string += "NINPUTS 2\n"
            
         if words[0] == "NINPUTS_XGPU" :
            output_string += "NINPUTS_XGPU 2\n"
                     
      else :
         output_string += line
         
         
   output_string_4096 = output_string[0:4096]         
   
   l = len(output_string_4096)
   print("INFO : header length = %d bytes" % (l))
   if l < 4096 :
      add_zeros = (4096-l)
      print("WARNING : header shorter than 4096 bytes -> adding %d zeros" % (add_zeros))
      last_char=output_string_4096[l-1]
      for i in range(0,add_zeros) :
         output_string_4096 += last_char
         

   l = len(output_string_4096)
   print("INFO : final header length = %d bytes" % (l))
   
   out_f.write(output_string_4096)      
   out_f.close()

def main() :   
   dada_header="header_template.txt"
   if len(sys.argv) > 0:
      dada_headar = sys.argv[1]
      
   out_file="header.txt"
   if len(sys.argv) > 1:
      out_file = sys.argv[2]
      
   (options, args) = parse_options(1)
   
   if options.frequency is None and options.channel is not None :
      options.frequency = options.channel * (400.00/512.00)
   
   print("######################################")
   print("PARAMETERS :")
   print("######################################")
   print("Dada header template = %s" % (dada_headar))
   print("Dada output header   = %s" % (out_file))
   if options.frequency is not None :
      print("Frequency            = %.4f MHz" % (options.frequency))
   else :
      print("WARNING : keeping the same frequency")
   print("######################################")
   
   convert_header_file( dada_headar, out_file , options )      

if __name__ == "__main__":
   main()
  
 