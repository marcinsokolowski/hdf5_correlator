// See also :
// h5dump program - part of linux or needs instalation of some hdf5 package 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <complex>

#include <H5Cpp.h>

// my :
#include <spectrometer.h>
#include <bg_date.h>
#include <bg_globals.h>
#include <cvalue_vector.h>
#include <calsol_values.h>

// 
#include <myfile.h>
#include <myparser.h>
#include <mystring.h>
#include <mystrtable.h>


#include "hdf5_commons.h"
// #include "eda1_tpm_coefficients.h"
#include "eda2_tpm_coefficients.h"
// #include "calsolutions.h"
// #include "mwaconfig.h"
// #include "defines.h"

int antenna1 = 0;
int antenna2 = 0;
int gDebugLevel = 0;
int gPol=0;

// saving statistics of zeros for hdf5 file :
string gZeroStatFile;

// global name for input hdf5 file :
string gInputHdf5Filename;


typedef struct {
    char re;  
    char im;  
} complex_t;


enum ePhaseNormalisation_t { eNoPhaseNorm=0, ePhaseNorm0_360=1, ePhaseNorm_m180_180=2 };
ePhaseNormalisation_t gPhaseNormalisation = ePhaseNorm0_360;

double gSign = 1.00;
int gMaxTimeSteps=-1;
int gDoPhaseCorrection=0;

int gRandomSeed = 477737;

// unix time :
int gFileUxTime = -1;

// number of samples for fine channalisation (FFT) :
int gNSamples=64;



/*// Sign to Test signs conventions :
double gSign = 1.00;    // sign of phase factor due to geometrical delay
double gCalSign = 1.00; // calibration sign

int gMaxTimeSteps=-1;

// if apply phase correction :
int gDoPhaseCorrection=0;

string gPointingString;
double gPointingAz_DEG = -1000;
double gPointingElev_DEG = -1000;

// saving data for Pulsar timing :
string gOutputTimeSeriesFileBin;
string gOutputTimeSeriesFileText;
string gOutputTimeSeriesFileBase;

*/

// current coefficients are for EDA2 and 48 antennas as for 20190610 :
#define gAntennaPhaseDiff gEDA2_AntennaPhaseDiff_MSok_Norm_m180_180_20190609


herr_t list_obj_iterate(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data)
{
        // string szTmp;
        // Beam2016Implementation* pBeamModelPtr = (Beam2016Implementation*) operator_data;  
        
//        if(pBeamModelPtr == nullptr)
//                throw std::runtime_error("The pointer to  Beam2016Implementation class in Beam2016Implementation::list_obj_iterate must not be null");

        if (name[0] == '.') {         /* Root group, do not print '.' */
        } else {
                switch (info->type) {
                    case H5O_TYPE_GROUP:
                    case H5O_TYPE_NAMED_DATATYPE:
                    default:
                        printf("GROUP : %s\n",name);
                        break;
                    case H5O_TYPE_DATASET:
                        // szTmp = name;
                        // pBeamModelPtr->m_obj_list.push_back(szTmp);
                        printf("DATASET : %s\n",name);
                        break;
                }
        }

        return 0;
}


void dump_data( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szOutFileName=NULL, int pol=0 )
{
   if( szOutFileName && strlen( szOutFileName) > 0 ){
      FILE* out_f = fopen( szOutFileName, "wb" );
      void* data_ptr = &(data[0]);
      // size_t fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
      long long written = fwrite( data_ptr, sizeof(data[0]), data.size(), out_f );
      double written_bytes = double(written)*sizeof(data[0]);
      double expected_write = double(data.size()) * sizeof(data[0]);
      fclose(out_f);
      printf("Written %.0f bytes to file %s, should write %d x %d = %d\n",written_bytes, szOutFileName, int(data.size()), int(sizeof(data[0])), expected_write );
      printf("Test values data[0] = %d / %d\n",data[0].re,data[0].im);
      printf("Test values data[1] = %d / %d\n",data[1].re,data[1].im);
      printf("Test values data[2] = %d / %d\n",data[2].re,data[2].im);
      printf("Test values data[10] = %d / %d\n",data[10].re,data[10].im);
   }else{
      printf("WARNING : no data dump to binary file is required - please verify !\n");
   }
   
   if( strlen( gZeroStatFile.c_str() ) > 0 ){
      std::vector<int> zeros_count , longest_zeros_in_row, current_zeros_in_row, start_of_current_zeros_in_row, start_of_longest_zeros_in_row;
      int n_inputs = n_ants*n_pols;
      zeros_count.assign( n_inputs , 0 );
      longest_zeros_in_row.assign( n_inputs , 0 );
      current_zeros_in_row.assign( n_inputs , 0 );
      start_of_current_zeros_in_row.assign( n_inputs , -1 );
      start_of_longest_zeros_in_row.assign( n_inputs , -1 );
            
            
      int time_count = data.size() / n_inputs;            
      bool bHeader = false;
      if ( !MyFile::DoesFileExist( gZeroStatFile.c_str() ) ){
         bHeader = true;
      }
      MyOFile out_f( gZeroStatFile.c_str() , "a+" );
      if( bHeader ){
         out_f.Printf(" # FILENAME MAX_ZERO_COUNT(on input) T1X T1Y T2X T2Y T3X T3Y T4X T4Y T5X T5Y T6X T6Y T7X T7Y T8X T8Y T9X T9Y T10X T10Y T11X T11Y T12X T12Y T13X T13Y T14X T14Y T15X T15Y T16X T16Y\n");
      }
      
      char szOutFileLongestChain[1024];
      sprintf(szOutFileLongestChain,"LongestChain_%s",gZeroStatFile.c_str());
      MyOFile out_f2( szOutFileLongestChain, "a+" );
      
      for(int t=0;t<time_count;t++){
         int index_offset = t*n_inputs;
       
         for(int inp=0;inp<n_inputs;inp++){
            if( data[index_offset+inp].re == 0 && data[index_offset+inp].im == 0 ){
               zeros_count[inp]++;
               current_zeros_in_row[inp]++;
               if( start_of_current_zeros_in_row[inp] < 0 ){
                  start_of_current_zeros_in_row[inp] = t;
               }
            }else{
               // reset current ZERO chain :
               current_zeros_in_row[inp] = 0;
               start_of_current_zeros_in_row[inp] = -1;
            }
           
            // update longest in row and start of longest in row:
            if( current_zeros_in_row[inp] >= longest_zeros_in_row[inp] ){
               longest_zeros_in_row[inp] = current_zeros_in_row[inp];
               start_of_longest_zeros_in_row[inp] = start_of_current_zeros_in_row[inp];
            }
         }  
      }
      
      int max_count = -1, max_count_inp = -1;
      for(int inp=0;inp<n_inputs;inp++){
         if( zeros_count[inp] > max_count ){
            max_count = zeros_count[inp];
            max_count_inp = inp;
         }
      }
      
      mystring outLine = gInputHdf5Filename.c_str();
      mystring outLine2 = gInputHdf5Filename.c_str();
      char szTmp[64];
      sprintf(szTmp," %d(%d) ",max_count,max_count_inp);
      outLine += szTmp;

      outLine += " ";
      outLine2 += " ";
      for(int inp=0;inp<n_inputs;inp++){
         sprintf(szTmp,"%d(%d) ",zeros_count[inp],longest_zeros_in_row[inp]);                  
         outLine += szTmp;         
         
         sprintf(szTmp,"%d(%d) ",longest_zeros_in_row[inp],start_of_longest_zeros_in_row[inp]);
         outLine2 += szTmp;
      }
      
      out_f.Printf("%s\n",outLine.c_str());
      
      out_f2.Printf("%s\n",outLine2.c_str());
   }
}

int get_ant_data(std::vector< complex_t >& data, int ant_idx, int n_ants, int n_pols, int n_samples, std::vector< std::vector< complex_t > >& out_ant_data, int pol=0 )
{
/*    for( std::vector< std::vector< complex_t > >::iterator it=out_ant_data.begin();it!=out_ant_data.end();it++){
        it->clear();
    }*/
    out_ant_data.clear();
    
    int ants_pols = n_ants * n_pols;    
    int start_index = n_pols*ant_idx + pol; // start at given antenna index (only X polarisation implemented yet ), was ant_idx which was ok if antenna was really antenna*2pols = 0, 2, 4, 6, 8 ...
    printf("get_ant_data : start_index = %d\n",start_index);
    while( start_index < data.size() ){
        std::vector< complex_t > samples;
        
        // add n_samples for fine channalisation :
        for(int i=0;i<n_samples;i++){
            int idx = start_index + i*ants_pols;
            
            if( idx < data.size() ){
                samples.push_back( data[ idx ] );
            }
        }    

                
        if( samples.size() == n_samples ){
            // add samples to output list only if matching the required number :
            out_ant_data.push_back( samples );
        }            
        
        start_index += n_samples*ants_pols;
    }
    
    return out_ant_data.size();
}


int get_ant_data(std::vector< complex_t >& data, int ant_idx, int n_ants, int n_pols, std::vector< complex_t >& out_ant_data, int pol=0, int max_timesteps=-1 )
{
    out_ant_data.clear();
    
    int ants_pols = n_ants * n_pols;    
    int idx = n_pols*ant_idx + pol; // start at given antenna index (only X polarisation implemented yet )
    printf("get_ant_data : start_index = %d\n",idx);
    while( idx < data.size() ){
        out_ant_data.push_back( data[ idx ] );
        
        idx += ants_pols;
        
        if( max_timesteps > 0 ){
           if( out_ant_data.size() >= max_timesteps ){
              break;
           }
        }
    }
    
    return out_ant_data.size();
}


// was to normalise to 0 - 360, but now changed to -180 - 180 
double normalise_phase_m180_180( double& phase_deg ){
     while( phase_deg < -180 ){
         phase_deg += 360.00;
     }

     while( phase_deg > 180.0 ){
         phase_deg -= 360.00;
     }
     
     return phase_deg;
}


double normalise_phase_0_360( double& phase_deg ){
     while( phase_deg < 0 ){
         phase_deg += 360.00;
     }

     while( phase_deg > 360.0 ){
         phase_deg -= 360.00;
     }
     
     return phase_deg;
}

/*
   enum ePhaseNormalisation_t { eNoPhaseNorm=0, ePhaseNorm0_360=1, ePhaseNorm_m180_180=2 };
   ePhaseNormalisation_t gPhaseNormalisation = ePhaseNorm0_360;
*/
double normalise_phase( double& phase_deg )
{
   switch( gPhaseNormalisation )
   {   
       case ePhaseNorm0_360 :
          return normalise_phase_0_360( phase_deg );
       case ePhaseNorm_m180_180 :
          return normalise_phase_m180_180( phase_deg );       
       default :
          return phase_deg;
   }

   return phase_deg;
}


int is_in_list( vector<int>& ant_list, int ant )
{
    for(int i=0;i<ant_list.size();i++){
        if( ant == ant_list[i] ){
            return i;
        }
    }

    return -1;
}

std::complex<double> phase_factor_func( double freq_mhz, double delay_picosec )
{
   double delay_sec = delay_picosec * 1e-12;
   double freq_hz   = freq_mhz * 1e6;
   double oversampling = (32.00/27.00);
   
   double phase_rad = gSign * 2.00*M_PI*freq_mhz*delay_picosec*(1e-6); // to avoid numerical errors 

   // test with respect to antenna 0 :
//   double phase_rad_ant0 = 2.00*M_PI*freq_mhz*(12099.42)*(1e-6);
//   phase_rad -= phase_rad_ant0;   
   
   
   std::complex<double> phase_coeff( cos(phase_rad), sin(phase_rad) );
//   std::complex<double> phase_coeff = std::exp(1i * phase_rad);
   
   return phase_coeff;
}

void correlate_antennas( std::vector< complex_t >& data, int n_ants, int n_pols, int n_samples, int pol=0 )
{
   std::vector< std::vector< complex_t > > data_per_ant1, data_per_ant2;   
   printf("Selecting data for antenna = %d\n",antenna1);fflush(stdout);
   int n_blocks1 = get_ant_data( data, antenna1, n_ants, n_pols, n_samples, data_per_ant1, pol );
   printf("Selected %d blocks of %d samples for antenna %d (total %d samples)\n",n_blocks1,n_samples,antenna1,(n_blocks1*n_samples));
   
   printf("Selecting data for antenna = %d\n",antenna2);fflush(stdout);
   int n_blocks2 = get_ant_data( data, antenna2, n_ants, n_pols, n_samples, data_per_ant2, pol );
   printf("Selected %d blocks of %d samples for antenna %d (total %d samples)\n",n_blocks2,n_samples,antenna2,(n_blocks2*n_samples));

   
   // initialise FFT  buffers:
   std::complex<float>* data1_in = new std::complex<float>[n_samples];
   std::complex<float>* data2_in = new std::complex<float>[n_samples];
   std::complex<float>* spectrum1_reim = new std::complex<float>[n_samples];
   std::complex<float>* spectrum2_reim = new std::complex<float>[n_samples];
   std::complex<float>* corr_reim = new std::complex<float>[n_samples];
   double* spectrum1 = new double[n_samples];
   double* spectrum2 = new double[n_samples];
   
   vector<double> mean_spectrum1,mean_spectrum2;
   vector< std::complex<float> > mean_corr;            // mean of RE/IM from correlation
   vector< double > mean_phase, sum2_phase, rms_phase; // mean and rms of the phase 
   int mean_count1=0,mean_count2=0,mean_corr_count=0;
   mean_spectrum1.assign( n_samples, 0.00 );
   mean_spectrum2.assign( n_samples, 0.00 );
   mean_corr.assign( n_samples, 0.00 );
   
   mean_phase.assign( n_samples, 0.00 );
   sum2_phase.assign( n_samples, 0.00 );
   rms_phase.assign( n_samples, 0.00 );

   if( n_blocks1 == n_blocks2 ){
      printf("OK : number of blocks from antennas %d and %d are the same -> ok\n",antenna1,antenna2);   
   }else{
      printf("ERROR : number of blocks from antennas %d and %d is different -> cannot correlate antennas\n",antenna1,antenna2);
      exit(-1);
   }


   int antenna1_index = antenna1; // now antenna1 is really antenna index not antenna_index*2pols !
   double delta_phase1_rad = gAntennaPhaseDiff[antenna1_index][1] * (M_PI/180.00);
   int antenna2_index = antenna2; // now antenna2 is really antenna index not antenna_index*2pols !
   double delta_phase2_rad = gAntennaPhaseDiff[antenna2_index][1] * (M_PI/180.00);

   printf("Antenna-%02d phase difference = %.2f [deg]\n",antenna1_index,delta_phase1_rad*(180.0/M_PI));
   printf("Antenna-%02d phase difference = %.2f [deg]\n",antenna2_index,delta_phase2_rad*(180.0/M_PI));
   sleep(1);

      
   int out_count1=0,out_count2=0;
   for(int block=0;(block<data_per_ant1.size() && (block<gMaxTimeSteps || gMaxTimeSteps<0));block++){
       if( (block % 1000) == 0 ){
           printf("Block %d / %d ...\n",block,(int)(data_per_ant1.size()));fflush(stdout);
       }
   
       std::vector< complex_t >& block_samples1 = data_per_ant1[block];
       std::vector< complex_t >& block_samples2 = data_per_ant2[block];

       for(int s=0;s<n_samples;s++){
           data1_in[s] = std::complex<float>( block_samples1[s].re , block_samples1[s].im );
           data2_in[s] = std::complex<float>( block_samples2[s].re , block_samples2[s].im );

           
 
           // corr_reim[s] = ( data1_in[s] ) * std::conj( data2_in[s] );
           
           if( block == 0 && gDebugLevel>=2 ){
               printf("DEBUG : %d / %d -> %.1f %.1f\n",block_samples1[s].re , block_samples1[s].im, data1_in[s].real() , data1_in[s].imag() );
           }
           
       }   
       
       double norm = 1.00;
       CSpectrometer::doFFT( data1_in, n_samples, spectrum1, spectrum1_reim, out_count1, norm );
       CSpectrometer::doFFT( data2_in, n_samples, spectrum2, spectrum2_reim, out_count2, norm );
       
       if( out_count1 != out_count2 ){
          printf("ERROR : outputs from doFFT have different number of samples %d != %d\n",out_count1,out_count2);
          exit(-1);
       }
       
       for(int ch=0;ch<out_count1;ch++){
           if(  gDoPhaseCorrection > 0 ){              
               std::complex<float> phase1_factor = std::complex<float>( cos(delta_phase1_rad), sin(delta_phase1_rad) );
               spectrum1_reim[ch] *= phase1_factor;

               std::complex<float> phase2_factor = std::complex<float>( cos(delta_phase2_rad), sin(delta_phase2_rad) );
               spectrum2_reim[ch] *= phase2_factor;                              
               
               printf("WARNING : phase calibration applied !!!\n");
           }
       
           corr_reim[ch] = ( spectrum1_reim[ch] ) * std::conj( spectrum2_reim[ch] );
       }
       
       for(int ch=0;ch<out_count1;ch++){
           mean_spectrum1[ch] += spectrum1[ch];
           mean_count1++;
       }
           
       for(int ch=0;ch<out_count2;ch++){
           mean_spectrum2[ch] += spectrum2[ch];
           mean_count2++;                      
       }
       
       for(int ch=0;ch<out_count2;ch++){
           mean_corr[ch] += corr_reim[ch];
           
           double re = corr_reim[ch].real();
           double im = corr_reim[ch].imag(); 
           // double phase_deg = atan2( im , re )*(180.00/M_PI);                      
           double phase_deg = std::arg( corr_reim[ch] ) * ( 180.00 / M_PI );
           // double phase_rad = std::arg( corr_reim[ch] );
           normalise_phase( phase_deg );
                      
           mean_phase[ch] += phase_deg;
           sum2_phase[ch] += phase_deg*phase_deg;


           if( ch == 0 && gDebugLevel >= 2 ){
               printf("DEBUG : phase_deg = %.4f [deg]\n",phase_deg);
           }
       }
       mean_corr_count++;              
   }
      
   delete [] data1_in;
   delete [] spectrum1;   
   delete [] spectrum1_reim;

   delete [] data2_in;
   delete [] spectrum2;   
   delete [] spectrum2_reim;

   char szOutFile1[128],szOutFile2[128],szOutCorr[128],szOutPhase[128];
   
   sprintf(szOutFile1,"test_spectrum%03d.txt",antenna1);
   sprintf(szOutFile2,"test_spectrum%03d.txt",antenna2);
   sprintf(szOutCorr,"corr_antenna_%03d_%03d.txt",antenna1,antenna2);
   sprintf(szOutPhase,"corr_phase_%03d_%03d.txt",antenna1,antenna2);
   
   FILE* out_f = fopen( szOutFile1 , "w");
   for(int ch=0;ch<out_count1;ch++){ // was (n_samples/2)
       mean_spectrum1[ch] / mean_count1;
       fprintf(out_f,"%d %.4f\n",ch,mean_spectrum1[ch]);
   }
   fclose(out_f);     
   
   out_f = fopen( szOutFile2 , "w");
   for(int ch=0;ch<out_count1;ch++){ // was (n_samples/2)
       mean_spectrum2[ch] / mean_count2;
       fprintf(out_f,"%d %.4f\n",ch,mean_spectrum2[ch]);
   }
   fclose(out_f);

   out_f = fopen( szOutCorr , "w" );
   FILE* out_phase_f = fopen( szOutPhase , "w");

   printf("Number of averaged spectra = %d\n",mean_corr_count);   
   fprintf(out_phase_f,"# CHANNEL <PHASE_DEG> RMS(PHASE_DEG) | OLD = PHASE( MEAN(RE/IM) )\n");
   for(int ch=0;ch<out_count1;ch++){ // was (n_samples/2)
       printf("%d : ( %.5f , %.5f ) / %d = ( %.5f , %.5f )\n",ch,mean_corr[ch].real(),mean_corr[ch].imag(),mean_corr_count,(mean_corr[ch].real() / mean_corr_count),(mean_corr[ch].imag() / mean_corr_count));
       mean_corr[ch] = std::complex<float>( mean_corr[ch].real() / mean_corr_count , mean_corr[ch].imag() / mean_corr_count );
              
       
       // std::complex<float> sum2_mean( sum2_corr[ch].real() / mean_corr_count , sum2_corr[ch].imag() / mean_corr_count );
       // rms_corr[ch]  = sqrt( sum2_mean - mean_corr[ch]*mean_corr[ch] );       
       
       // double phase_rad = atan2( mean_corr[ch].imag() , mean_corr[ch].real() );
       double phase_rad = std::arg( mean_corr[ch] );
       double phase_deg = (180.00/M_PI)*phase_rad;              
       normalise_phase( phase_deg );
       
       mean_phase[ch] = ( mean_phase[ch] / mean_corr_count );
       rms_phase[ch] = sqrt( ((sum2_phase[ch])/mean_corr_count) - (mean_phase[ch])*(mean_phase[ch]) );
   
       fprintf(out_f,"%d %.4f %.4f %.4f %.4f %.6f %.6f\n",ch,phase_deg,std::norm( mean_corr[ch] ),mean_corr[ch].real(),mean_corr[ch].imag(),mean_phase[ch],rms_phase[ch]);
       fprintf(out_phase_f,"%d %.4f %.4f | OLD = %4f\n",ch,mean_phase[ch],rms_phase[ch],phase_deg);
   }
   fclose(out_f);
   fclose(out_phase_f);
   
   delete [] corr_reim;

}

void usage()
{
   printf("read_hdf5_example HDF5_FILE -a ANT1 -b ANT2\n");
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "a:b:";   
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);
   
      switch (opt) {
         case 'a':
            antenna1 = atol( optarg );
            break;

         case 'b':
            antenna2 = atol( optarg );
            break;

         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
   
} 

void print_parameters()
{
   printf("##############################################\n");
   printf("PARAMETERS:\n");
   printf("##############################################\n");
   printf("Antenna1 = %d\n",antenna1);
   printf("Antenna2 = %d\n",antenna2);
   printf("##############################################\n");
   
   // 
   printf("Usage example [0][1] = %.4f , [15][1] = %.4f\n",gAntennaPhaseDiff[0][1],gAntennaPhaseDiff[15][1]);
}

int main(int argc,char* argv[])
{
   std::string filename = "channel_cont_0_20190307_25356_0.hdf5";
   if ( argc>=2 ){
       filename = argv[1];
   }
   gInputHdf5Filename = filename.c_str();
   parse_cmdline( argc , argv );
   print_parameters();
   
   // random number initialisation :
   srand( gRandomSeed );


   H5::H5File file( filename.c_str(), H5F_ACC_RDONLY ); 
   gFileUxTime = hdf5filename2uxtime( filename.c_str() );
   
//   H5::DataSet dataset = file.openDataSet( "chan" );
//   H5::H5Group group = H5Gopen(file_id, "/chan_", H5P_DEFAULT);
//   H5::H5Attribute attr = openAttribute

   int count = file.getObjCount();
   hid_t file_id = file.getId();
   printf("Number of objects in file %s is %d , file_id = %x\n",filename.c_str(), count, (int)file_id );
   herr_t status =  H5Ovisit (file_id, H5_INDEX_NAME, H5_ITER_NATIVE, list_obj_iterate, NULL);
   if( status < 0 ) {
       printf("ERROR\n");
       // throw std::runtime_error("H5Ovisit returned with negative value which indicates a critical error");
   }

   H5::DataSet dataset = file.openDataSet( "chan_/data" );
   printf("Successfully opened dataset chan_/data\n");
   
   H5::DataSpace dataspace = dataset.getSpace();
   int rank = dataspace.getSimpleExtentNdims();
   
   hsize_t dims_out[2];
   dataspace.getSimpleExtentDims( dims_out, NULL);
   dataspace.selectAll();
   H5T_class_t type_class = dataset.getTypeClass();
   printf("Data dimensions = %d x %d, data_type = %d\n",(int)(dims_out[0]),(int)(dims_out[1]),type_class);
   
   H5::DataType data_type = dataset.getDataType();
   hid_t data_type_id = data_type.getId();
   hid_t native_data_type_id = H5Tget_native_type( data_type_id , H5T_DIR_ASCEND);
//   printf("Data type = %s\n", (const char*)data_type.getTag() );
   printf("Data type size = %d vs. sizeof(complex_t) = %d , id = %d , native_data_type_id = %d\n",(int)(data_type.getSize()),(int)sizeof(complex_t),(int)data_type_id,(int)native_data_type_id);
   
//   hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof (complex_t));
//   H5Tinsert (complex_id, "real", HOFFSET(complex_t,re), H5T_NATIVE_CHAR);
//   H5Tinsert (complex_id, "imag", HOFFSET(complex_t,im), H5T_NATIVE_CHAR);
        
   std::vector< complex_t > data(dims_out[0]*dims_out[1]);
   printf("data.size = %d , sizeof(data[0]) = %d\n",(int)data.size(),(int)sizeof(data[0]));
   H5::DataSpace memspace( rank, dims_out );
//   dataset.read( data.data(), H5::PredType::NATIVE_CHAR, memspace, dataspace );

   H5::CompType ctype(sizeof(complex_t));
   ctype.insertMember( "real", HOFFSET(complex_t,re), H5::PredType::NATIVE_CHAR );
   ctype.insertMember( "imag", HOFFSET(complex_t,im), H5::PredType::NATIVE_CHAR );
   dataset.read( data.data(), ctype , memspace, dataspace );

   if( gDebugLevel >= 2 ){   
       for(int i=0;i<dims_out[1];i++){
           printf("%d / %d\n",data[i].re,data[i].im);
       }   
   }
   

   int n_samples = gNSamples; // number of samples to use for fine channalisation 
   int n_fine_channels = n_samples / 2;
   int ants_pols = dims_out[1];
   int n_ants = dims_out[1] / 2; // 32 / 2pols -> 16 antennas 
   int n_pols = 2;
   int samples_all = dims_out[0];
   int samples_per_pol = samples_all / ants_pols;
   
   printf("Action = correlate_antennas( data , %d , %d , %d, %d )\n",n_ants, n_pols, n_samples, gPol );
   correlate_antennas( data, n_ants, n_pols, n_samples, gPol );

//   H5::FileAccPropList prop_list = file.getAccessPlist();
}


