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
// #include <spectrometer.h>
// #include <bg_date.h>
// #include <bg_globals.h>

// 
// #include <myparser.h>
// #include <mystring.h>
// #include <mystrtable.h>

#include "hdf5_commons.h"
#include "eda1_tpm_coefficients.h"
#include "eda2_tpm_coefficients.h"

// current coefficients are for EDA2 and 48 antennas as for 20190610 :
#define gAntennaPhaseDiff gEDA2_AntennaPhaseDiff_MSok_Norm_m180_180_20190609

typedef struct {
    char re;  
    char im;  
} complex_t;


enum eActionType_t { eCorrelateAntennas=1, eBeamform=2, eDumpData=3, eQuickCal=4, eOptimiseQuickCal=5, eBeamformTest=6 };

enum ePhaseNormalisation_t { eNoPhaseNorm=0, ePhaseNorm0_360=1, ePhaseNorm_m180_180=2 };
ePhaseNormalisation_t gPhaseNormalisation = ePhaseNorm0_360;

int antenna1 = 0;
int antenna2 = 1;
int gDebugLevel = 0;
int gPol=0;
string gOutFileName;
vector<int> gAntennaListToProcess;
string      gAntennaListStr;

double gTestPhaseDeg = -1e8;
int gUseRandomPhase=0;
int gMaxTimeSteps=-1;

int gRandomSeed = 477737;

// if apply phase correction :
int gDoPhaseCorrection=0;
int gDoCalibrate=0;

// number of samples for fine channalisation (FFT) :
int gNSamples=64;

// unix time :
int gFileUxTime = -1;

// pointing parameters :
mystring gMeanDelaysString;
vector<double> gMeanDelays;

// phase offsets instead of pointing delays from HydA calibration for example 
mystring gPhaseOffsetsString;
vector<double> gPhaseOffsets;

// calibration phases :
mystring gCalPhaseOffsetsString;

// reference antenna for delays or phase offsets :
int gRefAnt=-1;

// Sign to Test signs conventions :
double gSign = 1.00;    // sign of phase factor due to geometrical delay
double gCalSign = 1.00; // calibration sign

// 
double gFreqMHz = (400.00/512.00)*204.00;

/*typedef struct {
    char re;  
    char im;  
} complex_t;*/

eActionType_t gActionType = eCorrelateAntennas;
int gIterations=0;
double gOptimisePhaseStep = 1.00;
double gOptimiseRadius    = -1.00; // disabled

string gSlopeFitFile;
vector<cValue> gSlopeParameters;

// execute in loop:
int gExecuteLoop = 0;

/* 
Mapping of TPM-23 inputs to EDA tiles :
see also : /home/msok/bighorns/software/analysis/scripts/shell/eda/eda2/pointing.py and ~/Desktop/EDA2/logbook/eda2_pointing.odt

Antenna 0 , RX 1    EDA-tile 0 
Antenna 1 , RX 2    EDA-tile 1
Antenna 2 , RX 3    EDA-tile 2 
Antenna 3 , RX 4    EDA-tile 3 
Antenna 8 , RX 5    EDA-tile 4 
Antenna 9 , RX 6    EDA-tile 5 
Antenna 10 , RX 7   EDA-tile 6 
Antenna 11 , RX 8   EDA-tile 7 
Antenna 15 , RX 9   EDA-tile 8 
Antenna 14 , RX 10  EDA-tile 9 
Antenna 13 , RX 11  EDA-tile 10 
Antenna 12 , RX 12  EDA-tile 11 
Antenna 7 , RX 13   EDA-tile 12
Antenna 6 , RX 14   EDA-tile 13
Antenna 5 , RX 15   EDA-tile 14
Antenna 4 , RX 16   EDA-tile 15
*/


// QUICK CAL :
void beamform_fit_phase_offsets( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, const char* szOutFileName=NULL, int pol=0, int n_iter=0 );
void beamform_optimise_phase_offsets( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, std::vector< double >& phase_offsets, const char* szOutFileName=NULL, int pol=0, int iteration=0, double phase_step = 1.00 );


// EDA1 / TPM mapping :
int rx2tile_mapping[16] = { 0, 1 , 2 , 3 , 8 , 9 , 10 , 11 , 15 , 14 , 13 , 12 , 7 , 6 , 5 , 4 }; 
int aavsant2edaindex_mapping[16] = { 0, 1, 2, 3, 15, 14, 13, 12, 4, 5, 6, 7, 11, 10, 9, 8 };

int tpmrx2edatileindex( int eda_tile_index )
{
//   int eda_tile_index = (rx-1);
   
   return rx2tile_mapping[eda_tile_index];   
}

int aavsant2edaindex( int aavs_index )
{
    return aavsant2edaindex_mapping[ aavs_index ];
}

int get_random( int lowerbound=0, int upperbound=15 )
{
   double r = double(rand())/double(RAND_MAX+1.0);
   int ret = lowerbound+(int)((upperbound - lowerbound)*r);
   
   if( ret >= 0 && ret <= upperbound ){
      return ret;
   }
   
   printf("ERROR ret = %d - outside of %d - %d range\n",ret,lowerbound,upperbound);   
   return 0;
}

string ListToStr(  vector<int>& list )
{
   string out_list;
   char szTmp[16];
   for(int i=0;i<list.size();i++){
      sprintf(szTmp,"%d",list[i] );
   
      out_list += szTmp;
      out_list += ",";
   }
   
   return out_list;
}

int ParseCommaList( char* szList, vector<int>& out_list, const char * sep )
{
   char* saveptr=NULL;
   char* str1=szList;   
   
   for (str1 = szList; ; str1 = NULL) {
      char* token = strtok_r(str1, sep, &saveptr);
      if (token == NULL)
         break;

      // trim right :         
      char szToken[64];
      strcpy(szToken,token);
      int i=strlen(szToken)-1;
      while(szToken[i]==' '){
         szToken[i]='\0';
         i--;
      }   
      
      char* ptr = szToken;
      i=0;
      while(ptr[0]==' ' && i<strlen(szToken)){
         ptr++;
         i++;
      }   


      out_list.push_back(atol(ptr));
   }
      
   return out_list.size();
}
 



const char* get_action_desc( eActionType_t action )
{
    switch ( action )     
    {
       case eCorrelateAntennas :
           return "Correlate antennas";
           break;
      
       case eBeamform :
           return "Beamform";
           break;

       case eBeamformTest :
           return "BeamformTest";
           break;
           
       case  eDumpData :
          return "Dump_Data";
          break;
       
       case  eQuickCal :
          return "QuickCal";
          break;

       case  eOptimiseQuickCal :
          return "Optimise_QuickCal";
          break;
       
       default :
           return "unknown";
    }
    
    return "Unknown";
}

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

void beamform( std::vector< complex_t >& data, int n_ants, int n_pols, int n_samples, const char* szInFileName, const char* szOutFileName=NULL, int pol=0, std::vector<double>* p_phase_offsets=NULL )
{
   int do_phase_correction = gDoPhaseCorrection;
   int do_calibrate        = gDoCalibrate;
   double calibrate_sign   = gCalSign;

   if( p_phase_offsets ){
      mystring szTestOffsets;
      for(int i=0;i<p_phase_offsets->size();i++){
          gAntennaPhaseDiff[i][1] = (*p_phase_offsets)[i];

          szTestOffsets << (*p_phase_offsets)[i] << " , ";
      }

      gUseRandomPhase = 0;
      gTestPhaseDeg = -1e8;
      gMeanDelays.clear();
      gPhaseOffsets.clear();
      do_phase_correction = 1;
      do_calibrate = 1;
      calibrate_sign = 1.00;

      printf("INFO : using externally provided set of %d phases : %s\n",int(p_phase_offsets->size()),szTestOffsets.c_str());
   }



   if( !szOutFileName || strlen(szOutFileName) == 0 ){
      szOutFileName = "beamformed.txt";
   }

   std::vector< std::vector< complex_t > >* antenna_data[MAX_ANTS]; // MAX_ANTS does not work here ???
   for(int i=0;i<MAX_ANTS;i++){ // MAX_ANTS does not work here ???
       antenna_data[i] = NULL;
   }
   
   vector< std::complex<double> > beamformed_data_accum; // beamformed data - coherently added 16 EDA tiles 
   vector<double> mean_spectrum;    
   int beamformed_count = 0;
   std::complex<float>* data_in = new std::complex<float>[n_samples];
    
   beamformed_data_accum.assign( n_samples/2, 0.00 );
   mean_spectrum.assign( n_samples, 0.00 );

   std::complex<float>* spectrum_reim = new std::complex<float>[n_samples];
   double* spectrum = new double[n_samples];
   
   
   int ants_to_process = n_ants; // mainly for debugging when I wanted to only check single antenna : 
   int block_count = -1;
   for(int ant=0;ant<ants_to_process;ant++){
        if( gAntennaListToProcess.size() ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                printf("WARNING : antenna %d not in the list -> skipped\n",ant);
                continue;
            }
        }
   
        std::vector< std::vector< complex_t > >* data_per_ant = new std::vector< std::vector< complex_t > >;
        antenna_data[ant] = data_per_ant;

        printf("Selecting data for antenna = %d\n",ant);fflush(stdout);
        int n_blocks = get_ant_data( data, ant, n_ants, n_pols, n_samples, (*data_per_ant), pol );
        printf("Selected %d blocks of %d samples for antenna %d (total %d samples)\n",n_blocks,n_samples,antenna1,(n_blocks*n_samples));                
        
        // check if all antennas have the same number of blocks -> otherwise error !
        if( block_count > 0 ){
            if( data_per_ant->size() == block_count ){
                printf("\tNumber of blocks = %d -> ok (same as all the others)\n",(int)(data_per_ant->size()));
            }else{
                printf("\tERROR : number of blocks = %d != previous ( = %d )\n",(int)(data_per_ant->size()),block_count);
                exit(-1);
            }
        }else{
            block_count = data_per_ant->size();
        }
    }    


    vector< std::complex<double> > beamformed_data; // beamformed data of all antennas (or EDA tiles) for given time step
    beamformed_data.assign( n_samples, 0.00 );

    for(int block=0;block<block_count;block++){
        if( (block % 1000) == 0 ){
           printf("Block %d / %d ...\n",block,block_count);fflush(stdout);
        }

        // clear values to ZERO :    
        beamformed_data.assign( n_samples, 0.00 );
        
        for(int ant=0;ant<ants_to_process;ant++){
            if( gAntennaListToProcess.size() ){
                if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                    printf("WARNING : antenna %d not in the list -> skipped\n",ant);
                    continue;
                }
            }

        
            double delta_phase_rad = gAntennaPhaseDiff[ant][1] * (M_PI/180.00);
            /*if( gPhaseOffsets.size() > 0 ){
                if( ant < gPhaseOffsets.size() ){
                   // MEAN delays were calculated for EDA tiles based on layout and setup in the pointing scripts. This is not how they are connected to the TPM 
                   // so mapping has to be applied here to convert from TPM input Antenna index to actual EDA tile index to get the proper delays
                   // int eda_tile_index = tpmrx2edatileindex( ant );

                   // WARNING : here we assume phase offsets to be passed by ANTENNA from images (no mapping required) :
                   int eda_tile_index = ant;

                   delta_phase_rad = gSign * gPhaseOffsets[eda_tile_index] * (M_PI/180.00);
                 }else{
                      printf("ERROR : number of phase offsets = %d , requesting delay index %d -> check number of parameters provided to -P option, should be %d\n",int(gPhaseOffsets.size()),ant,ants_to_process);
                      exit(-1);
                 }
            }*/

            
            if( block==0 ){
                printf("Antenna-%02d phase difference = %.2f [deg]\n",ant,delta_phase_rad*(180.0/M_PI));                
            }
            std::vector< std::vector< complex_t > >& data_per_ant = *(antenna_data[ant]);
            
            std::vector< complex_t >& block_samples = data_per_ant[block];
            memset( data_in, '\0', sizeof(data_in[0])*n_samples );
            for(int s=0;s<n_samples;s++){
               data_in[s] = std::complex<float>( block_samples[s].re , block_samples[s].im );
            }
            
            int out_count=0;
            double norm = 1.00;
            
            // WARNING : this causes massive memory leaks !!!
            CSpectrometer::doFFT( data_in, n_samples, spectrum, spectrum_reim, out_count, norm );
            
            if(  gDoPhaseCorrection > 0 ){
               if( block == 0 ){
                   printf("Progress : applying beamforming phases to antenna = %d ( using delta_phi = %.2f rad = %.2f deg ), fine_channels = %d\n",ant,delta_phase_rad,delta_phase_rad*(180.00/M_PI),out_count);
               }
               std::complex<float> phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) );


               for(int ch=0;ch<out_count;ch++){
                   if( gSlopeParameters.size() ){
                       // WARNING : it does not help, I've already tested using just the central channel and the result was the same !
                       //           the problem is not really in the slope here ...
                       cValue& fit = gSlopeParameters[ant];
                       double delta_phase_deg = fit.y + fit.z * ch;
                       delta_phase_rad = gCalSign * delta_phase_deg * (M_PI/180.00);
                       phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) );
                       
                       if( block == 0 && ch==(out_count/2) ){
                          printf("\tFine channelised phase at the center %d channel = %.4f [deg]\n",(out_count/2),delta_phase_deg);
                       }
                   }
               
                   spectrum_reim[ch] *= phase_factor;
               }
            }   
            
            // add to beamformed data 
            for(int ch=0;ch<out_count;ch++){
                beamformed_data[ch] += spectrum_reim[ch];
            }
        }
 
        // FILE* out_f = fopen("test.txt","a+");
        int middle_ch = beamformed_data.size()/2;
        for(int ch=0;ch<beamformed_data.size();ch++){
             beamformed_data[ch] = std::complex<double>( beamformed_data[ch].real() / n_ants , beamformed_data[ch].imag() / n_ants );
             
             beamformed_data_accum[ch] += beamformed_data[ch];             
             
             double power = std::norm( beamformed_data[ch] ); // was abs
             
             mean_spectrum[ch] += power;
             
             //if( ch == 7 ){
             //   fprintf(out_f,"%d %.4f\n",beamformed_count,power);
             //}
        }     
        // fclose(out_f);
        beamformed_count++;
 
    }
    
    FILE* beamformed_f = fopen( szOutFileName,"a+");
    double total_power = 0.00;
    for(int ch=0;ch<beamformed_data_accum.size();ch++){
         double mag = std::norm( beamformed_data_accum[ch] ); // was abs
//        double mag = beamformed_data_accum[ch].real()*beamformed_data_accum[ch].real() + beamformed_data_accum[ch].imag()*beamformed_data_accum[ch].imag();        
        double phase_deg = std::arg(  beamformed_data_accum[ch] );        
        total_power += mean_spectrum[ch];

        if( ch == 16 ){    
//           fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,mag,phase_deg,mean_spectrum[ch]/beamformed_count,gFileUxTime);
//           fprintf(beamformed_f,"%d %.4f %.4f %.4f\n",ch,mag,phase_deg,mean_spectrum[ch]/beamformed_count);
        }
    }
    // test with total power over all channels
    printf("DEBUG : beamformed_data_accum.size() = %d\n",beamformed_data_accum.size());     
    fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,0.00,0.00,total_power,gFileUxTime);
    fclose( beamformed_f );

    // clean memory :
    for(int ant=0;ant<ants_to_process;ant++){
        if( antenna_data[ant] ){
            delete antenna_data[ant];
        }
    }    
    delete [] data_in;
    delete [] spectrum_reim;
    delete [] spectrum;
}


double beamform2( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, const char* szOutFileName=NULL, int pol=0, 
                vector< std::complex<double> > * p_beamformed_data_out=NULL, 
                std::vector<double>* p_phase_offsets=NULL,
                int debug_level=1, 
                int do_average_mean = 1 )
{  
   // to read all samples of given antenna 
   int time_steps = data.size() / (n_ants*n_pols);
   printf("beamform2 start : Number of expected time steps = %d\n",time_steps);
   if( gMaxTimeSteps > 0 ){
       printf("WARNING : number of time steps overwritten with external %d value\n",gMaxTimeSteps);
       time_steps = gMaxTimeSteps;
   }
   
   int do_phase_correction = gDoPhaseCorrection;
   int do_calibrate        = gDoCalibrate;
   double calibrate_sign   = gCalSign;
   
   if( p_phase_offsets ){
      mystring szTestOffsets;
      for(int i=0;i<p_phase_offsets->size();i++){
          gAntennaPhaseDiff[i][1] = (*p_phase_offsets)[i];
          
          szTestOffsets << (*p_phase_offsets)[i] << " , ";
      }
      
      gUseRandomPhase = 0;
      gTestPhaseDeg = -1e8;
      gMeanDelays.clear();
      gPhaseOffsets.clear();
      do_phase_correction = 1;
      do_calibrate = 1;
      calibrate_sign = 1.00;
      
      printf("INFO : using externally provided set of %d phases : %s\n",int(p_phase_offsets->size()),szTestOffsets.c_str());
   }      
   
   if( !szOutFileName || strlen(szOutFileName) == 0 ){
      szOutFileName = "beamformed2.txt";
   }

//   vector< std::vector< complex_t > > antenna_data;  // MAX_ANTS does not work here ???
//   antenna_data.reserve( n_ants );
   std::vector< complex_t > antenna_data[MAX_ANTS]; // 48 works ok MAX_ANTS not !
   
   std::complex<double> beamformed_data_accum; // beamformed data - coherently added 16 EDA tiles 
   double mean_spectrum = 0.00;    
   int beamformed_count = 0;
    
   int ants_to_process = n_ants; // mainly for debugging when I wanted to only check single antenna : 
   int block_count = -1;
   int n_ant_processed = 0;
   for(int ant=0;ant<ants_to_process;ant++){ // loop over AAVS antenna 
        if( gAntennaListToProcess.size() ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                printf("WARNING : antenna %d not in the list -> skipped\n",ant);
                continue;
            }
        }
   
        printf("Selecting data for antenna = %d\n",ant);fflush(stdout);
        int n_blocks = get_ant_data( data, ant, n_ants, n_pols, antenna_data[ant], pol, gMaxTimeSteps );
        printf("Selected %d blocks of %d samples for antenna %d (total %d samples), selected number of samples = %d ( data[time=0] = %d / %d )\n",n_blocks,time_steps,ant,(n_blocks*time_steps),(int)(antenna_data[ant].size()),antenna_data[ant][0].re,antenna_data[ant][0].im);
        
        if( time_steps != antenna_data[ant].size() ){
            printf("ERROR in code : expected %d time_steps, but read %d -> cannot continue\n",time_steps,(int)(antenna_data[ant].size()));
            exit(-1);
        }        
        std::vector< complex_t >& ant_data = antenna_data[ant];
//        printf("DEBUG1 : ant = %d, t = %d / %d -> vis = %d / %d\n",ant,0, ant_data.size(), ant_data[0].re , ant_data[0].im );
        
        
        if(  do_phase_correction > 0 ){
           double delta_phase_rad = 0.00;
           if( do_calibrate > 0 ){
              delta_phase_rad = gAntennaPhaseDiff[ant][1] * (M_PI/180.00);
              //double delta_phase_rad = 0.00;
              if ( gTestPhaseDeg > -1e5 ){
                  delta_phase_rad = gTestPhaseDeg*(M_PI/180.00); // TEST 
              }
              if( gUseRandomPhase > 0 ){
                 if( gUseRandomPhase == 1 ){
                     int random_idx = get_random(0, 15);
                     delta_phase_rad = gAntennaPhaseDiff[random_idx][1] * (M_PI/180.00);
                     printf("Randomised phase difference = %.2f [deg] from antenna = %d used\n",gAntennaPhaseDiff[random_idx][1],random_idx);
                 }
                 if( gUseRandomPhase == 2 ){
                     double random_phase = get_random(0,360);
                     random_phase = random_phase - 180;
                     printf("Randomised phase difference = %.1f [deg] used\n",random_phase);
                     delta_phase_rad = random_phase * (M_PI/180.00);
                 }
              }
              
              delta_phase_rad = calibrate_sign * delta_phase_rad;
              printf("Progress : applying beamforming phases to antenna = %d ( using delta_phi = %.2f rad = %.2f deg\n",ant,delta_phase_rad,delta_phase_rad*(180.00/M_PI));
           }
           

           std::complex<float> phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) ); // should be opposite sign to what I got from "calibration" as I cross-correlated 
                                                                                                                 // ANT_0 and ANT_N which means multiplied Vis_0 * Conj(Vis_N) and took phase of this 
                                                                                                                 // so now to "undo" the phase I need to multiply by exp(-i * phi)

//           if( ant_data.size() ){
//               printf("DEBUG2 : ant = %d, t = %d / %d -> vis = %d / %d\n",ant,0, ant_data.size(), ant_data[0].re , ant_data[0].im );
//           }

           for(int t=0;t<ant_data.size();t++){
              std::complex<float> tmp_complex = std::complex<float>( ant_data[t].re , ant_data[t].im );
              tmp_complex *= phase_factor;
              
              if( gMeanDelays.size() > 0 ){
                  if( ant < gMeanDelays.size() ){
                      // MEAN delays were calculated for EDA tiles based on layout and setup in the pointing scripts. This is not how they are connected to the TPM 
                      // so mapping has to be applied here to convert from TPM input Antenna index to actual EDA tile index to get the proper delays
//                      int eda_tile_index = tpmrx2edatileindex( ant );
//                      int eda_tile_index = ant;
                      int eda_tile_index = aavsant2edaindex( ant );
                  
                      std::complex<double> pointing_phase_coeff = phase_factor_func( gFreqMHz, gMeanDelays[eda_tile_index] );
                  
                      tmp_complex *= pointing_phase_coeff;
                  }else{
                      printf("ERROR : number of mean delays = %d , requesting delay index %d -> check number of parameters provided to -D option, should be %d\n",int(gMeanDelays.size()),ant,ants_to_process);
                      exit(-1);
                  }
              }else{
                 // only allow -D or -P options 
                 if( gPhaseOffsets.size() > 0 ){
                     if( ant < gPhaseOffsets.size() ){
                      // MEAN delays were calculated for EDA tiles based on layout and setup in the pointing scripts. This is not how they are connected to the TPM 
                      // so mapping has to be applied here to convert from TPM input Antenna index to actual EDA tile index to get the proper delays
                      // int eda_tile_index = tpmrx2edatileindex( ant );

                      // WARNING : here we assume phase offsets to be passed by ANTENNA from images (no mapping required) :
                      int eda_tile_index = ant;
                  
                      double phase_rad = gSign * gPhaseOffsets[eda_tile_index] * (M_PI/180.00);
                      std::complex<double> pointing_phase_coeff = std::complex<double>( cos(phase_rad) , sin(phase_rad) );
    
                      tmp_complex *= pointing_phase_coeff;
                  }else{
                      printf("ERROR : number of phase offsets = %d , requesting delay index %d -> check number of parameters provided to -P option, should be %d\n",int(gPhaseOffsets.size()),ant,ants_to_process);
                      exit(-1);
                  }
              }
              }
              
              ant_data[t].re = tmp_complex.real();
              ant_data[t].im = tmp_complex.imag();
           }
        }   
        
        n_ant_processed++;
    }
        


    vector< std::complex<double> > beamformed_data; // beamformed data of all antennas (or EDA tiles) for given time step
    beamformed_data.assign( time_steps, 0.00 );

    for(int time_step=0;time_step<time_steps;time_step++){
        if( debug_level > 0 ){
           if( (time_step % 1000) == 0 ){
              printf("Time step %d / %d ...\n",time_step,time_steps);fflush(stdout);
           }
        }

        for(int ant=0;ant<ants_to_process;ant++){
            if( gAntennaListToProcess.size() ){
                if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                    // printf("WARNING : antenna %d not in the list -> skipped\n",ant);
                    continue;
                }
            }

        
            std::vector< complex_t >& ant_data = antenna_data[ant];
            
            // add to beamformed data 
            beamformed_data[time_step] += std::complex<float>( ant_data[time_step].re , ant_data[time_step].im );
        }
    }

    for(int time_step=0;time_step<beamformed_data.size();time_step++){
       if( do_average_mean > 0 ){
          beamformed_data[time_step] = std::complex<double>( beamformed_data[time_step].real() / n_ants , beamformed_data[time_step].imag() / n_ants ); 
       }
//       beamformed_data[time_step] = std::complex<double>( beamformed_data[time_step].real() , beamformed_data[time_step].imag() ); 
             
       beamformed_data_accum += beamformed_data[time_step];
             
       double power = std::norm( beamformed_data[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] );
// TEST        double power = std::abs( beamformed_data[time_step] );
       mean_spectrum += power;
    }     
    
    FILE* beamformed_f = fopen( szOutFileName,"a+");
    double mag = std::norm( beamformed_data_accum ) / beamformed_data.size(); // was abs
//        double mag = beamformed_data_accum[ch].real()*beamformed_data_accum[ch].real() + beamformed_data_accum[ch].imag()*beamformed_data_accum[ch].imag();        
    double phase_deg = std::arg(  beamformed_data_accum );        

    double final_power = mean_spectrum/beamformed_data.size();    
//    double final_power = mean_spectrum;
    fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,mag,phase_deg,final_power,gFileUxTime);
    fclose( beamformed_f );
    
    printf("Saved results to %s file , processed %d antennas\n",szOutFileName,n_ant_processed);
    
    if( p_beamformed_data_out ){
       p_beamformed_data_out->clear();
       
       for(int time_step=0;time_step<beamformed_data.size();time_step++){
           p_beamformed_data_out->push_back( beamformed_data[time_step] );
       }
    }
    
    return mean_spectrum;
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
   
   sprintf(szOutFile1,"test_spectrum%02d.txt",antenna1);
   sprintf(szOutFile2,"test_spectrum%02d.txt",antenna2);
   sprintf(szOutCorr,"corr_antenna_%02d_%02d.txt",antenna1,antenna2);
   sprintf(szOutPhase,"corr_phase_%02d_%02d.txt",antenna1,antenna2);
   
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
   printf("hdf5_correlator HDF5_FILE -a ANT1 -b ANT2 -f -o OUTFILE_NAME -l ANTENNA_LIST -p O/1 -n PHASE_NORMALISATION -g TASK\n");
   printf("\n\nDefault action = %d (%s)\n\n",gActionType,get_action_desc(gActionType));
   printf("\t-f : enables phase correction based on coefficients in the table\n");
   printf("\t-l ANTENNA_LIST : coma separated antenna list (default ALL)\n");
   printf("\t-p POL : polarisation [default %d]\n",gPol);   
   printf("\t-n PHASE_NORM : phase normalisation, 0 : no norm, 1 : 0-360 deg, 2 : -180 - 180 deg [default %d]\n",gPhaseNormalisation);
   printf("\t-x MAX_TIMESTEPS : maximum number of timesteps to use [default %d]\n",gMaxTimeSteps);
   printf("\t-B : enabled beamforming\n");
   printf("\t-r TYPE : use random phase, TYPE = 1 - of different antenna in the set, 2 - just random value from -180 - 180 deg range\n");
   printf("\t-D 0,0,0,0,0,0,0,0,0,0,0,0 : list of mean delays in picoseconds to be applied to point in a particular direction\n");
   printf("\t-P 0,0,0,0,0,0,0,0,0,0,0,0 : list of phase offsets [in degrees] (this uses beamform function which can apply a slope)\n");
   printf("\t-X list_of_phases          : list of calibration phases [in degrees] (use this to provide list of calibration coefficients for each antenna)\n");
   printf("\t-S MULTIPLIER : mainly to test sign convention, but it will scale phase offsets by this number [default %.2f - SHOULD BE ONE !]\n",gSign);
   printf("\t-C MULTIPLIER : mainly to test sign convention for cal. soluations, but it will scale phase offsets by this number [default %.2f - SHOULD BE ONE !]\n",gCalSign);
   printf("\t-R REF_ANT    : delays or offsets with respect to antenna REF_ANT [default %d]\n",gRefAnt);
   printf("\t-d dump data to file, use -o to specify file name\n");
   printf("\t-c DO_CALIBRATE [default %d]\n",gDoCalibrate);
   printf("\t-g TASK : quick_cal, opt_quick_cal\n");
   printf("\t-i NUMBER_OF_ITERATIONS [default %d]\n",gIterations);
   printf("\t-s OPTIMISE_PHASE_STEP [default %.2f degrees]\n",gOptimisePhaseStep);
   printf("\t-z OPTIMISE RADIUS (around current best phase offset) [default %.2f < 0 -> full range -180 - 180 degrees]\n",gOptimiseRadius);
   printf("\t-F SLOPE_FIT_FILE - 3 column file with : | ANT INTERCEPT SLOPE | fits [ default not specified / ignored ]\n");
   printf("\t-L execute task in loop using the first antenna and all the other antennas [default not enabled = %d]\n",gExecuteLoop);
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "a:b:Bo:l:t:p:n:fx:r:D:P:S:c:R:C:dg:X:i:s:z:TF:N:L";
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

         case 'i':
            gIterations = atol( optarg );
            break;

         case 'N':
            gNSamples = atol( optarg );
            break;

         case 'B':
            gActionType = eBeamform;
            break;

         case 'T':
            gActionType = eBeamformTest;
            break;

         case 'd' :
            gActionType = eDumpData;
            break;

         case 'p' :
            printf("TEST optarg = %s\n",optarg);fflush(stdout);
            gPol = atol( optarg );
            break;

         case 'f' :
            gDoPhaseCorrection = 1;
            break;

         case 'c' :
            gDoCalibrate = atol( optarg );
            break;

         case 'o' :
            gOutFileName = optarg;
            break;

         case 't' :
            gTestPhaseDeg = atof(optarg);
            break;

         case 'n' :
            gPhaseNormalisation = (ePhaseNormalisation_t)atol(optarg);
            break;

         case 'x' :
            gMaxTimeSteps = atol(optarg);
            break;

         case 'r' :
            gUseRandomPhase = atol(optarg);
            break;

         case 's' :
            gOptimisePhaseStep = atof(optarg);
            break;

         case 'R' :
            gRefAnt = atol( optarg );
            break;

         case 'D' :
            gMeanDelaysString = optarg;
            break;

         case 'P' :
            gPhaseOffsetsString = optarg;
            break;

         case 'z' :
            gOptimiseRadius = atof( optarg );
            break;

         case 'X' :
            gCalPhaseOffsetsString = optarg;
            break;

         case 'S' :
            gSign = atof( optarg );
            break;

         case 'C' :
            gCalSign = atof( optarg );
            break;

         case 'L' :
            gExecuteLoop = 1;
            break;

         case 'F' :
            gSlopeFitFile = optarg;
            printf("optarg = %s\n",optarg);
            break;

         case 'l' :
            if ( optarg[0] != '-' ){
                gAntennaListStr = optarg;
                ParseCommaList( optarg, gAntennaListToProcess, "," );
            }else{
                gAntennaListStr = "";
                gAntennaListToProcess.clear();
            }
            break;

         case 'g' :
            if( strcasecmp(optarg,"quick_cal")==0 || strcasecmp(optarg,"quickcal")==0 ){
               gActionType = eQuickCal;
               printf("ACTION = QUICK_CAL\n");
            }else{
               if( strcasecmp(optarg,"optimise_quick_cal")==0 || strcasecmp(optarg,"optimise_quickcal")==0 ){
                  gActionType = eOptimiseQuickCal;
                  printf("ACTION = OPTIMISE QUICK_CAL\n");
               }else{
                  printf("ERROR : UNKNOWN ACTION REQUIRED !!!\n");
               }
            }
            break;
            
         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
   
   if( strlen(gMeanDelaysString.c_str()) > 0 && strlen(gPhaseOffsetsString.c_str()) > 0 ){
       printf("ERROR : -D and -P cannot be used in the same time !\n");
       exit(-1);
   }
   
   if( strlen(gMeanDelaysString.c_str()) > 0 ){
      MyParser pars = gMeanDelaysString.c_str();
      CMyStrTable delays_strings;
      pars.GetItems( delays_strings );
      
      gMeanDelays.clear();
      for(int i=0;i<delays_strings.size();i++){
         gMeanDelays.push_back( atof(delays_strings[i].c_str()) );
      }
      
      if( gRefAnt >= 0 ){
          if( gRefAnt < gMeanDelays.size() ){
              for(int i=0;i<gMeanDelays.size();i++){
                  gMeanDelays[i] = gMeanDelays[i] - gMeanDelays[gRefAnt];
              }
          }else{
              printf("Reference antenna %d outside allowed range (max = %d)\n",gRefAnt,int(gMeanDelays.size()));
          }
      }
   }

   if( strlen(gPhaseOffsetsString.c_str()) > 0 ){
      MyParser pars = gPhaseOffsetsString.c_str();
      CMyStrTable strings;
      pars.GetItems( strings );
      printf("Parsed %d phase delays\n",int(strings.size()));
      
      gPhaseOffsets.clear();
      for(int i=0;i<strings.size();i++){
         gPhaseOffsets.push_back( atof(strings[i].c_str()) );
      }

      if( gRefAnt >= 0 ){
          if( gRefAnt < gPhaseOffsets.size() ){
              for(int i=0;i<gPhaseOffsets.size();i++){
                  gPhaseOffsets[i] = gPhaseOffsets[i] - gPhaseOffsets[gRefAnt];
              }
          }else{
              printf("Reference antenna %d outside allowed range (max = %d)\n",gRefAnt,int(gPhaseOffsets.size()));
          }
      }

   }
   
   if( strlen(gCalPhaseOffsetsString.c_str()) > 0 ){
       MyParser pars = gCalPhaseOffsetsString.c_str();
       CMyStrTable strings;
       pars.GetItems( strings );
       printf("Parsed %d phase delays\n",int(strings.size()));
       
       for(int i=0;i<strings.size();i++){
          gAntennaPhaseDiff[i][1] = atof( strings[i].c_str() );
       }

   }

   if( strlen(gSlopeFitFile.c_str()) > 0 ){
       int read = read_file( gSlopeFitFile.c_str(), gSlopeParameters, 0, 0, 1, 3 );
       printf("Read %d values from file %s\n",gSlopeFitFile.c_str());
/*       for(int i=0;i<gSlopeParameters.size();i++){
          cValue& val = gSlopeParameters[i];
          
          printf("\t%d %.4f %.4f\n",int(val.x),val.y,val.z);
       }      */
   }
} 

void print_parameters()
{
   mystring szTestString,szTestString2,szCalPhaseOffsets;
   for(int i=0;i<gMeanDelays.size();i++){
      char szTmp[64];
      sprintf(szTmp,"%.2f,",gMeanDelays[i]);
      
      szTestString += szTmp;
   }
   for(int i=0;i<gPhaseOffsets.size();i++){
      char szTmp[64];
      sprintf(szTmp,"%.2f,",gPhaseOffsets[i]);
      
      szTestString2 += szTmp;
   }

   for(int i=0;i<MAX_ANTS;i++){
      char szTmp[64];
      sprintf(szTmp,"%.2f,",gAntennaPhaseDiff[i][1]);
      
      szCalPhaseOffsets += szTmp;
   }

   printf("##############################################\n");
   printf("PARAMETERS:\n");
   printf("##############################################\n");
   printf("Antenna1 = %d\n",antenna1);
   printf("Antenna2 = %d\n",antenna2);      
   printf("Execute loop = %d\n",gExecuteLoop);
   printf("N samples for fine channaliser FFT = %d\n",gNSamples);
   printf("Polarisation = %d\n",gPol);
   printf("Phase normalisation = %d\n",gPhaseNormalisation);
   printf("Phase correction = %d\n",gDoPhaseCorrection);   
   printf("Apply calibration coefficients = %d\n",gDoCalibrate);
   printf("Action   = %s ( execute in loop = %d )\n",get_action_desc(gActionType),gExecuteLoop);
   printf("Outfile  = %s\n",gOutFileName.c_str());
   printf("Antenna list  = %s -> %s\n",gAntennaListStr.c_str(),ListToStr( gAntennaListToProcess ).c_str() );
   printf("Test angle = %.2f [deg]\n",gTestPhaseDeg);
   printf("Max time steps = %d\n",gMaxTimeSteps);
   printf("Random phase   = %d\n",gUseRandomPhase);
   printf("%d mean delays    = %s (%s)\n",int(gMeanDelays.size()),gMeanDelaysString.c_str(),szTestString.c_str());
   printf("%d phase offsets  = %s (%s)\n",int(gPhaseOffsets.size()),gPhaseOffsetsString.c_str(),szTestString2.c_str());
   printf("Calibration phases = %s\n",szCalPhaseOffsets.c_str());
   printf("Reference antenna = %d\n",gRefAnt);
   printf("Sign geometry     = %.2f\n",gSign);
   printf("Sign cal. sol.    = %.2f\n",gCalSign);
   printf("gSlopeFitFile     = %s\n",gSlopeFitFile.c_str());
   if( gSlopeParameters.size() > 0 ){
       for(int i=0;i<gSlopeParameters.size();i++){
          cValue& val = gSlopeParameters[i];

          printf("\t\t%d %.4f %.4f\n",int(val.x),val.y,val.z);
       }
   }
   printf("OPTIMISATIONS PARAMETERS:\n");
   printf("\tIterations        = %d\n",gIterations);
   printf("\tPhase step        = %.2f [degrees]\n",gOptimisePhaseStep);
   printf("\tRadius            = %.2f [degrees]\n",gOptimiseRadius);
   printf("##############################################\n");
   
   // 
   printf("Usage example [0][1] = %.4f , [15][1] = %.4f\n",gAntennaPhaseDiff[0][1],gAntennaPhaseDiff[15][1]);
}

/*time_t get_unixtime_from_local_string_localfunc( const char* szDTM )
{
   struct tm local_time_tm;
        memset( &local_time_tm, '\0', sizeof(struct tm));

        // temporary correction due to fact that strptime does not fill fields :
        // tm_isdst = 1, tm_gmtoff = -10800,  tm_zone = 0x85e85a8 "CLST"
        // thus not working exactly good ... , but this is now 
   // filling current values of this field , which may not work for past 
        // and future dates ...
        time_t ut_time=get_dttm();
        localtime_r( &ut_time , &local_time_tm );

   strptime( szDTM, "%Y%m%d_%H%M%S", &local_time_tm );
        time_t ret = mktime( &local_time_tm );
        return ret;             
}


double hdf5filename2uxtime( const char* filename, const char* format="channel_cont_0_%d_%d_0.hdf5" )
{
   int dtm = 0;
   int local_seconds = 0;
   
   if( sscanf( filename, format, &dtm, &local_seconds ) != 2 ){
      printf("ERROR : could not parse date and time from file name %s\n",filename);
      return -1;
   }
   
   char szMidnight[128];
   sprintf(szMidnight,"%d_000000", dtm );
   time_t uxtime = get_unixtime_from_local_string_localfunc( szMidnight );
   
   printf("TEST uxtime = %d, local_seconds = %d\n",uxtime,local_seconds);
   
   return (uxtime + local_seconds);
}*/

int main(int argc,char* argv[])
{
   std::string filename = "channel_cont_0_20190307_25356_0.hdf5";
   if ( argc>=2 ){
       filename = argv[1];
   }
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
   
   printf("Action = dump data\n");
   dump_data( data, n_ants, n_pols, gOutFileName.c_str(), gPol );   
}


