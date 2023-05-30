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

// local includes:
#include "calsolutions.h"
#include "mwaconfig.h"
#include "defines.h"

#define MAX_ANTS 50
// current coefficients are for EDA2 and 48 antennas as for 20190610 :
// #define gAntennaPhaseDiff gEDA2_AntennaPhaseDiff_MSok_Norm_m180_180_20190609
double gBeamformingCoeff[MAX_ANTS][2];
#define gAntennaPhaseDiff gBeamformingCoeff
int do_phase_correction=0;
int do_calibrate=0;
int calibrate_sign=1;

typedef struct {
    char re;  
    char im;  
} complex_t;

string gAntennaLocationsFile;
MWAConfig gConfig;

string gInputFilename;

enum eActionType_t { eCorrelateAntennas=1, eBeamform=2, eDumpData=3, eQuickCal=4, eOptimiseQuickCal=5, eBeamformTest=6 };

enum ePhaseNormalisation_t { eNoPhaseNorm=0, ePhaseNorm0_360=1, ePhaseNorm_m180_180=2 };
ePhaseNormalisation_t gPhaseNormalisation = ePhaseNorm0_360;

// global name for input hdf5 file :
string gInputHdf5Filename;


int antenna1 = 0;
int antenna2 = 1;
int gDebugLevel = 0;
int gPol=0;
string gPolName;
string gOutFileName;
vector<int> gAntennaListToProcess;
string      gAntennaListStr;
int         gNumberOfAntennasToProcess=-1; // this is to test beamforming and how power increases as I am adding more antennas !

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


// list of flagged antennas :
std::vector<int> gFlaggedAntennaList;

// gain amplitude:
double gGainAmplitude = 1.00;
CValueVector gGainAmpVsUxtime;
string gGainAmpVsUxtime_FileName;

// Cal.Sol. phase vs. unixtime 
string gCalSolPhase_vs_unixtime_FileName;
CalSolValues gCalSolPhase_vs_unixtime;

bool gSteveOrderRemoveAutos = false;

// using fitted cal. solutions :
bool gUseFittedCalSolutionAmplitude = false;
CCalSolFits gFittedCalSolutions; // amplitudes 
bool gInverseCalSolAmplitude = false;

CCalSolFits gFittedCalSolutionsPhase; // fitted phases 
double gFitUxTimeStart = 1568352407.00;
double gFitUxTimeEnd   = 1568352407.00 + 33879; // fits as in 20190913_eda2_channeliseddata_ch204_24h_Lightcurve_vs_calsol_FitPhasePoly.odt and 20190913_eda2_channeliseddata_ch204_24h_Lightcurve_vs_calsol_UPDATE.odt

string gPointingString;
double gPointingAz_DEG = -1000;
double gPointingElev_DEG = -1000;


// saving data for Pulsar timing :
string gOutputTimeSeriesFileBin;
string gOutputTimeSeriesFileText;
string gOutputTimeSeriesFileBase;
bool   gSaveTxtTimeSeries = false;

// saving statistics of zeros for hdf5 file :
string gZeroStatFile;

// save power of each antenna vs. time :
string gBaseNamePowerPerAnt; // empty - means don't do this 

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
   for(int i=0;i<int(list.size());i++){
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
      while(ptr[0]==' ' && i < int( strlen(szToken)) ){
         ptr++;
         i++;
      }   


      out_list.push_back(atol(ptr));
   }
      
   return int(out_list.size());
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


// see  : 20190913_eda2_channeliseddata_ch204_24h_Lightcurve_vs_calsol_UPDATE.odt 
// /home/msok/Desktop/EDA2/data/2019_09_13_24hours_ch204/calsolutions/24hours/SunCal/only_continous/amplitude/mean_gain_vs_time/
double calsolution_amplitude( double uxtime, int antenna_index, bool bInverse=false )
{
   if( gUseFittedCalSolutionAmplitude ){
      double x = (uxtime - gFitUxTimeStart); // WARNING hardcoded value of Time0 (subtracted from all unix times) 
   
      /*if( x < 0 ){
         return 1.00;
      }*/
   
      if( int(gFittedCalSolutions.size()) > 0 ){
         // if fit per antenna provided use it here :
         double ret = gFittedCalSolutions.value( x, antenna_index );
         if( gDebugLevel >= 3 ){
            printf("DEBUG : fitted solutions per antenna returned gain = %.4f for antenna = %d at x = %.4f\n",ret,antenna_index,x);
         }
      
         if( bInverse ){
            return (1.00/ret);
         }else{
            return ret;
         }
      }
   
/*    OLD hardcoded fit for 20190913 (X pol.)  
      if( x>=0 && x<=34000 ){
          double ret = +0.92261057056020945310592651367187500000000000000000 * pow(x,0) +0.00004197645466774702072143554687500000000000000000 * pow(x,1) +0.00000000392538291541649764226917795895133167505264 * pow(x,2) -0.00000000000111845686519874366253272934823570494700 * pow(x,3) +0.00000000000000007145466544168238824806903947660375 * pow(x,4) -0.00000000000000000000184691940080214429002934791552 * pow(x,5) +0.00000000000000000000000001351579410076390896011092 * pow(x,6) +0.00000000000000000000000000000019691228579901170100 * pow(x,7) -0.00000000000000000000000000000000000271086534647991 * pow(x,8);

          if( bInverse ){
             return (1.00/ret);
          }else{
             return ret;
          }  
      }*/
   }
   
   return 1.000;
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
      printf("Written %.0f bytes to file %s, should write %d x %d = %.1f bytes\n",written_bytes, szOutFileName, int(data.size()), int(sizeof(data[0])), expected_write );
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
            
            
      int time_count = int(data.size()) / n_inputs;            
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
    while( start_index < int( data.size() ) ){
        std::vector< complex_t > samples;
        
        // add n_samples for fine channalisation :
        for(int i=0;i<n_samples;i++){
            int idx = start_index + i*ants_pols;
            
            if( idx < int( data.size() ) ){
                samples.push_back( data[ idx ] );
            }
        }    

                
        if( int(samples.size()) == n_samples ){
            // add samples to output list only if matching the required number :
            out_ant_data.push_back( samples );
        }            
        
        start_index += n_samples*ants_pols;
    }
    
    return int(out_ant_data.size());
}


int get_ant_data(std::vector< complex_t >& data, int ant_idx, int n_ants, int n_pols, std::vector< complex_t >& out_ant_data, int pol=0, int max_timesteps=-1 )
{
    out_ant_data.clear();
    
    int ants_pols = n_ants * n_pols;    
    int idx = n_pols*ant_idx + pol; // start at given antenna index (only X polarisation implemented yet )
    printf("get_ant_data : start_index = %d\n",idx);
    while( idx < int(data.size()) ){
        out_ant_data.push_back( data[ idx ] );
        
        idx += ants_pols;
        
        if( max_timesteps > 0 ){
           if( int(out_ant_data.size()) >= max_timesteps ){
              break;
           }
        }
    }
    
    return int(out_ant_data.size());
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
    for(int i=0;i<int(ant_list.size());i++){
        if( ant == ant_list[i] ){
            return i;
        }
    }

    return -1;
}

std::complex<double> phase_factor_func( double freq_mhz, double delay_picosec )
{
//   double delay_sec = delay_picosec * 1e-12;
//   double freq_hz   = freq_mhz * 1e6;
//   double oversampling = (32.00/27.00);
   
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
//   int do_phase_correction = gDoPhaseCorrection;
//   int do_calibrate        = gDoCalibrate;
//   double calibrate_sign   = gCalSign;

   if( p_phase_offsets ){
      mystring szTestOffsets;
      for(int i=0;i<int(p_phase_offsets->size());i++){
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
        if( int(gAntennaListToProcess.size()) > 0 ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                if( gDebugLevel > 1 ){printf("WARNING : antenna %d not in the list -> skipped\n",ant);}
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
            if( int(data_per_ant->size()) == block_count ){
                printf("\tNumber of blocks = %d -> ok (same as all the others)\n",(int)(data_per_ant->size()));
            }else{
                printf("\tERROR : number of blocks = %d != previous ( = %d )\n",(int)(data_per_ant->size()),block_count);
                exit(-1);
            }
        }else{
            block_count = int(data_per_ant->size());
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
            if( int(gAntennaListToProcess.size()) > 0 ){
                if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                    if( gDebugLevel > 1 ){ printf("WARNING : antenna %d not in the list -> skipped\n",ant); } 
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
                   if( int(gSlopeParameters.size()) > 0 ){
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
        // int middle_ch = int(beamformed_data.size())/2;
        for(int ch=0;ch<int(beamformed_data.size());ch++){
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
    for(int ch=0;ch<int(beamformed_data_accum.size());ch++){
//         double mag = std::norm( beamformed_data_accum[ch] ); // was abs
//        double mag = beamformed_data_accum[ch].real()*beamformed_data_accum[ch].real() + beamformed_data_accum[ch].imag()*beamformed_data_accum[ch].imag();        
//        double phase_deg = std::arg(  beamformed_data_accum[ch] );        
        total_power += mean_spectrum[ch];

        if( ch == 16 ){    
//           fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,mag,phase_deg,mean_spectrum[ch]/beamformed_count,gFileUxTime);
//           fprintf(beamformed_f,"%d %.4f %.4f %.4f\n",ch,mag,phase_deg,mean_spectrum[ch]/beamformed_count);
        }
    }
    // test with total power over all channels
    printf("DEBUG : beamformed_data_accum.size() = %d\n",int(beamformed_data_accum.size()));     
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


double avg_power( std::vector< complex_t >& ant_data )
{
   double sum = 0;
   if( int(ant_data.size()) <= 0 ){
      return 0.00;
   }
   
   for(int i=0;i<int(ant_data.size());i++){
      double real = ant_data[i].re;
      double imag = ant_data[i].im;
      
      sum += (real*real + imag*imag);
   }
   
   return sum / int(ant_data.size());
}

// _antennae
bool calc_geometric_pointing_delays( vector<double>& geometric_delays, double freq_mhz )

{
   if( strlen(gPointingString.c_str())>0 && gPointingAz_DEG > -360 && gPointingElev_DEG > -360 && gConfig.NAntennae() > 0 ){
      double freq_hz = freq_mhz*1e6;
      printf("Calculating pointing delays for pointing direction (AZ,ELEV) = (%.4f,%.4f) [deg] at %.2f [Hz]\n",gPointingAz_DEG,gPointingElev_DEG,freq_hz);
      
      double za_deg = 90.00 - gPointingElev_DEG;
      double za_rad = za_deg*(M_PI/180.00);
//      double az_rad = gPointingAz_DEG*(M_PI/180.00);
      double phi_rad = (90.00 - gPointingAz_DEG)*(M_PI/180.00); // phi is from X axis which is in antenna locations coordinate system East-West (azimuth is in North-East !)
      
      double pointing_vector[3];
      pointing_vector[0] = sin( za_rad )*cos( phi_rad );
      pointing_vector[1] = sin( za_rad )*sin( phi_rad );
      pointing_vector[2] = cos( za_rad );
      
      
      geometric_delays.clear();
      for(int ant=0;ant<int(gConfig.NAntennae());ant++){
         const MWAAntenna& ant_info = gConfig.Antenna( ant );
         
         // ant_info.position[0];
         // dot product of baseline_vector * antenna_position 
         double dot_product = pointing_vector[0]*ant_info.position[0] + pointing_vector[1]*ant_info.position[1] + pointing_vector[2]*ant_info.position[2];
         double delay_sec = dot_product / SPEED_OF_LIGHT;
         double delay_ns  = delay_sec*1e9;
         double phase_rad = 2*M_PI*freq_hz*delay_sec; // + SIGN FOR SURE !!! WARNING : - will make it work with azimuth from South towards West, but
                                                      //   + makes it proper azimuth from North towards East !!!
         double phase_deg = phase_rad*(180.00/M_PI);
         
         geometric_delays.push_back( phase_deg );         
         printf("%s : (%.3f,%.3f,%.3f) [m] -> delay = %.2f [ns] , phase = %.2f [deg]\n",ant_info.name.c_str(),ant_info.position[0],ant_info.position[1],ant_info.position[2],delay_ns,phase_deg);
      }
      
      return true;
   }

   return false;
}

double beamform2( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, const char* szOutFileName=NULL, int pol=0, 
                vector< std::complex<double> > * p_beamformed_data_out=NULL, 
                std::vector<double>* p_phase_offsets=NULL,
                int debug_level=1, 
                int do_average_mean = 1 )
{  
//   do_average_mean = 0;
   printf("beamform2 : do_average_mean = %d\n",do_average_mean);
   string szPol = "XX";
   if ( pol == 1 ){
      szPol = "YY";
   }
   
   if( strlen( gPolName.c_str() ) ){
      printf("INFO : external polarisation provided = |%s|",gPolName.c_str());
      szPol = gPolName.c_str();
   }else{
      printf("WARNING : no polarisation name provided -> this may cause problems for stations like EDA2 which has polarisations swapped\n");
      sleep(300);
   }
   
   
   // check if pointing direction different than zenith is specified :
   vector<double> geometric_delays;
   if( strlen(gPointingString.c_str())>0 && gPointingAz_DEG > -360 && gPointingElev_DEG > -360 && gConfig.NAntennae() > 0 ){
      // bool calc_geometric_pointing_delays( vector<double>& geometric_delays, double freq_mhz )
      if( !calc_geometric_pointing_delays( geometric_delays, gFreqMHz ) ){
         printf("ERROR : could not calculate geometric delays for pointing string = %s -> (az,elev) = (%.4f,%.4f) [deg]\n",gPointingString.c_str(),gPointingAz_DEG,gPointingElev_DEG);
         exit(-1);
      }
   }

   double gain_amplitude = gGainAmplitude;
   if( int(gGainAmpVsUxtime.size()) > 0 && gFileUxTime>0 ){
      // MEAN gain vs. time , mean for all antennas :
      cValue* pGain = gGainAmpVsUxtime.find_closest_value( gFileUxTime );
      
      if( pGain ){
         gain_amplitude = pGain->y;
         printf("INFO : using gain = %.4f for uxtime = %d\n",gain_amplitude,gFileUxTime);
      }else{
         printf("WARNING : could not find gain amplitude for unixtime = %d\n",gFileUxTime);
      }
/*   }else{
      if( gUseFittedCalSolutionAmplitude ){
         gain_amplitude = calsolution_amplitude( gFileUxTime );
      }*/
   }

   // to read all samples of given antenna 
   int time_steps = int( data.size() ) / (n_ants*n_pols);
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
      for(int i=0;i<int(p_phase_offsets->size());i++){
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
   
   // try to use fitted phase :
   if( int(gFittedCalSolutionsPhase.size()) > 0 && gFileUxTime>0 ){
       // double ret = gFittedCalSolutions.value( x, antenna_index );
       double x = (gFileUxTime - gFitUxTimeStart);
       if( x>=0 && x<=gFitUxTimeEnd ){
          printf("x = %.2f within fit range (0,%.1f) -> using fitted values ( flagged antennas = %d )\n",x,gFitUxTimeEnd,int(gFlaggedAntennaList.size()));
          for(int ant=0;ant<n_ants;ant++){
             double fitted_value = gFittedCalSolutionsPhase.value( x, ant );
             printf("\tPhase( %.1f ) = %.2f (was %.2f)\n",x,fitted_value,gAntennaPhaseDiff[ant][1]);
             gAntennaPhaseDiff[ant][1] = fitted_value;
          }
       }
   }
   
   if( !szOutFileName || strlen(szOutFileName) == 0 ){
      szOutFileName = "beamformed2.txt";
   }

//   vector< std::vector< complex_t > > antenna_data;  // MAX_ANTS does not work here ???
//   antenna_data.reserve( n_ants );
   std::vector< complex_t > antenna_data[MAX_ANTS]; // 48 works ok MAX_ANTS not !
   
   std::complex<double> beamformed_data_accum; // beamformed data - coherently added 16 EDA tiles 
   double mean_spectrum = 0.00;    
//   int beamformed_count = 0;
    
   int ants_to_process = n_ants; // mainly for debugging when I wanted to only check single antenna : 
   if ( gNumberOfAntennasToProcess > 0 ){
      printf("Forcing only %d first antennas to be processed\n",gNumberOfAntennasToProcess);
      ants_to_process = gNumberOfAntennasToProcess;
   }
//   int block_count = -1;
   int n_ant_processed = 0;
   
   cCalSolsVsUxtime* pCalSolPhases = NULL;
   if( int(gCalSolPhase_vs_unixtime.size()) ){
      pCalSolPhases = gCalSolPhase_vs_unixtime.find_closest_value( gFileUxTime );
      
      if( pCalSolPhases ){
         printf("INFO : using calsol. phases for unixtime = %d\n",gFileUxTime);
      }else{
         printf("WARNING : could not find cal.sol. phases for unixtime = %d\n",gFileUxTime);
      }
   }
   
   for(int ant=0;ant<ants_to_process;ant++){ // loop over AAVS antenna 
        if( int(gAntennaListToProcess.size()) ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                if( gDebugLevel > 1 ){ printf("WARNING : antenna %d not in the list -> skipped\n",ant); }
                continue;
            }
        }
        if( int(gFlaggedAntennaList.size()) > 0 ){
            if( is_in_list( gFlaggedAntennaList , ant ) >= 0 ){
                if( gDebugLevel > 1 ){ printf("WARNING : antenna %d is flagged -> skipped\n",ant); }
                continue;
            }
        }
   
        printf("Selecting data for antenna = %d\n",ant);fflush(stdout);
        int n_blocks = get_ant_data( data, ant, n_ants, n_pols, antenna_data[ant], pol, gMaxTimeSteps );
        printf("Selected %d blocks of %d samples for antenna %d (total %d samples), selected number of samples = %d ( data[time=0] = %d / %d )\n",n_blocks,time_steps,ant,(n_blocks*time_steps),(int)(antenna_data[ant].size()),antenna_data[ant][0].re,antenna_data[ant][0].im);
        
        if( time_steps != int(antenna_data[ant].size()) ){
            printf("ERROR in code : expected %d time_steps, but read %d -> cannot continue\n",time_steps,(int)(antenna_data[ant].size()));
            exit(-1);
        }        
        std::vector< complex_t >& ant_data = antenna_data[ant];
        
        if( strlen(gBaseNamePowerPerAnt.c_str()) > 0 ){
           // if set save power from each antenna to file :
           char szAntFileName[256];
           sprintf(szAntFileName,"%s%03d_%s.txt",gBaseNamePowerPerAnt.c_str(),ant,szPol.c_str());
           FILE* out_ant_f = fopen( szAntFileName , "a+" );
           double ant_power = avg_power( ant_data );
           fprintf(out_ant_f,"%d %.8f\n",(int)gFileUxTime,ant_power);
           fclose(out_ant_f);
        }
        
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
           
           
           // TODO : add using closest calibration solution phase in time !
           if( pCalSolPhases ){
              if( ant < int(pCalSolPhases->antenna_cal_solutions.size()) ){
                 delta_phase_rad = (pCalSolPhases->antenna_cal_solutions[ant])*(M_PI/180.00);
              }else{
                 printf("ERROR : number of calibration solutions = %d is smaller than antenna index = %d\n",int(pCalSolPhases->antenna_cal_solutions.size()),ant);
                 exit(-1);
              }
           }

           std::complex<float> phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) ); // should be opposite sign to what I got from "calibration" as I cross-correlated 
                                                                                                                 // ANT_0 and ANT_N which means multiplied Vis_0 * Conj(Vis_N) and took phase of this 
                                                                                                                 // so now to "undo" the phase I need to multiply by exp(-i * phi)
//           if( ant_data.size() ){
//               printf("DEBUG2 : ant = %d, t = %d / %d -> vis = %d / %d\n",ant,0, ant_data.size(), ant_data[0].re , ant_data[0].im );
//           }

           if( int(geometric_delays.size()) > 0 && ant < int(geometric_delays.size()) ){
              // this means that pointing direction != ZENITH is specified :
              double pointing_direction_phase_deg = geometric_delays[ant];
              double pointing_direction_phase_rad = pointing_direction_phase_deg * (M_PI/180.00);
              
              std::complex<float> phase_pointing_factor = std::complex<float>( cos(pointing_direction_phase_rad), sin(pointing_direction_phase_rad) );  // - sign is already in there see : calc_geometric_pointing_delays
              phase_factor = phase_factor * phase_pointing_factor;
           }
           

           for(int t=0;t<int(ant_data.size());t++){
              std::complex<float> tmp_complex = std::complex<float>( ant_data[t].re , ant_data[t].im );
              tmp_complex *= phase_factor;
              
              if( int(gMeanDelays.size()) > 0 ){
                  if( ant < int(gMeanDelays.size()) ){
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
                 if( int(gPhaseOffsets.size()) > 0 ){
                     if( ant < int(gPhaseOffsets.size()) ){
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
    
    vector< double > autos;
    autos.assign( time_steps, 0.00 ); // to accumulate auto-correlations see page 8 (equations 41 - 46) in Steve Ord 2019 PASA paper 
    
    for(int time_step=0;time_step<time_steps;time_step++){
        if( debug_level > 0 ){
           if( (time_step % 1000) == 0 ){
              printf("Time step %d / %d ...\n",time_step,time_steps);fflush(stdout);
           }
        }

        for(int ant=0;ant<ants_to_process;ant++){
            if( int(gAntennaListToProcess.size())>0 ){
                if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                    // printf("WARNING : antenna %d not in the list -> skipped\n",ant);
                    continue;
                }
            }

            if( int(gFlaggedAntennaList.size()) > 0 ){
                if( is_in_list( gFlaggedAntennaList , ant ) >= 0 ){
                   if( gDebugLevel > 1 ){ printf("WARNING : antenna %d is flagged -> skipped\n",ant); }
                   continue;
                }
            }

        
            std::vector< complex_t >& ant_data = antenna_data[ant];
            
            // add to beamformed data 
            std::complex<float> tmp( ant_data[time_step].re , ant_data[time_step].im );
            if( gUseFittedCalSolutionAmplitude ){
               gain_amplitude = calsolution_amplitude( gFileUxTime, ant, gInverseCalSolAmplitude  );
            }
            tmp = tmp * std::complex<float>( gain_amplitude , 0.00 );
            beamformed_data[time_step] += tmp;            

            // accumulate auto-correlations see page 8 (equations 41 - 46) in Steve Ord 2019 PASA paper             
            autos[time_step] += std::norm( tmp ); // https://en.cppreference.com/w/cpp/numeric/complex/norm = absolute square 
        }
        
    }
    
    FILE* out_timeseries_f = NULL;
    FILE* out_timeseries_text_f = NULL;
    if( strlen(gOutputTimeSeriesFileBin.c_str()) ){
       out_timeseries_f = fopen( gOutputTimeSeriesFileBin.c_str() , "wb" );
    }
    if( strlen(gOutputTimeSeriesFileText.c_str()) && gSaveTxtTimeSeries ){
       out_timeseries_text_f = fopen( gOutputTimeSeriesFileText.c_str() , "w" );
    }

    for(int time_step=0;time_step<int(beamformed_data.size());time_step++){
/*       if( do_average_mean > 0 ){ // && 0 ){ - testing for Steve Ord's version - gives values << 0 !!!???
          beamformed_data[time_step] = std::complex<double>( beamformed_data[time_step].real() / n_ants , beamformed_data[time_step].imag() / n_ants ); 
          autos[time_step] = autos[time_step] / n_ants;
       }*/
//       beamformed_data[time_step] = std::complex<double>( beamformed_data[time_step].real() , beamformed_data[time_step].imag() ); 
             
       beamformed_data_accum += beamformed_data[time_step];

       // 20190916 - testing subraction of auto-correlations as in the Steve Ord paper to get same as correlator output !             
       double power = std::norm( beamformed_data[time_step] ); // Steve Ord:  - (autos[time_step]); // Must be power !!! std::abs( beamformed_data[time_step] ); 1) Returns the squared magnitude of the complex number z.       
       if( do_average_mean > 0 ){
          power = power / n_ants;
       }
       
       if( gSteveOrderRemoveAutos ){
          if( do_average_mean > 0 ){
             power = power - autos[time_step] / n_ants; // autos[time_step] is already std:norm (see lines above)
          }else{
             power = power - autos[time_step];
          }          
       }
       
       if( out_timeseries_text_f ){
          fprintf( out_timeseries_text_f , "%.8f %.8f\n",double(gFileUxTime+(time_step*(1.08/1000000.0))),power);
       }
       
// TEST        double power = std::abs( beamformed_data[time_step] );
       mean_spectrum += power;
    }     

    // mean power :    
// 2020-06-12 - this line was responsible for very low numbers in the resulting power - this is not acceptable plots looked like Quantised ...    
//    mean_spectrum = mean_spectrum / beamformed_data.size();
    
    char szOutFile[1024];
    sprintf(szOutFile,"meanpower_vs_time_pol%d.txt",gPol);
    FILE* out_mean_power_f = fopen(szOutFile,"a+");
    fprintf( out_mean_power_f , "%.8f %.8f %.8f %.8f %d\n",gPointingAz_DEG,gPointingElev_DEG,double(gFileUxTime+((double(beamformed_data.size())/2.00)*(1.08/1000000.0))),(mean_spectrum/int(beamformed_data.size())),int(beamformed_data.size()));
    fclose( out_mean_power_f );
    
    if( out_timeseries_f ){
       fwrite( &(beamformed_data[0]), beamformed_data.size(), sizeof( std::complex<double> ), out_timeseries_f );
       fclose( out_timeseries_f );
    }
    if( out_timeseries_text_f ){
       fclose( out_timeseries_text_f );
    }

    
    FILE* beamformed_f = fopen( szOutFileName,"a+");
    double mag = std::norm( beamformed_data_accum ) / int(beamformed_data.size()); // was abs
//        double mag = beamformed_data_accum[ch].real()*beamformed_data_accum[ch].real() + beamformed_data_accum[ch].imag()*beamformed_data_accum[ch].imag();        
    double phase_deg = std::arg(  beamformed_data_accum );        

    double final_power = mean_spectrum/int(beamformed_data.size());    
//    double final_power = mean_spectrum;
    if( gNumberOfAntennasToProcess > 0 ){
         fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,mag,phase_deg,final_power,gFileUxTime,gNumberOfAntennasToProcess);
    }else{
        fprintf(beamformed_f,"%s %s %.4fdeg %.4f %.4f %.4f %d\n",szInFileName,ListToStr(gAntennaListToProcess).c_str(),gTestPhaseDeg,mag,phase_deg,final_power,gFileUxTime);
    }
    fclose( beamformed_f );
    
    printf("Saved results to %s file , processed %d antennas\n",szOutFileName,n_ant_processed);
    
    if( p_beamformed_data_out ){
       p_beamformed_data_out->clear();
       
       for(int time_step=0;time_step<int(beamformed_data.size());time_step++){
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
   for(int block=0;(block<int(data_per_ant1.size()) && (block<gMaxTimeSteps || gMaxTimeSteps<0));block++){
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
           
//           double re = corr_reim[ch].real();
//           double im = corr_reim[ch].imag(); 
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
       mean_spectrum1[ch] = mean_spectrum1[ch] / mean_count1;
       fprintf(out_f,"%d %.4f\n",ch,mean_spectrum1[ch]);
   }
   fclose(out_f);     
   
   out_f = fopen( szOutFile2 , "w");
   for(int ch=0;ch<out_count1;ch++){ // was (n_samples/2)
       mean_spectrum2[ch] = mean_spectrum2[ch] / mean_count2;
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
   printf("beamform_era INPUT_FILE\n");
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "h";   
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);
   
      switch (opt) {
         case 'h' :
            usage();
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
   printf("##############################################\n");
   
}

int main(int argc,char* argv[])
{
   std::string filename = "data.volts";
   if ( argc>=2 ){
       filename = argv[1];
   }
   gInputFilename = filename.c_str();
   parse_cmdline( argc , argv );
   print_parameters();
   

/*   int n_samples = gNSamples; // number of samples to use for fine channalisation 
   int n_fine_channels = n_samples / 2;
   int ants_pols = 1; // dims_out[1];
   int n_ants = MAX_ANTS; // dims_out[1] / 2; // 32 / 2pols -> 16 antennas 
   int n_pols = 1;
   int samples_all = 0; // dims_out[0];
   int samples_per_pol = samples_all / ants_pols;
*/
   
}


void SavePhaseOffsets( std::vector<double> phase_offsets, int n_iter=-1, int pol=0 )
{
   char szFile[1024];   
   char szPol='X';
   if( pol > 0 ){
      szPol='Y';
   }   
   
   if( strlen( gPolName.c_str() ) ){
      printf("INFO : external polarisation provided = |%s|",gPolName.c_str());
      szPol = gPolName[0];
   }else{
      printf("WARNING : no polarisation name provided -> this may cause problems for stations like EDA2 which has polarisations swapped\n");
      sleep(300);
   }


   char szPrefixIter[64];
   strcpy(szPrefixIter,"initial");
   if( n_iter >= 0 ){
      sprintf(szPrefixIter,"iter%02d",n_iter);
   }

   sprintf(szFile,"%s_phase_vs_antenna_%c.txt",szPrefixIter,szPol);
   
   
   FILE* out_f = fopen(szFile,"w");
   mystring szPhaseOffsets;
   for(int i=0;i<int(phase_offsets.size());i++){
       printf("Antenna = %d phase offset = %.4f [deg]\n",i,phase_offsets[i]);
       fprintf(out_f,"%d %.4f\n",i,phase_offsets[i]);

       char szTmp[64];
       sprintf(szTmp,"%.2f,",phase_offsets[i]);
       szPhaseOffsets << szTmp;
   }
   printf("Iteration = %d : -X %s\n",n_iter,szPhaseOffsets.c_str());
   
   fclose(out_f);

}


void beamform_fit_phase_offsets( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, const char* szOutFileName, int pol, int n_iter )
{  
   // to read all samples of given antenna 
   int time_steps = int(data.size()) / (n_ants*n_pols);
   printf("Number of expected time steps = %d\n",time_steps);
   if( gMaxTimeSteps > 0 ){
       printf("WARNING : number of time steps overwritten with external %d value\n",gMaxTimeSteps);
       time_steps = gMaxTimeSteps;
   }
   
   
   
   if( !szOutFileName || strlen(szOutFileName) == 0 ){
      szOutFileName = "beamformed2.txt";
   }

//   vector< std::vector< complex_t > > antenna_data;  // MAX_ANTS does not work here ???
//   antenna_data.reserve( n_ants );
   std::vector< complex_t > antenna_data[MAX_ANTS]; // 48 works ok MAX_ANTS not !

// initialise phase offsets with just 0 degrees for antenna = 0    
   std::vector< double > phase_offsets;
   
//   std::complex<double> beamformed_data_accum; // beamformed data - coherently added 16 EDA tiles 
   // beamformed data using optimal delays (initialise with 0.00 for ant=0) :
   vector< std::complex<double> > beamformed_data_optimal;
   beamformed_data_optimal.assign( time_steps, 0.00 ); 
   
   
   double mean_spectrum = 0.00;    
//   int beamformed_count = 0;
    
   int ants_to_process = n_ants; // mainly for debugging when I wanted to only check single antenna : 
//   int block_count = -1;
   int n_ant_processed = 0;
   for(int ant=0;ant<ants_to_process;ant++){ // loop over antennas to find optimal phase shit for all of them :
        if( gAntennaListToProcess.size() ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                if( gDebugLevel > 1 ){ printf("WARNING : antenna %d not in the list -> skipped\n",ant); }
                continue;
            }
        }
   
        printf("Selecting data for antenna = %d\n",ant);fflush(stdout);
        int n_blocks = get_ant_data( data, ant, n_ants, n_pols, antenna_data[ant], pol, gMaxTimeSteps );
        printf("Selected %d blocks of %d samples for antenna %d (total %d samples), selected number of samples = %d\n",n_blocks,time_steps,antenna1,(n_blocks*time_steps),(int)(antenna_data[ant].size()));
        
        if( time_steps != int(antenna_data[ant].size()) ){
            printf("ERROR in code : expected %d time_steps, but read %d -> cannot continue\n",time_steps,(int)(antenna_data[ant].size()));
            exit(-1);
        }       
//        std::vector< complex_t >& ref_ant_data = antenna_data[0]; 
        std::vector< complex_t >& ant_data = antenna_data[ant];

        // std::complex<float> phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) );
        double phase_step = 1.00; // degree
        double phase_offset = -180.00; // degree;
        double max_mean_spectrum = -1e20;
        double max_phase_offset = 0.00;

                
        if( ant > 0 ){ // skip antenna 0 which is reference antenna                
           printf("Finding optimal phase for antenna = %d : phase range 0 - 360 in step of %.2f degrees\n",ant,phase_step);fflush(stdout);
           while( phase_offset <= 180.00 ){
               double phase_offset_rad = phase_offset * (M_PI/180.00);
               std::complex<float> phase_factor = std::complex<float>( cos(phase_offset_rad), sin(phase_offset_rad) );
        
               vector< std::complex<double> > beamformed_data; // beamformed data of all antennas (or EDA tiles) for given time step
               // initialise with already optimally beamformed data for ant-1 antennas :
               for(int i=0;i<int(beamformed_data_optimal.size());i++){
                  beamformed_data.push_back( beamformed_data_optimal[i] );
               }

               double mean_spectrum = 0.00;
               for(int time_step=0;time_step<time_steps;time_step++){
                   //if( (time_step % 1000) == 0 ){
                   //   printf("Time step %d / %d ...\n",time_step,time_steps);fflush(stdout);
                   //}
                   std::complex<float> tmp_complex = std::complex<float>( ant_data[time_step].re , ant_data[time_step].im );
                   tmp_complex *= phase_factor;
                                
                   // add to beamformed data 
                   beamformed_data[time_step] += std::complex<float>( tmp_complex.real() , tmp_complex.imag() );
                
                   double power = std::norm( beamformed_data[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] )
                   mean_spectrum += power;                
               }
            
               if( mean_spectrum > max_mean_spectrum ){
                   max_mean_spectrum = mean_spectrum;
                   max_phase_offset = phase_offset;
               }
               printf("\tAntenna = %d , test phase_offset = %.2f [deg], mean_spectrum = %.2f vs. max_mean_spectrum = %.2f at optimal_phase_offset = %.2f [deg]\n",ant,phase_offset,mean_spectrum,max_mean_spectrum,max_phase_offset);fflush(stdout);
             
               phase_offset += phase_step;
           }
        }
        
        
        printf("Antenna = %d : optimal phase offset = %.2f deg, mean_power = %.4f\n",ant,max_phase_offset,max_mean_spectrum);
        // found optimal phase now calculated beamformed sum for 2 antennas :
        phase_offsets.push_back( max_phase_offset );
        double max_phase_offset_rad = max_phase_offset*(M_PI/180.00);
        std::complex<float> phase_factor = std::complex<float>( cos(max_phase_offset_rad), sin(max_phase_offset_rad) );
        
        for(int time_step=0;time_step<time_steps;time_step++){
            std::complex<float> tmp_complex = std::complex<float>( ant_data[time_step].re , ant_data[time_step].im );
            tmp_complex *= phase_factor;
        
            beamformed_data_optimal[time_step] += std::complex<float>( tmp_complex.real() , tmp_complex.imag() );
        }
        
        
        n_ant_processed++;
    }
    
    mean_spectrum = 0.00;
    for(int time_step=0;time_step<time_steps;time_step++){
       double power = std::norm( beamformed_data_optimal[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] )
       mean_spectrum += power;
    }
    printf("Optimal power = %.4f , mean = %.4f\n",mean_spectrum,mean_spectrum/time_steps);

    SavePhaseOffsets( phase_offsets, -1, pol );    

    int iter=0;    
    printf("Improving with %d iterations ...\n",n_iter);
    while( iter < n_iter ){
       printf("Running iteration = %d\n",iter);
    
       beamform_optimise_phase_offsets( data, n_ants, n_pols, szInFileName, phase_offsets, szOutFileName, pol, iter );
       
       mean_spectrum = 0.00;
       for(int time_step=0;time_step<time_steps;time_step++){
          double power = std::norm( beamformed_data_optimal[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] )
          mean_spectrum += power;
       }
       printf("Optimal power = %.4f , mean = %.4f\n",mean_spectrum,mean_spectrum/time_steps);

       SavePhaseOffsets( phase_offsets, iter, pol );
       
       iter++;
    }
}

/*
void beamform2( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, const char* szOutFileName=NULL, int pol=0, 
                vector< std::complex<double> > * p_beamformed_data_out=NULL, 
                std::vector<double>* p_phase_offsets=NULL )

*/

void beamform_optimise_phase_offsets( std::vector< complex_t >& data, int n_ants, int n_pols, const char* szInFileName, std::vector< double >& phase_offsets, const char* szOutFileName, 
                                      int pol, int iteration, double phase_step )
{  
   // to read all samples of given antenna 
   int time_steps = data.size() / (n_ants*n_pols);
   printf("Number of expected time steps = %d\n",time_steps);
   if( gMaxTimeSteps > 0 ){
       printf("WARNING : number of time steps overwritten with external %d value\n",gMaxTimeSteps);
       time_steps = gMaxTimeSteps;
   }
   
   
   std::vector< complex_t > antenna_data[MAX_ANTS]; // 48 works ok MAX_ANTS not !
   std::vector< double > start_phase_offsets  = phase_offsets;
   std::vector< double > test_phase_offsets   = phase_offsets;
   
   // getting optimally beamformed data :
   vector< std::complex<double> > beamformed_data_optimal;
   int debug_level = 0;
   double start_mean_power = beamform2( data, n_ants, n_pols, szInFileName, szOutFileName, pol, &beamformed_data_optimal, &test_phase_offsets, debug_level, 0 /*do_average_mean*/ );
   printf("Inital value of mean power (before optimisation started) = %.4f\n",start_mean_power);

//   beamform2( data, n_ants, n_pols, szInFileName, szOutFileName, pol, &beamformed_data_optimal, &phase_offsets );
   
   
   double mean_spectrum = 0.00;    
//   int beamformed_count = 0;
    
   int ants_to_process = n_ants; // mainly for debugging when I wanted to only check single antenna : 
//   int block_count = -1;
   int n_ant_processed = 0;
   for(int ant=0;ant<ants_to_process;ant++){ // loop over antennas to find optimal phase shit for all of them :
        if( gAntennaListToProcess.size() ){
            if( is_in_list( gAntennaListToProcess , ant ) < 0 ){
                if( gDebugLevel > 1 ){ printf("WARNING : antenna %d not in the list -> skipped\n",ant); }
                continue;
            }
        }
   
        printf("Selecting data for antenna = %d\n",ant);fflush(stdout);
        int n_blocks = get_ant_data( data, ant, n_ants, n_pols, antenna_data[ant], pol, gMaxTimeSteps );
        printf("Selected %d blocks of %d samples for antenna %d (total %d samples), selected number of samples = %d\n",n_blocks,time_steps,antenna1,(n_blocks*time_steps),(int)(antenna_data[ant].size()));
        
        if( time_steps != int(antenna_data[ant].size()) ){
            printf("ERROR in code : expected %d time_steps, but read %d -> cannot continue\n",time_steps,(int)(antenna_data[ant].size()));
            exit(-1);
        }       
//        std::vector< complex_t >& ref_ant_data = antenna_data[0]; 
//        std::vector< complex_t >& ant_data = antenna_data[ant];

        // std::complex<float> phase_factor = std::complex<float>( cos(delta_phase_rad), sin(delta_phase_rad) );
        double phase_offset = -180.00; // degree;
        double phase_offset_end = +180.00;
        double max_mean_spectrum = start_mean_power;
        double max_phase_offset = start_phase_offsets[ant];
        int    max_found = 0;
        std::vector< double > tested_mean_powers;
 
        if( gOptimiseRadius > 0 ){
            phase_offset = start_phase_offsets[ant] - gOptimiseRadius;
            phase_offset_end = start_phase_offsets[ant] + gOptimiseRadius;
            
            printf("Antenna = %d : optimising in range [%.4f , %.4f] degrees in steps of %.4f degrees\n",ant,phase_offset,phase_offset_end,phase_step);
        }
                
        if( ant > 0 ){ // skip antenna 0 which is reference antenna                
           printf("Finding optimal phase for antenna = %d : phase range 0 - 360 in step of %.2f degrees\n",ant,phase_step);fflush(stdout);           
           while( phase_offset <= phase_offset_end ){
//               double phase_offset_rad = phase_offset * (M_PI/180.00);
//               std::complex<float> phase_factor = std::complex<float>( cos(phase_offset_rad), sin(phase_offset_rad) );
        
               vector< std::complex<double> > beamformed_data; // beamformed data of all antennas (or EDA tiles) for given time step
               test_phase_offsets[ant] = phase_offset;
               
               // beamform using test phase_offet for antenna = ant :
               beamform2( data, n_ants, n_pols, szInFileName, szOutFileName, pol, &beamformed_data, &test_phase_offsets, debug_level, 0 /* do_average_mean */ );
               

               double mean_spectrum = 0.00;
               for(int time_step=0;time_step<time_steps;time_step++){
                   double power = std::norm( beamformed_data[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] )
                   mean_spectrum += power;                
               }
               tested_mean_powers.push_back( mean_spectrum );
            
               if( mean_spectrum > max_mean_spectrum ){
                   max_mean_spectrum = mean_spectrum;
                   max_phase_offset  = phase_offset;
                   max_found         = 1;
               }
               printf("\tAntenna = %d , test phase_offset = %.2f [deg], mean_spectrum = %.2f vs. max_mean_spectrum = %.2f at optimal_phase_offset = %.2f [deg]\n",ant,phase_offset,mean_spectrum,max_mean_spectrum,max_phase_offset);fflush(stdout);
             
               phase_offset += phase_step;
           }

           if( max_found > 0 ){           
              phase_offsets[ant] =  max_phase_offset;
              test_phase_offsets[ant] = max_phase_offset; // setting test offset for the tested antenna to current offset 
           }else{
              test_phase_offsets[ant] = start_phase_offsets[ant];
              printf("WARNING : for antenna = %d could not find a phase offset improving maximum power ( max_mean_spectrum = %.4f ) restoring original value of offset = %.2f deg\n",ant,max_mean_spectrum,test_phase_offsets[ant]);
           }
        }else{
            max_mean_spectrum = start_mean_power;
        }
        
        
        if( max_found > 0 ){
           printf("Antenna = %d : optimal phase offset = %.2f deg, mean_power = %.4f (original phase offset = %.4f deg)\n",ant,max_phase_offset,max_mean_spectrum,start_phase_offsets[ant]);
        }else{
           printf("WARNING : Antenna = %d : could not find optimal phase offset improving mean_power = %.4f -> keeping original offset = %.4f deg (vs. ? = %.4f deg ? )\n",ant,max_mean_spectrum,phase_offsets[ant],start_phase_offsets[ant]);
           mystring szTestedMeanPowers;
           char szTmp[128];
           for(int i=0;i<int(tested_mean_powers.size());i++){
              sprintf(szTmp,"%.4f,",tested_mean_powers[i]);
              szTestedMeanPowers << szTmp;
           }
           printf("\t\tTested mean powers = %s\n",szTestedMeanPowers.c_str());
        }
        
        n_ant_processed++;
    }
    
    // calculate optimally beamformed data using final phase offsets :
    beamform2( data, n_ants, n_pols, szInFileName, szOutFileName, pol, &beamformed_data_optimal, &phase_offsets );
    
    mean_spectrum = 0.00;
    for(int time_step=0;time_step<time_steps;time_step++){
       double power = std::norm( beamformed_data_optimal[time_step] ); // Must be power !!! std::abs( beamformed_data[time_step] )
       mean_spectrum += power;
    }
    printf("Optimal power = %.4f , mean = %.4f\n",mean_spectrum,mean_spectrum/time_steps);

    mystring szPhaseOffsets;            
    for(int i=0;i<int(phase_offsets.size());i++){
       printf("Antenna = %d phase offset = %.4f [deg] (iteration = %d)\n",i,phase_offsets[i],iteration);
       
       char szTmp[64];
       sprintf(szTmp,"%.2f,",phase_offsets[i]);
       szPhaseOffsets << szTmp;
    }
    printf("-X %s\n",szPhaseOffsets.c_str());
}

