// See also :
// h5dump program - part of linux or needs instalation of some hdf5 package 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <complex>

// my :
#include <spectrometer.h>
#include <bg_date.h>

// 
#include <myparser.h>
#include <mystring.h>
#include <mystrtable.h>
#include <myfile.h>

#define MAX_ANTS 16

enum eActionType_t { eCorrelateAntennas=1, eBeamform=2 };

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
int gDoCalibrate=1;

// unix time :
int gFileUxTime = -1;

// pointing parameters :
mystring gMeanDelaysString;
vector<double> gMeanDelays;

// phase offsets instead of pointing delays from HydA calibration for example 
mystring gPhaseOffsetsString;
vector<double> gPhaseOffsets;

// reference antenna for delays or phase offsets :
int gRefAnt=-1;

// Sign to Test signs conventions :
double gSign = 1.00;
double gCalSign = 1.00;

// 
double gFreqMHz = (400.00/512.00)*204.00;

double gDelayAnt0=0.00;

string gSource;

typedef struct {
    char re;  
    char im;  
} complex_t;

eActionType_t gActionType = eCorrelateAntennas;

int rx2tile_mapping[MAX_ANTS] = { 0, 1 , 2 , 3 , 8 , 9 , 10 , 11 , 15 , 14 , 13 , 12 , 7 , 6 , 5 , 4 }; 
int aavsant2edaindex_mapping[MAX_ANTS] = { 0, 1, 2, 3, 15, 14, 13, 12, 4, 5, 6, 7, 11, 10, 9, 8 };


int tpmrx2edatileindex( int eda_tile_index )
{
//   int eda_tile_index = (rx-1);
   
   return rx2tile_mapping[eda_tile_index];   
}

int aavsant2edaindex( int aavs_index )
{
    return aavsant2edaindex_mapping[ aavs_index ];
}

double phase_factor_func( double freq_mhz, double delay_picosec )
{
   double delay_sec = delay_picosec * 1e-12;
   double freq_hz   = freq_mhz * 1e6;
   double oversampling = (32.00/27.00);
   
   // TEST :
//   freq_mhz = freq_mhz*oversampling;
   
   double phase_rad = gSign * 2.00*M_PI*freq_mhz*delay_picosec*(1e-6); // to avoid numerical errors 
// TEST :   phase_rad *= oversampling;

   // test with respect to antenna 0 :
//   double phase_rad_ant0 = 2.00*M_PI*freq_mhz*(12099.42)*(1e-6);
//   phase_rad -= phase_rad_ant0;   
   
   
   std::complex<double> phase_coeff( cos(phase_rad), sin(phase_rad) );
//   std::complex<double> phase_coeff = std::exp(1i * phase_rad);
   
   return phase_rad*(180.00/M_PI);
}

double delay_picosec_from_phase_deg( double freq_mhz, double phase_deg )
{
   double phase_rad = phase_deg*(M_PI/180.00);
   double delay_picosec = phase_rad / (gSign * 2.00*M_PI*freq_mhz)*1000000.00;
   
   return delay_picosec;
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


void usage()
{
   printf("delay2phase mean_delay_picosec_file.txt calibration_phase.txt -s SOURCE_NAME\n");
   printf("\t-S MULTIPLIER : mainly to test sign convention, but it will scale phase offsets by this number [default %.2f - SHOULD BE ONE !]\n",gSign);
   printf("\t-r reference_antenna [default %d]\n",gRefAnt);
   printf("\t-R value of delay to subtract [default %.4f]\n",gDelayAnt0);
   printf("\t-C CAL. MULTIPLIER : mainly to test sign convention for calibration phase [default %.2f - SHOULD BE ONE !]\n",gCalSign);
   printf("\t-s SOURCE_NAME     : name of source [default %s]\n",gSource.c_str());
   printf("\t-f FREQ_MHZ        : default = %.2f MHz\n",gFreqMHz);
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hS:r:C:s:R:f:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
           case 'f':
              gFreqMHz = atof( optarg );
              break;

           case 'r':
              gRefAnt = atol( optarg );
              break;

           case 'R':
              gDelayAnt0 = atof( optarg );
              break;

           case 's':
              gSource = optarg;
              break;

           case 'S' :
               gSign = atof( optarg );
               break;

           case 'C' :
               gCalSign = atof( optarg );
               break;


           case 'h' :
               usage();
               break;

         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }   
} 

int read_calfile( const char* calfile, vector<double>& cal_phase )
{
   cal_phase.clear();
   
   MyFile infile( calfile );   
   const char* pLine = NULL;

   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' ){
         continue;
      }

      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );

      int aavs_ant_index  = atol( items[0].c_str() );
      double phase_deg    = atof( items[1].c_str() );

      printf("AAVS-ANTENNA-%d %.4f [deg]\n",aavs_ant_index,phase_deg);
      cal_phase.push_back( phase_deg );
   }

   
}

void print_parameters()
{
   printf("##############################################\n");
   printf("PARAMETERS:\n");
   printf("##############################################\n");
   printf("Sign      = %.2f\n",gSign);
   printf("CAL. Sign = %.2f\n",gCalSign);
   printf("Reference antenna = %d\n",gRefAnt);
   printf("Reference mean delay = %.4f [ps]\n",gDelayAnt0);
   printf("Source    = %s\n",gSource.c_str());
   printf("Frequency = %.2f [MHz]\n",gFreqMHz);
   printf("##############################################\n");
   
   // 
}

int main(int argc,char* argv[])
{
   std::string filename = "mean_delay_vs_EDA-BF-INDEX.txt";
   if ( argc>=2 ){
       filename = argv[1];
   }
   std::string calfile = "phase_vs_AAVS-ANT-INDEX_GalTransit.txt";
   if ( argc>=3 ){
       calfile = argv[2];
   }
   std::string pointed_file = "phase_vs_AAVS-ANT-INDEX_HydATransit.txt";
   if ( argc>=4 ){
       pointed_file = argv[3];
   }

   
   parse_cmdline( argc , argv );
   print_parameters();
   
   vector<double> cal_phase;
   read_calfile( calfile.c_str(), cal_phase );
   
   vector<double> point_phase;
   read_calfile( pointed_file.c_str(), point_phase );
   
   MyFile infile( filename.c_str() );   
   const char* pLine = NULL;
   int ant_counter = 0;

   int good_counter = 0;   
   while( pLine = infile.GetLine( TRUE ) ){
      if( mystring::get_first_non_white( pLine )=='#' ){
         continue;
      }
                        
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
      
      int eda_bf_index  = atol( items[0].c_str() );
      char sign = items[0][0];
//      printf("sign = %c\n",sign);
      switch( sign ){
         case 'A' :
             eda_bf_index = 10;
             break;
         case 'B' :
             eda_bf_index = 11;
             break;
         case 'C' :
             eda_bf_index = 12;
             break;
         case 'D' :
             eda_bf_index = 13;
             break;
         case 'E' :
             eda_bf_index = 14;
             break;
         case 'F' :
             eda_bf_index = 15;
             break;
      
         default :
//             printf("ERROR : unknown number %s\n", items[0].c_str() );
             break;
      }
      const char* eda_bf_str = items[0].c_str();
      double mean_delay_picosec = atof( items[1].c_str() ); // - 7383.26; // wrt ANT=0
      if( ant_counter == 0 && fabs(gDelayAnt0)<0.0001 ){
         gDelayAnt0 = mean_delay_picosec;
         printf("SET reference delay to %.4f picoseconds\n",gDelayAnt0);
      }
      if( gRefAnt == 0 ){
          mean_delay_picosec = mean_delay_picosec - gDelayAnt0;
      }
      double phase_deg  = phase_factor_func( gFreqMHz, mean_delay_picosec );  // EDA tile meanoffset -> phase [degrees]
      
      // EDA tile index -> AAVS-ANTENNA 
      int aavs_antenna_idx  = tpmrx2edatileindex( eda_bf_index );
//      int aavs_antenna_idx  = eda_bf_index; // TEST 
      double cal_phase_tile = gCalSign * cal_phase[aavs_antenna_idx];
      printf("\n");
      printf("EDA tile index = %d\n",eda_bf_index);
      printf("\tAAVS antenna index = %d -> cal_phase = %.4f [deg]\n",aavs_antenna_idx,cal_phase_tile);
      
      double pointing_phase = cal_phase_tile + phase_deg;
      double norm_pointing_phase = pointing_phase;
      normalise_phase_m180_180( norm_pointing_phase );
      
      double eda_tile_point_phase = point_phase[aavs_antenna_idx];

// TEST :
//      eda_tile_point_phase -= 200;
      
      double diff_raw = (norm_pointing_phase - eda_tile_point_phase);
      double diff     = diff_raw;
      normalise_phase_m180_180( diff );
      
      double corr_delay = delay_picosec_from_phase_deg( gFreqMHz, diff );
      
      printf("\t%02d (aavs-ant = %02d): %+09.2f ps -> %+07.1f deg + cal_phase=%+06.1f [deg] = %+07.1f [deg] -> %+06.1f [deg] (normalised) vs. Source(%s) = %+06.1f [deg] -> DIFF = %+05.2f [deg] -> delay %.1f [ps]\n",
            eda_bf_index,aavs_antenna_idx,mean_delay_picosec,phase_deg,cal_phase_tile,pointing_phase,gSource.c_str(),norm_pointing_phase,eda_tile_point_phase,diff,corr_delay);
   
      if( fabs(diff) < 40 ){
          good_counter++;
      }
   
      ant_counter++;   
   }

   printf("Number [ps for grep :)] of agreeing (within 40deg) = %d\n",good_counter);      
   printf("360 deg is %.1f [ps], 10 deg is %.1f [ps], 50 deg is %.1f [ps]\n", delay_picosec_from_phase_deg( gFreqMHz, 360.0 ), delay_picosec_from_phase_deg( gFreqMHz, 10.0 ), delay_picosec_from_phase_deg( gFreqMHz, 50.0 ) );
}
