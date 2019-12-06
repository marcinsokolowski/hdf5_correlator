#include "mwaconfig.h"

/** @file mwaconfig.cpp
 * @brief Several functions to read the MWA meta data files.
 * 
 * Most functions were converted to C++ from corr2uvfits.c
 * written by Randall Wayth.
 * @author André Offringa offringa@gmail.com
 */

#define MWA_LATTITUDE -26.703319        // Array latitude. degrees North
#define MWA_LONGITUDE 116.67081         // Array longitude. degrees East
#define MWA_HEIGHT 377               // Array altitude. meters above sea level

#include <cstdlib>
#include <cmath>

#include <fstream>
#include <sstream>
#include <iostream>

#include <stdexcept>

#include "geometry.h"
// #include "lfiles_reader.h"

#define VEL_FACTOR  1.204         // the velocity factor of electic fields in RG-6 like coax.

extern int gDebugLevel;

double getJulianDay( int year, int month, int day, double hour )                      //Reference: "Astronomical Algorithms" by Jean Meeus
{
        //valid from 1900/3/1 to 2100/2/28
        if (month<=2) {month=month+12; year=year-1;}
        return (int)(365.25*year) + (int)(30.6001*(month+1)) - 15 + 1720996.5 + day + hour/24.0;
}


double getJulianDay( time_t ut_time )
{
        struct tm gmtm;
        gmtime_r(&ut_time,&gmtm);
        int year = gmtm.tm_year+1900;
        int month = gmtm.tm_mon + 1;

        double h = ( gmtm.tm_hour + double(gmtm.tm_min)/60.00 + double(gmtm.tm_sec)/3600.00 );
        double ret = getJulianDay( year, month, gmtm.tm_mday,h );
        return ret;
}


double getMJD( double ut_time )
{
     double fractional_second = ut_time - (time_t)ut_time;

     double jd = getJulianDay( (time_t)ut_time ) + (fractional_second / 86400.00);
     double mjd = jd - 2400000.5;

     return mjd;
}


void get_ymd_hms_ut( time_t ut_time, int& year, int& month, int& day,
                  int& hour, int& minute, int& sec )
{
   struct tm gmtm;
   gmtime_r( &ut_time , &gmtm );
   year = 1900 + gmtm.tm_year; // + 1900
   month = (gmtm.tm_mon+1);
   day = gmtm.tm_mday;
   hour = gmtm.tm_hour;
   minute = gmtm.tm_min;
   sec = gmtm.tm_sec;
}


time_t get_gmtime_from_string( const char* szGmTime )
{
   struct tm gmtime_tm;
   // sscanf( szGmTime, "%.4u%.2u%.2u_%.2u%.2u%.2u", &gmtime_tm.tm_year,&gmtime_tm.tm_mon,
   // &gmtime_tm.tm_mday,&gmtime_tm.tm_hour,&gmtime_tm.tm_min,&gmtime_tm.tm_sec);

   if( !strptime( szGmTime, "%Y%m%dT%H%M%S", &gmtime_tm ) ){
      if( !strptime( szGmTime, "%Y%m%d_%H%M%S", &gmtime_tm ) ){
         printf("ERROR : unrecognised time format |%s| -> please use either \%Y\%m\%dT\%H\%M\%S or \%Y\%m\%dT\%H\%M\%S\n",szGmTime);
         exit(-1);
      }
   }

   // gmtime_tm.tm_year -= 1900;
   // gmtime_tm.tm_mon--;
   time_t ret = timegm( &gmtime_tm );
        // time_t ret = timegm( &gmtime_tm );
   printf("Parsed UTC = %s -> ret = %d\n",szGmTime,ret);     
        
   return ret;
}


MWAHeader::MWAHeader() :
	nInputs(0),
	nScans(0),
	nChannels(32), // was 0
	correlationType(None),
	integrationTime(0.14155776), // 131072 * 1.08 usec / 1000000.000 [seconds]
	centralFrequencyMHz(159.375), // was 0.00
	freqChannel(204),
	bandwidthMHz(0.78125),
/*	raHrs(186.53066716/15.00), // was -99 // 186.532 -26.70331900
	decDegs(-26.70331900), // was -99
	haHrsStart(0.0), // was -99 
	refEl((M_PI/180.00)*90.00),*/
	raHrs(-1000.0),
	decDegs(-1000.0),
	haHrsStart(-1000.0),
	refEl(M_PI*0.5),
	refAz(0.0),
/*	year(2019),
	month(9),
	day(13),
	refHour(5),
	refMinute(11),
	refSecond(47.0),
	dateFirstScanMJD( getMJD( time_t(1568351507.00)) ),*/
        year(1970),
        month(1),
        day(1),
        refHour(0),
        refMinute(0),
        refSecond( 0.00 ),
        dateFirstScanMJD( 0.00 ),	
	conjugate(false),
	geomCorrection(true),
	fieldName(),
	polProducts("XXXYYXYY")
{
}

MWAHeaderExt::MWAHeaderExt() :
	gpsTime(0), observerName("Unknown"), projectName("Unknown"),
	gridName("Unknown"), mode("Unknown"), filename("Unknown"),
	hasCalibrator(false), hasGlobalSubbandGains(false),
	centreSBNumber(0),
	//fiberFactor(VEL_FACTOR),
	tilePointingRARad(0.0), tilePointingDecRad(0.0),
	dateRequestedMJD(0.0)
{
	for(size_t i=0; i!=16; ++i) delays[i] = 0;
	for(size_t i=0; i!=24; ++i) subbandGains[i] = 0;
}

/*double MWAHeader::GetStartDateMJD() const
{
	return Geometry::GetMJD(year, month, day, refHour, refMinute, refSecond);
}*/

/*double MWAHeader::GetDateFirstScanFromFields() const
{
	return 0.5*(integrationTime/86400.0) + Geometry::GetMJD(
		year, month, day, refHour, refMinute, refSecond);
}*/

void MWAConfig::ReadAntennaPositions(const std::string& filename, bool bConvertToXYZ) {
	std::ifstream file(filename.c_str());
	if(!file.good())
		throw std::runtime_error(std::string("Could not open ") + filename);
	const double arrayLattitudeRad = MWA_LATTITUDE*(M_PI/180.0);

	size_t antennaIndex = 0;
  /* Scan through lines. Convert east, north, height units to XYZ units */
  while(file.good()) {
		std::string line;
		std::getline(file, line);
		if(!line.empty() && line[0]!='#')
		{
			if(antennaIndex >= _antennae.size())
				_antennae.push_back(MWAAntenna());
			MWAAntenna& antenna = _antennae[antennaIndex];
			
			std::istringstream lineStr(line);
			double east, north, height;
			lineStr >> antenna.name >> east >> north >> height;
			if(lineStr.fail())
				throw std::runtime_error("Parsing antenna file failed on line: " + line);
			
			antenna.position[0] = east;
			antenna.position[1] = north;
			antenna.position[2] = height;
			
                        if( bConvertToXYZ ){
   			    Geometry::ENH2XYZ_local(east, north, height, arrayLattitudeRad, antenna.position[0], antenna.position[1], antenna.position[2]);
                        }
			
			antenna.stationIndex = antennaIndex;			
			antenna.tileNumber = atol( antenna.name.c_str()+3 );
			/*if(antenna.name.size() > 4 && antenna.name.substr(0, 4) == "Tile")
				antenna.tileNumber = atoi(antenna.name.substr(4).c_str());
			else
				antenna.tileNumber = 0;*/
			++antennaIndex;
		}
  }
  std::cout << "Read positions for " << antennaIndex << " tiles.\n";
}

int MWAConfig::polCharToIndex(char pChar)
{
	switch(pChar)
	{
		case 'X': case 'x':
		case 'R': case 'r':
		case 'I': case 'i':
			return 0;
		case 'Y': case 'y':
		case 'L': case 'l':
			return 1;
		default:
			throw std::runtime_error(std::string("Unknown pol char: ") + pChar);
	}
}

double MWAConfig::ArrayLattitudeRad()
{
	return MWA_LATTITUDE*(M_PI/180.0);
}

double MWAConfig::ArrayLongitudeRad()
{
	return MWA_LONGITUDE*(M_PI/180.0);
}

double MWAConfig::ArrayHeightMeters()
{
	return MWA_HEIGHT;
}

size_t MWAConfig::CentreSubbandNumber() const
{
	return round(ChannelFrequencyHz(_header.nChannels/2) / 1280000.0 + 0.5);
}

double MWAConfig::VelocityFactor()
{
	return VEL_FACTOR;
}


int MWAConfig::read_mapping_file( const char* filename )
{
   if( !filename || strlen(filename) == 0 ){
      printf("ERROR in InputMapping::read_mapping_file wrong filename parameter\n");
      return -1;
   }
   
/*   if( !CLFilesReader::DoesFileExist( filename ) ){
      printf("ERROR in InputMapping::read_mapping_file file %s does not exist\n",filename);
      return -1;
   }*/

   InputMapping input;

/* MyIFile in( filename );   
   while( pLine = in.GetLine( TRUE ) ){
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );

   }*/

   int max_antenna = -1;
   m_Inputs.clear();   
   FILE* map_file = fopen(filename,"r");
   char buff[1024];
   int lSize = 1024;
   
   while (1) {
      if(fgets(buff,lSize,map_file)==0)
         break;
      if(buff[0]=='#')
         continue;      

      char* ptr=NULL;
      char* search_ptr=buff;
      int col=0;    
      while( (ptr = strtok(search_ptr," \t")) ){
         search_ptr = NULL;
         if( 0 ){
            printf("ptr = %s\n",ptr);
         }

         switch( col ){
            case 0 :
                input.input = atol(ptr);
                break;

            case 1 :
                input.antenna = atol(ptr);
                break;
                
            case 2 :
                input.pol = ptr[0];
                break;     
                
            case 3:
                input.delta = atol(ptr);
                break;     

            case 4:
                input.flag = atol(ptr);
                break;     
                
            default :
                if( gDebugLevel>= 1 ){
                   printf("WARNING : column = %d not expected in file %s\n",col,filename);
                }
         }
         
         col++;
      }
      
      if( input.antenna > max_antenna ){
         max_antenna = input.antenna;
      }
      
      m_Inputs.push_back( input );         
   }   
   fclose(map_file);

   int antenna_count = (max_antenna+1);   
   printf("Read mapping for %d inputs (n_antennas = %d)\n",m_Inputs.size(),antenna_count);
   // m_CorrelatorMatrix
   int n_inputs = m_Inputs.size();
   std::vector<int> corr_line;
   corr_line.assign( n_inputs, -1 );
   m_CorrelatorMatrix.assign( n_inputs, corr_line );
   
   int index=0;
   for(int y=0;y<n_inputs;y++){
      for(int x=0;x<n_inputs;x++){
         if( y < x ){
            m_CorrelatorMatrix[y][x] = index;
            index++;
         }
      }
   }


   // fill antenna to input mapping :
   m_AntToInputMapping.clear();
   for(int ant=0;ant<antenna_count;ant++){
      sInputIndex ant2input_map;
   
      for(int i=0;i<m_Inputs.size();i++){
         InputMapping& input = m_Inputs[i];
         
         if( input.antenna == ant ){
            if( input.pol == 'X' || input.pol == 'x' ){
               ant2input_map.pol_x_input_index = i;
            }
            if( input.pol == 'Y' || input.pol == 'y' ){
               ant2input_map.pol_y_input_index = i;
            }            
         }
      }
      m_AntToInputMapping.push_back( ant2input_map );
   }
   
   return m_Inputs.size();
}

int MWAConfig::getInput( int antenna , int pol )
{
   InputMapping* pInput = getInputInfo( antenna , pol );
   if( pInput ){
      return pInput->input;      
   }
   
   return -1;
}

InputMapping* MWAConfig::getInputInfo( int antenna , int pol )
{
   char cPol = 'X';
   if( pol >= 1 ){
      cPol = 'Y';      
   }
   
   int inputIndex = antenna*2 + pol;
   if( inputIndex < m_Inputs.size() ){
      InputMapping& input = m_Inputs[ inputIndex ];
      
      if( input.antenna == antenna && input.pol == cPol ){
         return &( m_Inputs[ inputIndex ] );
      }else{
//         printf("ERROR : currently only very simple mapping is handled. Seems that the one requested is not implemented talk to MS\n");
         printf("WARNING : using experimental version of mapping ...\n");
                
         for(int i=0;i<m_Inputs.size();i++){
            InputMapping& input = m_Inputs[ i ];
            
            if( input.antenna == antenna && input.pol == cPol ){
               return (&input);           
            }
         }                  
      }
   }

   printf("ERROR: could not find information on antenna %d and polarisation = %d -> cannot continue with this error as results may be wrong !!!\n",antenna,pol);
   exit(-1);
   
   return NULL;
}


bool MWAConfig::IsAntennaFlagged( int antenna )
{
   InputMapping* pInputX = getInputInfo( antenna , 0 );
   InputMapping* pInputY = getInputInfo( antenna , 1 );
   
   if( pInputX && pInputY ){
      return ((pInputX->flag>0 ) || (pInputY->flag>0 ));
   }
   
   return true;
}

int MWAConfig::GetAntennaIndex( int input, int& pol )
{
   if( input < m_Inputs.size() ){
      int ant = m_Inputs[input].input;
      pol = 0;
      if ( m_Inputs[input].pol == 'X' || m_Inputs[input].pol == 'x' ){
         pol = 0;
      }
      if ( m_Inputs[input].pol == 'Y' || m_Inputs[input].pol == 'y' ){
         pol = 1;
      }
      
      return ant;
   }
   
   return -1;
}

int MWAConfig::GetBaselineIndex( int input1, int input2 )
{
   int baseline_index = -1;
   if( m_CorrelatorMatrix.size() > 0 && input2 < m_CorrelatorMatrix.size() && input1 < m_CorrelatorMatrix[0].size() ){
       if( input2 < input1 ){
          baseline_index = m_CorrelatorMatrix[input2][input1];
       }       
   }
   
   return baseline_index;
}
