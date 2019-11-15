#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

// #include <string>
// #include <vector>
// #include <complex>

// #include <H5Cpp.h>

// my :
// #include <spectrometer.h>
#include <bg_date.h>

// 
// #include <myparser.h>
// #include <mystring.h>
// #include <mystrtable.h>


time_t get_unixtime_from_local_string_localfunc( const char* szDTM )
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



double hdf5filename2uxtime( const char* filename, const char* format /*="channel_cont_0_%d_%d_0.hdf5"*/, const char* format2 )
{
   int dtm = 0;
   int local_seconds = 0;

   // ERROR : could not parse date and time from file name channel_cont_20190718_71727_0.hdf5   
   // ERROR : could not parse date and time from file name channel_cont_20190718_83942_0.hdf5 using formatter |channel_cont_0_%d_%d_0.hdf5|
   if( sscanf( filename, format, &dtm, &local_seconds ) != 2 ){      
      printf("WARNING : could not parse date and time from file name %s using formatter |%s| -> trying alternative |%s|\n",filename,format,format2);      
      if( sscanf( filename, format2, &dtm, &local_seconds ) != 2 ){ 
         printf("ERROR : could not parse date and time from file name %s using formatter |%s| -> unknown filename format, date/time unknown\n",filename,format2);
         return -1;
      }
   }
   
   char szMidnight[128];
   sprintf(szMidnight,"%d_000000", dtm );
   time_t uxtime = get_unixtime_from_local_string_localfunc( szMidnight );
   
   printf("TEST uxtime = %d, local_seconds = %d\n",(int)uxtime,local_seconds);
   
   return (uxtime + local_seconds);
}
