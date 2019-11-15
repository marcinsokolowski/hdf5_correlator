#ifndef _CAL_SOLUTIONS_H__
#define _CAL_SOLUTIONS_H__

// g++ calsolutions.cpp  -I$SRCDIR/cmn/baselib/ -D_UNIX -D_MAIN_ -o test_calsolutions $NDIR/slib/libbaselib.a $NDIR/slib/libmathlib.a -ldl `root-config --libs`
//
// g++ calsolutions.cpp  -I$SRCDIR/cmn/baselib/ -D_UNIX -c 

#include <vector>
#include <string>

class CCalSolFit : public std::vector<double> // coefficients of fitted polynomial (or any other function)
{
public :
   int    antenna_index;
   std::string fit_filename;
   double fit_chi2;
   double fit_residuals_rms;
   bool   antenna_ok;
   
   CCalSolFit();
   CCalSolFit( const CCalSolFit& right );
   
   CCalSolFit& operator=( const CCalSolFit& right );
   
   double value( double x );
   void   test( const char* outfile="test.txt", double start_x=0.00, double end_x=34000, double step=1.00 );
   
};

class CCalSolFits : public std::vector<CCalSolFit>
{
public :
   std::string filename;
 
   CCalSolFits();   
   int read( const char* filename="fit_file.txt" );
   int validate_fits( std::vector<int>* p_flagged_list=NULL, double max_chi2=1.00, double max_fit_rms=1.00, bool bFlagFromFitQuality=true );
   
   CCalSolFit* find_antenna( int antenna_index );
   double value( double x, int antenna_index );
};

#endif
