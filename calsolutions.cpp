#include "calsolutions.h"

#include <math.h>

#include <myfile.h>
#include <myparser.h>
#include <mystring.h>  
#include <mystrtable.h>
#include <basestructs.h>


CCalSolFit::CCalSolFit()
: antenna_index(-1),fit_chi2(-1),fit_residuals_rms(-1),antenna_ok(false)
{

}

CCalSolFit::CCalSolFit( const CCalSolFit& right )
{
   (*this) = right;
}

CCalSolFit& CCalSolFit::operator=( const CCalSolFit& right )
{
   (std::vector<double>&)(*this) = (std::vector<double>&)right;

   antenna_index = right.antenna_index;
   fit_filename  = right.fit_filename;
   fit_chi2      = right.fit_chi2;
   fit_residuals_rms = right.fit_residuals_rms;
   antenna_ok    = right.antenna_ok;
   
   return (*this);
}

double CCalSolFit::value( double x )
{
   double sum = 0.00;
   double x_power = 1.00;
   
   for(int i=0;i<size();i++){
      double par = (*this)[i];
      
      sum += x_power * par;
      
      x_power = x_power * x;
   }
   
   return sum;
}

void CCalSolFit::test( const char* outfile, double start_x, double end_x, double step )
{
   FILE* out_f = NULL;
   if( outfile ){
      out_f = fopen( outfile, "w" );
   }
   
   double x = start_x;
   while (x <= end_x ){
      double val = value( x );
      fprintf(out_f,"%.2f %.8f\n",x,val);      
      
      x += step;
   }
   
   if( out_f ){
      fclose(out_f);
   }
}
   
CCalSolFits::CCalSolFits()
{}


int CCalSolFits::read( const char* filename )
{
   MyFile infile(filename);
   const char* pLine=NULL;
   bool bFirstAntenna = true;
   int offset = 0;
   
   clear();
   while( pLine = infile.GetLine( TRUE ) ){
      bool bComment = false;
      if( mystring::get_first_non_white( pLine )=='#' ){
         bComment = true;
         continue;
      }
                        
      MyParser pars=pLine;
      CMyStrTable items;
      pars.GetItems( items );
      
      int ant_index = atol( items[0].c_str() );
      if( bFirstAntenna ){
         if( ant_index > 0 ){
            offset = ant_index;
         }
         bFirstAntenna = false;
      }
      
      CCalSolFit tmp;
      tmp.antenna_index = atol( items[0].c_str() ) - offset; // -1 to have them from 0 !!!
      tmp.fit_filename  = items[1].c_str();
      tmp.fit_chi2      = atof( items[2].c_str() );
      tmp.fit_residuals_rms = atof( items[3].c_str() );
      
      // read coefficients of the fitted function :
      for(int par=4;par<items.size();par++){
         tmp.push_back( atof(items[par].c_str()) );
      }      
      push_back( tmp );
   }      

   validate_fits();

   return size();
}


int CCalSolFits::validate_fits( std::vector<int>* p_flagged_list /*=NULL*/ , double max_chi2, double max_fit_rms, bool  bFlagFromFitQuality /*=true*/ )
{
   printf("CCalSolFits::validate_fits - validating fit results and potentially flagging bad antennas (if chi2 or rms exceed thresholds %.2f and %.2f respectively\n",max_chi2,max_fit_rms);
   int flagged_count = 0;
   
   for(int i=0;i<size();i++){
      CCalSolFit& ant_fit = (*this)[i];
      
      if( ant_fit.fit_chi2 > max_chi2 || fabs(ant_fit.fit_chi2)<0.00000001 ||ant_fit.fit_residuals_rms > max_fit_rms ){ // if chi2==0 or chi2>1 or rms>1 -> flag antenna
         ant_fit.antenna_ok = false;
         flagged_count++;
         
         if( p_flagged_list ){
            bool found = false;
            for(int k=0;k<p_flagged_list->size();k++){
               if( (*p_flagged_list)[k] == ant_fit.antenna_index ){
                  found = true;
                  break;
               }
            }
            
            if( !found ){
               if(  bFlagFromFitQuality ){
                  printf("WARNING : flagging antenna %d due to bad fit ( chi2 = %.2f , rms = %.2f )\n",ant_fit.antenna_index,ant_fit.fit_chi2,ant_fit.fit_residuals_rms);
                  p_flagged_list->push_back( ant_fit.antenna_index );
               }else{
                  printf("WARNING : antenna %d has bad fit ( chi2 = %.2f , rms = %.2f ), but flagging is not required !\n",ant_fit.antenna_index,ant_fit.fit_chi2,ant_fit.fit_residuals_rms);
               }
            }
         }
      }
   }
   
   return flagged_count;
}

CCalSolFit* CCalSolFits::find_antenna( int antenna_index ){
   for(int i=0;i<size();i++){
      if( (*this)[antenna_index].antenna_index == antenna_index ){
         return &((*this)[antenna_index]);
      }
   }
   
   return NULL;
}

double CCalSolFits::value( double x, int antenna_index )
{
   if( antenna_index>=0 && antenna_index < size() ){
      if( (*this)[antenna_index].antenna_index == antenna_index ){
         return (*this)[antenna_index].value( x );
      }
   }else{
      CCalSolFit* pAntFit = find_antenna( antenna_index );
      if( pAntFit ){
         if( pAntFit->antenna_index == antenna_index ){
            return pAntFit->value( x );
         }else{
            printf("ERROR in code : antenna index does not match (%d != %d)\n",pAntFit->antenna_index,antenna_index);
            exit(-1);
         }
      }
   }
   
   return (0.00/0.00); // NaN 
}


#ifdef _MAIN_
// g++ calsolutions.cpp  -I$SRCDIR/cmn/baselib/ -D_UNIX -D_MAIN_ -o test_calsolutions $NDIR/slib/libbaselib.a $NDIR/slib/libmathlib.a -ldl
int main( int argc,char* argv[] )
{
    CCalSolFits solutions;
    solutions.read();
    
    int ant_index = 0;
    if( argc >= 2 ){
       ant_index = atol( argv[1] );
    }
    char szOutTestFile[128];
    sprintf(szOutTestFile,"test_%03d.txt",ant_index);
    
    for(int i=0;i<solutions.size();i++){
       CCalSolFit& ant_fit = solutions[i];
       
       if( ant_fit.antenna_index == ant_index ){
          ant_fit.test( szOutTestFile );
       }       
    }
    
}
#endif

