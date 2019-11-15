// Marcin Sokolowski - read hdf5 from aavs correlator and save to text files 

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <complex>

#include <H5Cpp.h>

// my :
// #include <spectrometer.h>
 #include <bg_globals.h>

// 
#include <myparser.h>
#include <mystring.h>
#include <mystrtable.h>

#include "hdf5_commons.h"

std::string file_list = "corr_file_list";
// unix time :
int gFileUxTime = -1;


void usage()
{
   printf("hdf5corr2others hdf5_file_list\n");
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ha:";
   int opt;

   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'a':
//            antenna1 = atol( optarg );
            break;

         case 'h':
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
   printf("Correlated file list = %s\n",file_list.c_str());
   printf("##############################################\n");
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

int read_hdf5corr_file( const char* hdf5_corr_file )
{
   string filename = hdf5_corr_file;

   printf("\n\n");
   printf("Reading file %s\n",hdf5_corr_file);   
   H5::H5File file( filename.c_str(), H5F_ACC_RDONLY ); 
   gFileUxTime = hdf5filename2uxtime( filename.c_str() );
   int count = file.getObjCount();
   hid_t file_id = file.getId();
   printf("Number of objects in file %s is %d , file_id = %x\n",filename.c_str(), count, (int)file_id );
   herr_t status =  H5Ovisit (file_id, H5_INDEX_NAME, H5_ITER_NATIVE, list_obj_iterate, NULL);
   if( status < 0 ) {
       printf("ERROR\n");
       // throw std::runtime_error("H5Ovisit returned with negative value which indicates a critical error");
   }

   H5::DataSet dataset = file.openDataSet( "correlation_matrix/data" );
   printf("Successfully opened dataset correlation_matrix/data\n");


   H5::DataSpace dataspace = dataset.getSpace();
   int rank = dataspace.getSimpleExtentNdims();

   hsize_t dims_out[HDF5_CORR_DIMS];
   dataspace.getSimpleExtentDims( dims_out, NULL);
   dataspace.selectAll();
   H5T_class_t type_class = dataset.getTypeClass();
   printf("Data dimensions = %d x %d x %d x %d, rank = %d, data_type = %d\n",(int)(dims_out[0]),(int)(dims_out[1]),(int)(dims_out[2]),(int)(dims_out[3]),rank,type_class);


   H5::DataType data_type = dataset.getDataType();
   hid_t data_type_id = data_type.getId();
   hid_t native_data_type_id = H5Tget_native_type( data_type_id , H5T_DIR_ASCEND);
   printf("Data type size = %d vs. sizeof(std::complex<float>) = %d , id = %d , native_data_type_id = %d\n",(int)(data_type.getSize()),(int)sizeof(std::complex<float>),(int)data_type_id,(int)native_data_type_id);

   int one_pol_data_size = dims_out[2];
   int total_data_size = dims_out[0]*dims_out[1]*dims_out[2]*dims_out[3];
//   std::vector< std::complex<float> > data( dims_out[0] , dims_out[1] , dims_out[2], dims_out[3] );
   std::vector< std::complex<float> > data( total_data_size );
   data[0] = std::complex<float>( +1000.00, -1000 );
//   std::vector< vis_t > data( one_pol_data_size );     
//   std::vector< float > data( total_data_size*2 );
   printf("data.size = %d , sizeof(data[0]) = %d, total_data_size = %d\n",(int)data.size(),(int)sizeof(data[0]),total_data_size);



   H5::DataSpace memspace( rank, dims_out ); // rank

/*   H5::CompType ctype(sizeof(std::complex<float>));
   ctype.insertMember( "real", HOFFSET(std::complex<float>,re), H5::PredType::NATIVE_FLOAT ); // crash H5T_IEEE_F32LE ); // H5::PredType::NATIVE_FLOAT
   ctype.insertMember( "imag", HOFFSET(std::complex<float>,im), H5::PredType::NATIVE_FLOAT ); // crash H5T_IEEE_F32LE ); // H5::PredType::NATIVE_FLOAT*/

/*   H5::CompType ctype2(sizeof(vis_t));
   ctype2.insertMember( "vis", HOFFSET(vis_t,xx), ctype ); 
   ctype2.insertMember( "vis", HOFFSET(vis_t,xy), ctype );
   ctype2.insertMember( "vis", HOFFSET(vis_t,yx), ctype );
   ctype2.insertMember( "vis", HOFFSET(vis_t,yy), ctype );*/


   dataset.read( data.data(), H5::DataType(native_data_type_id), memspace, dataspace );

//   float* data_f = new float[total_data_size*2];
//   dataset.read( data_f, H5::PredType::NATIVE_FLOAT, memspace, dataspace );   
//   printf("ODO : %.4f\n",data_f[0]);
   
//   H5std_string string_buffer;
//   dataset.read( string_buffer, ctype, memspace, dataspace );
//   printf("string_buffer = %s\n",string_buffer);

//   std::vector< std::complex<float> > data2( one_pol_data_size*4 );
//   dataset.read( data2.data(), ctype, memspace, dataspace );
//   dataset.read( data.data(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace );
//   printf("Errno = %d\n",errno);

   printf("TEST DATA :\n");
   float max_val = -1e6;
   int test_idx=0;
   for(test_idx=0;test_idx<=4;test_idx++){
       printf("data[%d] = %.2f / %.2f\n",test_idx,data[test_idx].real(),data[test_idx].imag()); 
   }
//  printf("data[%d] = %.2f \n",test_idx,data[test_idx]);
   
   test_idx = one_pol_data_size;
   printf("data[%d] = %.2f / %.2f\n",test_idx,data[test_idx].real(),data[test_idx].imag()); 
//   printf("data[%d] = %.2f\n",test_idx,data[test_idx]);
   
   for(int i=0;i<total_data_size;i++){
      if( data[i].real() > max_val ){
         max_val = data[i].real();
//      if( data[i] > max_val ){         
//          max_val = data[i];
      }
   }
   printf("max_val = %.4f\n",max_val);
}


void dump_hdf5corr_files( const char* filelist, int antenna=0 )
{
    vector<string> list;
    int n_files = bg_read_list( filelist, list );
    
    for(int i=0;i<n_files;i++){
       read_hdf5corr_file( list[i].c_str() );
    }    
}

int main(int argc,char* argv[])
{
   if ( argc>=2 ){
       file_list = argv[1];
   }
   parse_cmdline( argc , argv );
   print_parameters();

   dump_hdf5corr_files( file_list.c_str() );
}

