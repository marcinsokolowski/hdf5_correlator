#ifndef _HDF5_COMMONS_H__
#define _HDF5_COMMONS_H__

// data sizes of hdf5 files :
#define HDF5_CORR_DIMS 4

time_t get_unixtime_from_local_string_localfunc( const char* szDTM );

double hdf5filename2uxtime( const char* filename, const char* format="channel_cont_0_%d_%d_0.hdf5", const char* format2="channel_cont_%d_%d_0.hdf5" );


#endif
