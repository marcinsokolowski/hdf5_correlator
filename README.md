# hdf5_correlator
Data converter, beamformer, correlator and possibly more for hdf5 data from SKA-Low prototype stations

# requirements :

  msfitslib ( https://github.com/marcinsokolowski/msfitslib )

  sudo apt-get install libnova-dev libnova fftw2 fftw-dev libnova-0.16-0 libnova-dev fftw3-dev libhdf5-dev libcfitsio-dev

# may require fix in CMakeLists.txt
  find /usr -name H5Cpp.h
  Edit CMakeLists.txt and set path manuall to whatever the above find has found :
    include_directories(${FITSIO_INCLUDE_DIR} /usr/include/hdf5/serial/)
  
# build :

  git clone git@github.com:marcinsokolowski/hdf5_correlator.git
  cd hdf5_correlator
  mkdir build
  cd build
  cmake ..
  make 
  sudo make install

  cd /usr/local/bin/
  sudo ln -s path/hdf5_correlator/scripts/hdf2uvfits_sun.sh 
  sudo ln -s path/hdf5_correlator/scripts/hdf2uvfits_zenith.sh 
  sudo ln -s path/hdf5_correlator/scripts/hdf5_to_bin_all.sh 
  # where path is location of hdf5_correlator clone 

# WARNING : hdf5 might be a problem, but I'll try to solve it ASAP 

# Example of reading HDF5 file in C++ and correlation of 2 specified
# antennas given in :

  cd examples/read_hdf5
  mkdir build
  cmake ..
  make 

  