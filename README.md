# hdf5_correlator
Data converter, beamformer, correlator and possibly more for hdf5 data from SKA-Low prototype stations

# requirements :

  msfitslib ( https://github.com/marcinsokolowski/msfitslib )

  cfitsio libnova
  
# build :

  git clone git@github.com:marcinsokolowski/hdf5_correlator.git
  cd hdf5_correlator
  mkdir build
  cd build
  cmake ..
  make 
  sudo make install

# WARNING : hdf5 might be a problem, but I'll try to solve it ASAP 
