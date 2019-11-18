COMMON_LIBS=-lmsfits  -lcfitsio -lnova

HDF5_DIR=/usr/include/hdf5/include/
# HDF5_LIB=/opt/caastro/ext/anaconda3/lib/libhdf5_cpp.so /opt/caastro/ext/anaconda3/lib/libhdf5.so
# HDF5_LIB=/usr/lib/x86_64-linux-gnu/libhdf5_cpp.so 


# Ubuntu16 - ASUS :
# HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so //usr/lib/x86_64-linux-gnu/libhdf5_serial.so.10
HDF5_LIB=/usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.so /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.so

all : hdf5_correlator delay2phase hdf5corr2others
	chmod +x *! *.sh scripts/*.sh
#	cp hdf5_correlator $(BIGHORNS)/bin/
#	cp delay2phase $(BIGHORNS)/bin/
#	cp hdf5corr2others $(BIGHORNS)/bin/
#	cp *.sh $(BIGHORNS)/bin/
#	cp *! $(BIGHORNS)/bin/

hdf5_correlator : main.cpp hdf5_commons.cpp eda2_tpm_coefficients.cpp calsolutions.cpp
	g++ main.cpp hdf5_commons.cpp eda1_tpm_coefficients.cpp eda2_tpm_coefficients.cpp calsolutions.cpp -o hdf5_correlator $(OPT) -I ${HDF5_DIR} ${HDF5_LIB} -lfftw3 $(COMMON_LIBS) $(OPT) -D_UNIX

delay2phase : delay2phase.cpp
	g++ delay2phase.cpp -o delay2phase  ${HDF5_LIB}  -lfftw3 $(COMMON_LIBS) $(OPT) -D_UNIX

# test version which works : hdf5corr2others_tests.cpp
hdf5corr2others : hdf5corr2others.cpp hdf5_commons.cpp eda1_tpm_coefficients.cpp eda2_tpm_coefficients.cpp
	g++ hdf5corr2others.cpp hdf5_commons.cpp eda1_tpm_coefficients.cpp eda2_tpm_coefficients.cpp $(OPT) -o hdf5corr2others -I ${HDF5_DIR} ${HDF5_LIB} -lfftw3 $(COMMON_LIBS) $(OPT) -D_UNIX
# DEBUG :
#	g++ hdf5corr2others.cpp hdf5_commons.cpp -g -o hdf5corr2others -I /home/msok/mwa_software/hdf5/hdf5-1.10.5/src/ -I /home/msok/mwa_software/hdf5/hdf5-1.10.5/c++/src/ -I /opt/pi/dev/pisys/daq/src//cmn/baselib/ -I ../fitslib/ -I ../speclib/ ../speclib/libspeclib.a -lfftw3 -L/opt/caastro/bighorns//lib -lfitslib /opt/pi/dev/pisys/daq/ndir//slib/libbaselib.a /opt/pi/dev/pisys/daq/ndir//slib/libmathlib.a  `root-config --libs`  -lcfitsio -lnova  -g -D_UNIX /home/msok/mwa_software/hdf5/hdf5-1.10.5/src/.libs/libhdf5.a -lz -ldl /home/msok/mwa_software/hdf5/hdf5-1.10.5//c++/src/.libs/libhdf5_cpp.a /home/msok/mwa_software/hdf5/hdf5-1.10.5/src/.libs/libhdf5.a

clean :
	rm -f hdf5_correlator
	rm -f delay2phase
	rm -f hdf5corr2others
	
	
		