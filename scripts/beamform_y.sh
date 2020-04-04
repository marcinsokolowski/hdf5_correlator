#!/bin/bash


delays="-22.43,19.51,0.00,6.42,-23.34,13.79,-35.23,29.49,63.43,35.36,25.36,-74.92,22.39,-36.29,6.18,-0.79,12.66,21.90,51.60,-30.20,28.64,36.39,-18.54,-41.12,-28.35,1.74,-13.90,-7.70,-1.81,52.14,0.00,30.52,13.43,-19.89,-26.27,21.53,-47.52,-24.26,-26.62,-47.42,9.53,10.71,1.58,-13.69,-59.88,5.45,-4.78,-33.25,-19.38,-32.22,63.57,31.87,22.98,-38.81,-29.25,21.09,11.85,-68.51,-30.09,-2.16,13.30,10.19,30.45,-38.31,15.14,-4.24,6.46,12.11,5.58,2.18,32.44,-18.61,18.70,25.48,16.26,30.03,-11.53,11.39,17.05,-3.95,25.91,14.54,22.24,17.52,-8.69,14.73,0.00,-17.98,-16.08,1.17,24.34,6.42,25.39,-36.94,-52.16,0.00,71.31,45.09,-16.58,-5.38,-35.17,-20.46,5.67,42.12,24.97,-5.77,25.95,9.20,-5.30,59.41,17.64,23.55,36.40,24.81,19.35,10.58,83.69,-57.64,-39.65,59.16,-21.15,38.34,1.48,-6.28,0.00,-32.88,14.18,119.23,6.74,18.29,26.61,22.11,11.10,0.00,26.77,22.57,-14.47,19.75,27.64,-9.48,22.86,52.12,-19.67,-18.34,-7.94,36.85,16.93,16.70,68.99,43.30,40.99,58.17,34.18,18.19,28.82,69.67,11.60,42.67,0.00,0.00,-43.32,-39.65,-38.87,28.91,-4.39,-11.23,41.47,-32.08,-17.60,12.68,-29.76,-46.42,-33.72,25.62,-4.05,18.97,25.96,-41.21,-22.56,-24.77,33.84,26.16,6.25,58.27,-32.76,0.64,12.46,3.09,-6.90,-15.85,-5.78,28.14,0.00,2.58,17.98,3.03,49.57,-58.78,14.58,-16.08,-4.64,10.68,-17.01,-16.69,-52.57,8.36,51.70,56.58,21.77,-21.53,-1.21,-31.19,-26.13,-6.81,-19.43,24.03,14.33,3.42,-22.13,-28.44,3.95,-19.68,18.87,-24.86,-14.70,2.75,-22.26,-21.11,6.23,0.00,9.47,24.35,-11.26,29.95,-32.98,23.41,-31.15,-9.46,-17.83,19.90,-19.46,-28.56,1.96,-15.18,26.95,1.98,-8.79,-11.97,-2.74,9.63,-5.69,7.39,-19.27,-5.97,-0.52,-26.25"
if [[ -n "$1" && "$1" != "-" ]]; then
   delays=$1
fi

dt=20200225
if [[ -n "$2" && "$2" != "-" ]]; then
   dt="$2"
fi

# tag="calsol20191129_ch204_HydA_SELECTED_PLUS"
tag="calsol20200224_ch126_Gero_X"
if [[ -n "$3" && "$3" != "-" ]]; then
   tag="$3"
fi


# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

# list passed :
# echo "nohup beamform_all.sh \"-C 1 -X $delays -K (0.00,75.392237) -E antenna_locations_eda2.txt\" ${dt}_24hours_eda2_256ant_${tag}.txt - hdf5_list > ${tag}.out 2>&1 &"
# nohup beamform_all.sh "-C 1 -X $delays -K (0.00,75.392237) -E antenna_locations_eda2.txt" ${dt}_24hours_eda2_256ant_${tag}.txt - hdf5_list > ${tag}.out 2>&1 &

# all hdf5 files :
ls *.hdf5 > hdf5_list


# (AZIM,ALT) = ( 207.0806 , 70.8457 ) [deg]
# az=200 # was error 226 
# az=153
# el=70.8457
# echo "beamform_all.sh \"-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1\" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1"
# beamform_all.sh "-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1

# (AZIM,ALT) = ( 207.0806 , 70.8457 ) [deg]
# az=207 # was error 226 
# az=153
# el=70.8457
# echo "beamform_all.sh \"-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1\" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1"
# beamform_all.sh "-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1


# az=207 # was error 226 
# el=30
# echo "beamform_all.sh \"-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1\" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1"
# beamform_all.sh "-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1

az=27.081946 # 360 - 153 
el=70.846771
freq_mhz=98.4375
echo "beamform_all.sh \"-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H $freq_mhz -p 1\" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1"
beamform_all.sh "-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H $freq_mhz -p 1" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1

# az=153 # 360 - 153 
# el=70.8457
# echo "beamform_all.sh \"-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1\" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1"
# beamform_all.sh "-C 1 -X $delays -K ($az,$el) -E antenna_locations_eda2.txt -e Geraldton_az${az}el${el}deg_timeseries_Y -H 98.4375 -p 1" 20200225_24hours_eda2_256ant_ch126_GeroBeam_az${az}el${el}deg_Y.txt - hdf5_list > az${az}el${el}deg_y.out 2>&1

