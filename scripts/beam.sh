#!/bin/bash

# beam.sh DATE channel pol delays tag 

dt=20190823 # WARNING : dt is here to be able to use it in call
if [[ -n "$1" && "$1" != "-" ]]; then
   dt="$1"
fi

ch=204
if [[ -n "$2" && "$2" != "-" ]]; then
   ch=$2
fi

pol=XX
if [[ -n "$3" && "$3" != "-" ]]; then
   pol=$3
fi

# delays="-95.35,17.18,0.00,-149.72,-166.31,-246.60,-211.97,-56.50,-136.87,-131.55,-2.16,67.56,10.75,-207.34,-151.58,-160.35,-91.02,-205.54,-182.94,-154.47,-19.52,-15.88,-202.25,-5.37,-181.92,-183.59,-12.15,-183.30,-27.88,-9.46,-174.90,-32.92,-196.60,-50.19,-160.58,-185.74,-4.95,-7.73,-208.19,-160.58,-17.80,-181.84,-200.26,-155.67,-16.46,-156.89,-133.35,-17.05,-172.59,-157.46,-144.83,-191.15,-26.48,-35.53,-167.67,-33.56,-164.47,-56.23,-26.57,-141.67,-193.32,-141.97,-173.17,-34.26,-17.65,-183.07,-7.90,-170.10,-53.78,-208.41,-36.22,-217.41,-180.73,-139.00,-174.44,-9.37,-39.18,-201.11,6.20,-39.48,-5.25,-182.20,-147.31,-31.78,-171.58,-215.64,-176.46,-174.61,12.01,-173.72,-24.10,-140.81,-187.53,-253.66,6.80,-133.26,-215.46,-138.35,-206.24,23.69,-155.91,-182.80,-167.76,-28.43,-173.75,32.31,-40.82,-163.76,-49.41,-121.54,-158.47,4.62,-85.04,-187.54,-176.69,-206.11,-149.20,-196.54,-162.97,-129.36,-32.44,-15.69,8.59,-203.84,-197.80,-8.56,-265.74,-249.35,-144.12,-162.06,-156.62,-1.97,-155.52,-34.30,-8.89,-163.90,-22.64,-154.50,-30.33,-5.78,-147.29,-151.99,-207.23,-43.46,-179.86,-53.57,-27.80,-144.81,-166.17,-154.02,-184.25,-25.09,-179.96,-156.53,5.09,-175.58,-1.49,-19.66,-39.24,-148.40,-34.01,-181.44,-65.19,-178.01,-181.57,-191.47,-145.95,-21.55,-158.68,-178.93,-171.37,-17.57,-142.66,-38.50,-55.52,11.96,-169.63,-40.82,-36.64,-20.97,-193.92,-180.71,-142.16,-130.70,-154.11,-23.93,-197.57,-158.46,-170.34,-174.27,-51.24,-131.43,-8.10,-31.61,-65.25,-5.98,-207.39,-15.04,2.31,-137.64,-164.13,-174.67,-188.38,-171.30,-156.82,5.03,-137.24,-163.66,-192.50,-28.76,11.61,-32.09,-139.94,-167.72,-213.74,-197.56,-228.33,-171.37,-93.65,-9.89,-8.59,-170.54,6.74,-31.47,17.60,-11.26,-184.02,-169.62,-73.96,-132.16,-170.34,16.90,-12.98,-52.44,-155.49,-47.76,-179.28,-206.09,-35.04,-31.27,-43.23,-231.23,-161.67,-187.54,-223.65,-174.20,-202.90,-107.85,-209.87,-31.49,-0.34,-178.70,-133.42,-212.77,5.20,-81.39"
delays=`get_last_calib.sh ${ch} ${pol} ${dt}` # uses last calibration, which does not have to be the best for the data do be beamformed use "get_last_calib.sh 204 XX 2019_12_11" to get a calibration for a specific date
echo "default delays = $delays"
if [[ -n "$4" && "$4" != "-" ]]; then
   delays=$4
fi

tag="GalCal"
if [[ -n "$5" && "$5" != "-" ]]; then
   tag="$5"
fi


# export PATH=~/msok/eda2/msok_scripts:$PATH
# export LD_LIBRARY_PATH=~/msok/eda2/lib:$LD_LIBRARY_PATH

echo "nohup beamform_all.sh \"-C 1 -X $delays\" ${dt}_24hours_ch${ch}_${tag}.txt  > ${tag}.out 2>&1 &"
nohup beamform_all.sh "-C 1 -X $delays" ${dt}_24hours_ch${ch}_${tag}.txt  > ${tag}.out 2>&1 &


# 20190812 cal. sols : 75.11,0.00,169.95,28.38,9.48,-28.50,-1.03,54.71,152.29,47.97,126.29,-63.03,-166.94,-148.61,-35.99,-50.46,51.25,-90.94,0.00,-147.12,-86.03,-8.34,-133.85,-157.40,-70.33,-39.96,-130.96,-125.25,-50.52,74.04,-29.49,0.00,32.42,-66.78,-147.42,-175.92,-177.75,66.23,138.61,-172.34,35.80,-99.99,-109.31,9.06,117.06,36.00,-51.83,-77.75,168.19,-68.53,-122.18,-143.60,50.38,-72.76,-68.02,-98.28,35.42,-38.38,-141.96,-127.16,21.60,106.35,50.15,-82.27,-23.71,161.50,-49.28,158.84,-57.19,97.25,-93.15,-0.36,112.09,-100.16,-43.06,-173.15,-58.69,109.21,-40.66,-128.50,38.25,0.00,-115.97,90.81,-147.72,152.40,-174.49,135.38,47.02,-73.73,-50.03,-139.04,-124.20,49.17,-44.27,-46.81,-133.28,-90.05,-149.31,78.35,121.65,-79.98,-132.30,-17.66,-45.80,141.60,110.25,160.21,78.86,-2.31,-116.00,179.39,122.07,39.04,-35.15,22.43,128.79,156.48,128.60,126.88,-93.06,-31.62,0.10,-3.08,118.77,-54.84,78.64,81.23,-107.74,106.37,-17.80,48.36,41.06,143.63,-23.62,-161.37,0.00,-168.19,51.13,160.06,-22.25,0.00,154.39,89.47,177.53,74.84,-87.99,-14.44,-111.58,147.16,149.32,-54.39,-8.35,119.81,-12.26,96.07,-99.03,-129.53,0.00,110.09,-69.08,74.71,142.55,-10.72,-106.84,28.46,-145.42,-66.60,-76.76,157.53,-172.66,100.21,-91.20,-93.06,157.45,-179.62,-11.89,71.69,11.92,-138.13,2.35,-15.14,35.29,-50.81,135.38,138.58,-148.42,-100.83,20.67,177.16,34.48,-8.55,0.00,155.95,103.08,151.52,-34.52,-9.05,-54.82,91.05,83.28,119.45,143.13,161.96,-101.02,106.39,-146.54,-132.20,136.47,15.63,0.00,-1.19,-69.33,-45.47,4.71,-93.00,-107.79,-52.48,-22.96,-23.41,23.21,0.00,37.33,64.61,-137.76,175.85,-10.03,59.93,18.00,118.86,-42.21,-93.53,86.01,0.00,62.21,-167.54,146.16,-51.34,162.80,-162.25,151.84,0.00,161.10,4.58,38.29,85.95,70.87,151.59,-64.06,100.37,156.05,23.49,0.00,-52.83,162.21,64.77
# 20190819 cal. sols : -93.88,0.00,0.00,-143.49,-160.65,163.62,-168.67,-113.97,-17.75,-123.22,-40.95,126.50,25.14,46.92,158.40,140.14,-45.57,173.17,0.00,117.55,178.80,-103.31,128.83,105.09,-168.25,-137.87,131.28,135.49,-148.42,-18.93,-123.31,0.00,-126.58,130.31,49.54,22.43,21.15,-92.03,-21.53,20.52,-133.63,90.58,85.71,-158.96,-49.78,-130.90,141.19,113.23,63.21,-176.29,133.50,111.60,-53.90,177.50,-174.76,151.58,-59.33,-137.19,119.35,144.56,-73.70,12.49,-46.48,-177.56,-118.16,68.04,-144.79,63.87,-155.87,-1.78,170.13,-97.94,13.75,163.45,-140.52,91.48,-157.26,17.23,-134.83,137.69,28.59,0.00,-125.40,79.36,-156.42,143.42,-150.65,125.31,34.61,-90.14,-64.27,-156.59,-140.57,36.37,-58.55,-59.07,146.32,-170.14,128.10,-3.04,40.86,-161.02,146.20,-96.88,-128.55,59.98,28.03,77.90,-1.19,-78.78,164.94,98.55,-100.48,174.60,99.07,155.66,-96.08,-69.39,-98.83,-98.57,41.71,103.07,134.66,131.00,-107.61,80.87,-144.32,-142.61,-130.40,90.95,-42.01,23.94,15.33,115.05,-46.93,175.16,0.00,164.84,29.72,137.26,-45.08,0.00,132.27,67.07,-132.87,123.07,-40.52,31.62,-60.38,-162.92,-162.31,-5.83,39.68,168.10,35.10,144.40,-50.82,-78.03,0.00,158.93,-103.45,40.06,107.31,-46.41,-140.05,-7.45,178.66,-101.14,-94.10,149.75,179.94,87.64,-99.33,-103.33,148.10,160.37,97.24,178.06,119.37,-30.17,113.18,96.83,144.49,56.90,-122.09,-112.53,-59.84,3.91,123.35,-74.58,140.26,93.87,0.00,16.35,-41.69,15.97,-148.57,-135.38,170.84,9.43,-72.70,-50.40,-2.90,-2.64,113.15,-41.55,74.14,71.85,-112.19,126.56,0.00,109.42,44.81,65.81,114.90,17.16,9.71,74.07,108.35,98.17,144.84,0.00,173.96,-161.93,-158.72,147.17,-31.77,44.77,-7.00,92.08,-68.66,-118.42,44.46,0.00,23.90,155.30,106.89,-87.57,123.45,159.30,-175.78,35.25,-165.53,39.15,76.54,122.49,105.58,-175.05,-29.62,133.50,-167.51,57.70,0.00,-14.78,-159.08,100.40
