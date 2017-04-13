#!/bin/bash

# create build directory
if ! [ -d "build" ]; then
  mkdir build
fi

# start cmake in build directory
cd build

# compile source
cmake .. && make

# run kalman Filter
./ExtendedKF ../data/sample-laser-radar-measurement-data-1-lidar.txt test_output.txt

# remove test output file
rm test_output.txt
