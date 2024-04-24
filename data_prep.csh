#!/bin/csh

set input_data = "Lg_efficiency.txt"
set output_data = "lg_test.txt"

### Compile ray_trace C++ program using C++11 standards
echo "Compiling the ray_trace program..."
g++ ray_trace.cpp -o ray_trace -std=c++11
echo "Compilation of ray_trace complete."

### Extract event location, station location, and efficiency from input_data
awk '{print $1, $2, $3, $4, $5}' $input_data > temp

### Run ray_trace to produce inputs for Bayesian logit analysis 
echo "Running the ray_trace program..." 
./ray_trace temp $output_data 2 2 > mesh.out
echo "Processing completed. Data prepared and stored in $output_data."

