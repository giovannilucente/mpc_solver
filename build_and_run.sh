#!/bin/bash

cd build
cmake ..
make
./mpc_solver
cd ..
python3 plot_trajectories.py
