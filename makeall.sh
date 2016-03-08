#!/bin/bash
cd util/Release
make all -j
cd ..
cd ..
cd Simulator/Release
make all -j
cd ..
cd ..
cd Main/Release
make clean
make all -j
cd ..
cd ..