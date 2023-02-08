#!/bin/bash
rm *exe
make clean
make
cd examples/
./runSimulation
cd ..
