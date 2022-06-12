#!/bin/bash
rm *exe
make clean
make
cd examples/noBrownian/
./runSimulation
cd ../..
