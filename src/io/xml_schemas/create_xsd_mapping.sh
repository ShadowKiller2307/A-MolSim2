#!/bin/bash

# get input file
file=$1

xsdcxx cxx-tree --std c++11 --hxx-suffix .h --cxx-suffix .cpp --generate-doxygen --generate-serialization $file
