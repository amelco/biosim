#!/bin/bash

gfortran -c variables.f90
gfortran variables.o BMG.f90 -o bmg
