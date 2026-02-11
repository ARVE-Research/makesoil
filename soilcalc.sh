#!/usr/bin/env bash

# this script prepares an empty output file for derived soil properties
# that are calculated using pedotransfer functions in soilcalc.f90

# create a soil input file using makesoil.sh before running this script

# load the netCDF-Fortran module (if necessary on your system)

# module load netCDF-Fortran

# update the program

make

# get the dimension sizes from a ncdump of the input file and create an empty output file

template=soildata_template_proj.cdl

infile=${1}   # first argument on the command line after the command

# file name of the output file is input_soilcalc.nc

tmp=${infile##*/}  # strip the directory path from the input file

output=${tmp%%.*}_soilcalc.nc  # strip the file suffix and add a new one

# extract the dimensions of the input file from an ncdump of the file header

xlen=`ncdump -h $infile | egrep -Eo 'x = [0-9]+' | grep -o '[0-9]*'`

ylen=`ncdump -h $infile | egrep -Eo 'y = [0-9]+' | grep -o '[0-9]*'`

# replace the placeholder variables in the metadata template with the actual file dimensions and generate the file

echo "generating output file $output with size $xlen x $ylen"

sed -e "s/XLEN/$xlen/g" -e "s/YLEN/$ylen/g" $template | ncgen -4 -o $output

# calculate the derived soil properties

./soilcalc $infile $output
