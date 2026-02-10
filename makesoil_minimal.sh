#!/usr/bin/env bash

# This script generates a netCDF map of derived soil properties file at resolutions of 30" (1km) and coarser.
# 
# The script uses the 1km SoilGrids rasters of the following soil physical properties:
#   sand (mass fractions)
#   clay
#   soil organic carbon
#   coarse fragments (volume fraction)
# it also uses a 250m raster of USDA soil name that informs the pedotransfer functions.
# 
# The following software is REQUIRED to run the script's programs: 
#    bc cURL GDAL GMT CDO NCO netCDF netCDF-Fortran

# ------------------------
# USER SETTINGS

# specify a directory for the output file NB this directory has to exist before specifying, 
# example:

outdir=/home/terraces/datasets/soils/WNA1km

# specify a target directory where the raw data is stored (or should be downloaded), 
# example:

datadir=/home/terraces/datasets/soils/soilgrids1km_v2020

# landmask file - this script requires a netCDF file containing the fraction of every gridcell that is land
# this can be generated using other datasets such as the G3WBM global water map :
# https://hydro.iis.u-tokyo.ac.jp/~yamadai/G3WBM/index.html

# -------------------------
# setup

landfrac=/home/terraces/datasets/hydrography/G3WBM/classfrac_WNA_1km.nc

wkt=/home/terraces/datasets/climate/NAwest1km/NAwest1km.wkt

gridspec=/home/terraces/datasets/climate/NAwest1km/NAwest1km.gridspec

proj="+proj=laea +lat_0=45 +lon_0=-115 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

bounds="-1190000 -1522000 1150000 1178000"

xres=1000
yres=1000

res="$xres $yres"

tmpdir=/dev/shm/$USER/makeclimate

outfile=$outdir/WNAsoils_v2026.nc

# -------------------------
# set to true if the raw data should be downloaded

getdata=false

# ------------------------
# 0) download the raw data (only if necessary)

if [ "$getdata" = true ]
then
  
  echo "Downloading raw data files from ISRIC"
  
  # 250m WRB code (575 MB):
  
  curl --output-dir $datadir -O https://files.isric.org/soilgrids/former/2017-03-10/data/TAXNWRB_250m_ll.tif
  
  # All other soil physical properties (1km Homolosine rasters about 190 MB each, need about 6 GB for all data)
  
  ./download_from_isric.sh $datadir
  
fi

# -----
# make sure the helper programs are up-to-date

make

# -----
# generate the output file

echo "making output file $outfile"

ncgen -4 -o $outfile soildata_minimal.cdl

# -----
# 5) paste USDA soil class into output

infile=$datadir/TAXOUSDA_250m_ll.tif

gdalwarp --quiet -overwrite -t_srs $wkt -te $bounds -wm 12G -multi -wo NUM_THREADS=16 -tr $res -tap -r mode -dstnodata 0 -of netCDF $infile $tmpdir/tmp.nc

# extrapolate to missing areas (urban and some water)

cdo -s setmisstonn $tmpdir/tmp.nc $tmpdir/tmp1.nc

# clip to landmask

./masklandmask_byte $landfrac $tmpdir/tmp1.nc

./pastesoilcode $tmpdir/tmp1.nc $outfile USDA

# ---------
# geodetic lon and lat

echo "generate coordinates"

gmt grd2xyz $tmpdir/tmp.nc | invproj -o $proj > $tmpdir/coords.bin 

/home/terraces/datasets/topography/WNA1km/pastecoords_bin $tmpdir/coords.bin $outfile

rm $tmpdir/coords.bin 

# -----
# 6) paste soil depth into output

# the original soil/regolith depth datasets are:
#   upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.tif 
#   upland_hill-slope_soil_thickness.tif 
#   hill-slope_valley-bottom.tif
# These files can be downloaded from NASA Earthdata at the following doi url:
# https://doi.org/10.3334/ORNLDAAC/1304
# I cannot find an easy way to automate this download as it requires a login

echo "make soil depth"

i=1
for infile in upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.tif upland_hill-slope_soil_thickness.tif hill-slope_valley-bottom.tif
do

  echo "reproject $i $infile"

  gdalwarp --quiet -overwrite -t_srs $wkt -te $bounds -wm 12G -multi -wo NUM_THREADS=16 -tr $res -tap -r average -of netCDF $datadir/$infile $tmpdir/depth$i.nc
  
  let i++

done

cp $tmpdir/depth3.nc $tmpdir/tmp.nc

./makethickness $tmpdir/depth1.nc $tmpdir/depth2.nc $tmpdir/depth3.nc $tmpdir/tmp.nc

# extrapolate to missing areas (urban and some water)

cdo -s setmisstonn $tmpdir/tmp.nc $tmpdir/tmp1.nc

# mask to landmask and paste into output

./mask-and-pack $tmpdir/tmp1.nc $outfile thickness


# $outfile

# -----
# 7) add coordinates

./pastecoords $tmpdir/tmp.nc $outfile

# -----
# 8) paste soil properties into file

for var in sand clay cfvo soc
do

  l=1

  for level in 0-5 5-15 15-30 30-60 60-100 100-200
  do
    
    infile=$datadir/$var"_"$level"cm_mean_1000.tif"
    
    gdalwarp --quiet -overwrite -t_srs $wkt -te $bounds -wm 12G -multi -wo NUM_THREADS=16 -tr $res -tap -r mode -of netCDF $infile $tmpdir/tmp.nc
    
    ncatted -a scale_factor,Band1,c,d,0.1 $tmpdir/tmp.nc
    
    # extrapolate basic input variables to missing areas (urban and some water)
    
    cdo -s setmisstonn $tmpdir/tmp.nc $tmpdir/tmp1.nc
    
    # clip to landmask
    
    ./masklandmask $landfrac $tmpdir/tmp1.nc
    
    # past into output

    ./ncpaste $tmpdir/tmp1.nc $outfile $var $l
  
    let l++    

  done
done

# -----
# 10) finish

rm $tmpdir/tmp.nc

# rm tmp*.nc  # to remove the intermediate soil depth files

echo "finished!"
