#!/usr/bin/bash

# This script generates a netCDF map of derived soil properties file at resolutions of 30" (1km) and coarser.
# 
# The script uses the 1km SoilGrids rasters of the following soil physical properties:
#   sand (mass fractions)
#   silt
#   clay
#   organic matter
#   coarse fragments (volume fraction)
# it also uses a 250m raster of WRB soil name that informs the pedotransfer functions.
# 
# The following software is REQUIRED to run the script's programs: 
#    bc cURL GDAL GMT CDO NCO netCDF netCDF-Fortran
# 
# The raw soil data could be downloaded in advance, otherwise a data download script is also provided

# ------------------------
# USER SETTINGS

# specify a directory for the output file NB this directory has to exist before specifying, 
# example:

# outdir=../global30minute
outdir=/home/terraces/datasets/soils/NA5km

# specify a target directory where the raw data is stored (or should be downloaded), 
# example:

datadir=/home/terraces/datasets/soils/soilgrids1km_v2020

# landmask file - this script requires a netCDF file containing the fraction of every gridcell that is land
# this can be generated using other datasets such as the G3WBM global water map :
# https://hydro.iis.u-tokyo.ac.jp/~yamadai/G3WBM/index.html

# landfrac=$datadir/classfrac_30m.nc

landfrac=/home/terraces/datasets/hydrography/G3WBM/NAlandfrac.nc

# set to true if the raw data should be downloaded

getdata=false

# specify the output projection, extent, and resolution 

# specify a map projection using an EPSG code, proj4 string, or external file

# proj="EPSG:4326"  # example: unprojected lon-lat

proj="/home/terraces/datasets/climate/climateNA/v730/prismdat/NAlaea.prj"

# specify the map extent and resolution

# extent="-180. -90.  180. 90."      #  <xmin> <ymin> <xmax> <ymax>

extent="-4350000. -3885000.  3345000. 3780000."      #  <xmin> <ymin> <xmax> <ymax>

res=5000.

# min=30.                           # target resolution in MINUTES for lat-lon or METERS for projected grids

if [ $proj == "EPSG:4326" ]
then
  res=`echo "$min / 60" | bc -l`    # convert to degrees
fi

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
# 2) decimate or project the WRB soil code raster into the target map domain and projection to retrieve the output file dimensions

infile=$datadir/TAXNWRB_250m_ll.tif

gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -dstnodata 0 -of netCDF $infile tmp.nc

# get the dimensions of the target file

bounds=( $(gmt grdinfo -C tmp.nc?Band1) )

xlen=${bounds[9]}
ylen=${bounds[10]}

xmin=${bounds[1]}
xmax=${bounds[2]}

ymin=${bounds[3]}
ymax=${bounds[4]}

echo $xlen $ylen $res $xmin $xmax $ymin $ymax

# -----
# 3) create output file based on the dimensions of the input

outfile=$outdir/soils.nc

echo creating $outfile

sed -e \
's/xlen/'$xlen'/g 
 s/ylen/'$ylen'/g
 s/xmin/'$xmin'/g
 s/xmax/'$xmax'/g 
 s/ymin/'$ymin'/g
 s/ymax/'$ymax'/g' \
soildata_template.cdl > soildata.cdl

ncgen -4 -o $outfile soildata.cdl

# -----
# 4) paste WRB code into output

# extrapolate to missing areas (urban and some water)

cdo -s setmisstonn tmp.nc tmp1.nc

# clip to landmask

./masklandmask_byte $landfrac tmp1.nc

./pastesoilcode tmp1.nc $outfile WRB

# -----
# 5) paste USDA soil class into output

infile=$datadir/TAXOUSDA_250m_ll.tif

gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -dstnodata 0 -of netCDF $infile tmp.nc

# extrapolate to missing areas (urban and some water)

cdo -s setmisstonn tmp.nc tmp1.nc

# clip to landmask

./masklandmask_byte $landfrac tmp1.nc $outfile

./pastesoilcode tmp1.nc $outfile USDA

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

  gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r average -of netCDF $datadir/$infile depth$i.nc
  
  let i++

done

./makethickness depth1.nc depth2.nc depth3.nc $outfile

# -----
# 7) add coordinates

./pastecoords tmp.nc $outfile

# -----
# 8) paste soil properties into file

for var in sand silt clay cfvo soc bdod
do

  l=1

  for level in 0-5 5-15 15-30 30-60 60-100 100-200
  do
    
    infile=$datadir/$var"_"$level"cm_mean_1000.tif"
    
    gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -of netCDF $infile tmp.nc
    
    ncatted -a scale_factor,Band1,c,d,0.1 tmp.nc
    
    # extrapolate basic input variables to missing areas (urban and some water)
    
    cdo -s setmisstonn tmp.nc tmp1.nc
    
    # clip to landmask
    
    ./masklandmask $landfrac tmp1.nc
    
    # past into output

    ./ncpaste tmp1.nc $outfile $var $l
  
    let l++    

  done
done

# -----
# 9) calculate derived soil properties

./soilcalc $outfile $outfile

# -----
# 10) finish

rm tmp.nc

# rm tmp*.nc  # to remove the intermediate soil depth files

echo "finished!"
