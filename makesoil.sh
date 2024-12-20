#!/opt/local/bin/bash

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
#    cURL GDAL GMT NCO netCDF netCDF-Fortran
# 
# The raw soil data could be downloaded in advance, otherwise a data download script is also provided

# ------------------------
# USER SETTINGS

# specify a directory for the output file NB this directory has to exist before specifying, 
# example:

# outdir=../global30minute
outdir=NA5km

# specify a target directory where the raw data is stored (or should be downloaded), 
# example:

datadir=/Volumes/Amalanchier/datasets/soils

# set to true if the raw data should be downloaded

getdata=false

# specify the output projection, extent, and resolution 

# specify a map projection using an EPSG code, proj4 string, or external file

# proj="EPSG:4326"  # example: unprojected lon-lat

proj="NAlaea.prj"

# specify the map extent and resolution

# extent="-180. -90.  180. 90."      #  <xmin> <ymin> <xmax> <ymax>

extent="-4350000. -3885000.  3345000. 3780000."      #  <xmin> <ymin> <xmax> <ymax>

res=5000.

min=30.                           # target resolution in MINUTES for lat-lon or METERS for projected grids

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

gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -of netCDF $infile tmp.nc

# get the dimensions of the target file

fileinfo=( $(gmt grdinfo -C tmp.nc?Band1) )

xlen=${fileinfo[9]}
ylen=${fileinfo[10]}

echo $xlen $ylen $res

# -----
# 3) create output file based on the dimensions of the input

outfile=$outdir/soils.nc

echo creating $outfile

sed -e "s/xlen/$xlen/g" -e "s/ylen/$ylen/g" soildata.cdl | ncgen -4 -o $outfile

# -----
# 4) paste WRB code into output

./pastesoilcode tmp.nc $outfile WRB

# 5) paste USDA soil class into output

infile=$datadir/TAXOUSDA_250m_ll.tif

gdalwarp --quiet -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -of netCDF $infile tmp.nc

./pastesoilcode tmp.nc $outfile USDA

# -----
# 6) paste soil depth into output

# infile=$datadir/hill-slope_valley-bottom.tif  # upland_valley-bottom_and_lowland_sedimentary_deposit_thickness.tif # upland_hill-slope_soil_thickness.tif
# 
# gdalwarp -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r average -of netCDF $infile tmp.nc

# infile=$datadir/average_soil_and_sedimentary-deposit_thickness.tif
# 
# gdalwarp -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r med -of netCDF $infile tmp.nc

# gdalwarp -overwrite -t_srs $proj -te $extent -wm 12G -multi -wo NUM_THREADS=16 -tr $res $res -tap -r mode -of netCDF $infile tmp.nc

# exit
# 
# ./pastesoilcode tmp.nc $outfile thickness

./makethickness $outfile

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

    ./ncpaste tmp.nc $outfile $var $l
  
    let l++    

  done
done

# -----
# 9) calculate derived soil properties

./soilcalc $outfile $outfile

# -----
# 10) finish

rm tmp.nc

echo "finished!"

