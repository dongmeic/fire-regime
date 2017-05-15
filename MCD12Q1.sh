
# NOTE: The projection doesn't work in ACISS

#!/bin/bash

# These are the inputs
Year=$1

if [ $HOSTNAME = "d136-228.uoregon.edu" ]
	then
		DATADIR=/Volumes/dongmeichen/MCD12Q1
		OUTPUTDIR=/Volumes/dongmeichen/output/LC_China
		echo "Reading local folder"
	else
		DATADIR=/home2/dongmeic/MCD12Q1
		OUTPUTDIR=/home2/dongmeic/fire/output/LC_China
		echo "Reading data in ACISS"
fi

Tiles=("h22v03" "h22v04" "h23v03" "h23v04" "h23v05" "h24v03" "h24v04" "h24v05" "h24v06" "h24v07" "h25v03" "h25v04" "h25v05" "h25v06" "h25v07" "h26v03" "h26v04" "h26v05" "h26v06" "h26v07" "h27v04" "h27v05" "h27v06" "h27v07" "h28v04" "h28v05" "h28v06" "h28v07" "h29v05" "h29v06" "h29v07")

# for each tile
NumberOfTiles=`echo ${#Tiles[@]} - 1 | bc -l`
Files=()
for tile in `seq 0 $NumberOfTiles`
do
	Files+=("`ls $DATADIR/MCD12Q1.A$Year???.${Tiles[$tile]}.051.?????????????.hdf`")
done

# extract the "Land Cover Type 1 (IGBP)" layer
layer_files=()
for tile in `seq 0 $NumberOfTiles`
do
	echo ${Files[$tile]}
        thelayer='HDF4_EOS:EOS_GRID:'${Files[$tile]}':MOD12Q1:Land_Cover_Type_1'

	gdalwarp -of GTiff -tr 500 500 $thelayer $OUTPUTDIR/tmp.$tile.tif  > /dev/null
         # echo $OUTPUTDIR/tmp.$tile.tif
	layer_files+=("$OUTPUTDIR/tmp.$tile.tif")
		
done	
#echo LAYER FILES : ${layer_files[@]}

# mosaic to create a map of china
outputfile="$OUTPUTDIR/tmp_LC_China_$Year.tif"
		
gdal_merge.py  -of GTiff -o "$outputfile" "${layer_files[@]}"

# reproject from sinusoidal to wgs84 latlong
finaloutputfile="$OUTPUTDIR/LC_China_$Year.tif"
gdalwarp -s_srs "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" -t_srs "+proj=latlong +datum=WGS84" -te 73 17.5 135.5 54 -tr 0.005 0.005 "$outputfile" "$finaloutputfile"

rm $OUTPUTDIR/*tmp*



