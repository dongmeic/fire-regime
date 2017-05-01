#!/bin/bash

# These are the inputs
Year=$1

DATADIR=/home2/dongmeic/fire/data/MCD64A1
OUTPUTDIR=/home2/dongmeic/fire/output
DoY_month=("001" "032" "060" "091" "121" "152" "182" "213" "244" "274" "305" "335") 	
Tiles=("h22v03" "h22v04" "h23v03" "h23v04" "h23v05" "h24v03" "h24v04" "h24v05" "h24v06" "h24v07" "h25v03" "h25v04" "h25v05" "h25v06" "h25v07" "h26v03" "h26v04" "h26v05" "h26v06" "h26v07" "h27v04" "h27v05" "h27v06" "h27v07" "h28v04" "h28v05" "h28v06" "h28v07" "h29v05" "h29v06" "h29v07")



# for each month
themonth=()
for DoY in `seq 0 11`
do
	# for each tile
	NumberOfTiles=`echo ${#Tiles[@]} - 1 | bc -l`
	Files=()
	for tile in `seq 0 $NumberOfTiles`
    	do
		Files+=("`ls $DATADIR/${Tiles[$tile]}/MCD64A1.A$Year${DoY_month[$DoY]}.${Tiles[$tile]}.051.?????????????.hdf`")

	done

	# extract the "burn date" layer
	layer_files=()
	for tile in `seq 0 $NumberOfTiles`
	do
		echo ${Files[$tile]}
		thelayer="HDF4_EOS:EOS_GRID:"${Files[$tile]}:"MOD_Grid_Monthly_500m_DB_BA:Burn Date"
		gdalwarp -of GTiff -tr 500 500 "$thelayer" $OUTPUTDIR/tmp.$tile.tif  > /dev/null
		echo $OUTPUTDIR/tmp.$tile.tif
		layer_files+=("$OUTPUTDIR"/"tmp.$tile.tif")
		
	done
	echo LAYER FILES : ${layer_files[@]}

	# mosaic to create a map of china
	outputfile=$OUTPUTDIR/tmp.BA_China_$Year_${DoY_month[$DoY]}.tif
		
	/packages/gdal/2.1.3/bin/gdal_merge.py  -of GTiff -o $outputfile ${layer_files[@]}

	# reproject from sinosuidal to wgs84 latlong
	themonth=`echo ${DoY} + 1 | bc`
	echo THE MONTH = $themonth

	finaloutputfile=$OUTPUTDIR/BA_China_$Year-$themonth.tif
	gdalwarp -t_srs "+proj=latlong +datum=WGS84" -te 73.125 17.476 135.703 53.956 -tr 0.0045 0.0045 $outputfile $finaloutputfile

	rm $OUTPUTDIR/*tmp*

done
