# to run year 2002 second half

#!/bin/bash

# These are the inputs
Year=$1

DATADIR=/Users/dongmeichen/CMG/hdf
OUTPUTDIR=/Users/dongmeichen/CMG/output/CloudCorrFirePix
#OUTPUTDIR=/Users/dongmeichen/CMG/output/CorrFirePix # 0
Types=("MOD" "MYD")

# for every 8 days
the8daysDoY=()
#for i in $(seq -f "%03g" 185 8 365)
for i in $(seq -f "%02g" 7 1 12)
do
	echo $i
	# for each type
	NumberOfTypes=`echo ${#Types[@]} - 1 | bc -l`
	Files=()
	for type in `seq 0 $NumberOfTypes`
    	do
		#Files+=("`ls $DATADIR/${Types[$type]}14C8H.$Year$i.005.01.hdf`")
		Files+=("`ls $DATADIR/${Types[$type]}14CMH.$Year$i.005.01.hdf`")

	done

	# extract the "CloudCorrFirePix" layer
	layer_files=()
	for type in `seq 0 $NumberOfTypes`
	do
		echo ${Files[$type]}
		thelayer="HDF4_SDS:UNKNOWN:"${Files[$type]}:"1"
		#thelayer="HDF4_SDS:UNKNOWN:"${Files[$type]}:"0"
		gdal_translate -of GTiff $thelayer $OUTPUTDIR/tmp.$type.tif  > /dev/null
		echo $OUTPUTDIR/tmp.$type.tif
		layer_files+=("$OUTPUTDIR"/"tmp.$type.tif")
		
	done
	echo LAYER FILES : ${layer_files[@]}

	# merge two layers
	outputfile=$OUTPUTDIR/tmp.CMG_$Year_$i.tif
		
	gdal_merge.py  -of GTiff -o $outputfile ${layer_files[@]}

	# project to wgs84 latlong
	the8daysDoY=`echo ${i} | bc`
	echo THE WEEK = $the8daysDoY

	finaloutputfile=$OUTPUTDIR/CMG_Global_$Year-$the8daysDoY.tif
	gdal_translate -a_srs "+proj=latlong +datum=WGS84" $outputfile $finaloutputfile
	
	rm $OUTPUTDIR/*tmp*

done