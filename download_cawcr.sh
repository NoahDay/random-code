#!/bin/bash

WGET='/usr/bin/wget'

YEAR='2010'
MONTH='01'
MONTHEND='02' # MONTH + 1
DAY='01'
StartSeq='0'
EndSeq='12'

NCSS='https://data-cbr.csiro.au/thredds/ncss/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded'
MODEL='ww3.glob_24m.'

VARS="var=dir&var=fp&var=hs"
Subset='disableLLSubset=on&disableProjSubset=on'
SPATIAL='horizStride=1'
TimeStride='timeStride=1'

for PlusMonth in `seq $StartSeq $EndSeq`; do

  MyTime=`gdate -d "$YEAR-$MONTH-$DAY +$PlusMonth months" +%Y-%m-%dT%H:%M:%SZ`
  echo "Downloading CAWCR:  "$MyTime
  TimeStart="time_start=$MyTime"
  MyTimeEnd=`gdate -d "$YEAR-$MONTHEND-$DAY +$PlusMonth months -1 days" +%Y-%m-%dT%H:%M:%SZ`
  TimeEnd="time_end=$MyTimeEnd"

  grid="gridded_ww3.glob_24m"
  YearOut="`echo $MyTime | cut -d '-' -f 1`"
  MonthOut="`echo $MyTime | cut -d '-' -f 2`"
  OutFile=$grid"."$YearOut$MonthOut".nc"
  TimeStart="`echo $TimeStart | cut -d 'T' -f 1`"
  TimeEnd="`echo $TimeEnd | cut -d 'T' -f 1`"
  TimeStartSuf='T00%3A00%3A00Z'
  TimeEndSuf='T23%3A00%3A00Z'
  TimeStart="$TimeStart$TimeStartSuf"
  TimeEnd="$TimeEnd$TimeEndSuf"


  URL="$NCSS/$MODEL$YearOut$MonthOut.nc?$VARS&$Subset&$SPATIAL&$TimeStart&$TimeEnd&$TimeStride&addLatLon=true&accept=netcdf"

#  URL="https://data-cbr.csiro.au/thredds/ncss/catch_all/CMAR_CAWCR-Wave_archive/CAWCR_Wave_Hindcast_aggregate/gridded/ww3.glob_24m.201612.nc?var=dir&var=fp&var=hs&disableLLSubset=on&disableProjSubset=on&horizStride=1&time_start=2016-12-01T00%3A00%3A00Z&time_end=2016-12-31T23%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf"

  if [ -s $OutFile ]; then
        echo "[warning] File $OutFile exists (skipping)"
  else
        wget -O $OutFile "$URL"
  fi
done
