#!/bin/bash

if [ $# -eq 0 ] 
then
    echo "sh irun.sh [runmun]"
    exit
fi

#for i in $(seq -f%04g 6 76)
for i in $(seq -f%04g $1 $1)
do
#  echo "./r2root /run/media/jing.li/May312018/user/1769_data1x/Run00${i}HFC/HFC.dat ../1769_root/skim00${i}.root CC.cal ../1769_skim/skim00${i}.dat"
  ./gsort ../1563_root/run${i}.root ../1563_gsort/run${i}.root
done
