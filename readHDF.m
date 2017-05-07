close all; clear; clc;
fileinfo = hdfinfo('/Volumes/dongmeic-22/fire/raw/MCD64A1/h22v03/MCD64A1.A2000306.h22v03.006.2017012024847.hdf')
fileinfo.Attributes(1)

fileinfo = hdfinfo('/Volumes/dongmeic-22/fire/raw/CMG/C5/8day/hdf/MOD14C8H.2000305.005.01.hdf')
fileinfo.SDS(1)