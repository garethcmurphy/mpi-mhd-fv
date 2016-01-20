#!/bin/bash


LIMIT=9
for ((a=1; a <= LIMIT ; a++))  # Double parentheses, and "LIMIT" with no "$".
do
  convert -crop 0x0 chombo_plt_cnt_000$a.hdf5.ppm chombo_plt_cnt_000$a.hdf5.jpeg
  convert -rotate 90 chombo_plt_cnt_000$a.hdf5.jpeg chombo_plt_cnt_000$a.hdf5.ppm
done                           # A construct borrowed from 'ksh93'.
LIMIT=99
for ((a=10; a <= LIMIT ; a++))  # Double parentheses, and "LIMIT" with no "$".
do
  convert -crop 0x0 chombo_plt_cnt_00$a.hdf5.ppm chombo_plt_cnt_00$a.hdf5.jpeg
  convert -rotate 90 chombo_plt_cnt_00$a.hdf5.jpeg chombo_plt_cnt_00$a.hdf5.ppm
done                           # A construct borrowed from 'ksh93'.


