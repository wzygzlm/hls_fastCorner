############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
open_project hls_fastCorner
set_top parseEvents
add_files fastFPGA_hls/src/fast_detector.cpp
add_files fastFPGA_hls/src/fast_detector.h
add_files -tb hls_fastCorner/src/mainapp.cpp
open_solution "solution4"
set_part {xc7z007sclg225-1}
create_clock -period 10 -name default
source "./hls_fastCorner/solution4/directives.tcl"
csim_design -compiler gcc
csynth_design
cosim_design
export_design -format ip_catalog
