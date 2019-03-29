############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
open_project hls_fastCorner
set_top tempSorted
add_files hls_fastCorner/src/fastCorner.cpp
add_files hls_fastCorner/src/fast_detector.cpp
add_files hls_fastCorner/src/temp.cpp
add_files -tb hls_fastCorner/src/test.cpp
open_solution "solution1"
set_part {xc7z007sclg225-1}
create_clock -period 10 -name default
#source "./hls_fastCorner/solution1/directives.tcl"
csynth_design
