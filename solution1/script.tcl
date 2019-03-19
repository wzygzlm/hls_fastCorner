############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
open_project hls_fastCorner
set_top insertion_cell_sort
add_files /home/minliu/github_prjs/hls_abmof/hls_fastCorner/src/fast_detector.cpp
add_files /home/minliu/github_prjs/hls_abmof/hls_fastCorner/src/fast_detector.h
add_files hls_fastCorner/src/insertion_cell_sort.h
add_files hls_fastCorner/src/sortHW.cpp
add_files -tb hls_fastCorner/src/sortTB.cpp
open_solution "solution1"
set_part {xc7z007sclg225-1}
create_clock -period 10 -name default
source "./hls_fastCorner/solution1/directives.tcl"
csim_design
csynth_design
cosim_design
export_design -format ip_catalog
