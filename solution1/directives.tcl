############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
set_directive_interface -mode s_axilite -register "parseEvents"
set_directive_interface -mode ap_fifo -depth 500 "parseEvents" data
set_directive_interface -mode ap_fifo -depth 500 "parseEvents" eventSlice
set_directive_loop_tripcount -min 12 -max 16 "FastDetectorisFeature/FastDetectorisFeature_label3"
set_directive_loop_tripcount -min 4 -max 8 "FastDetectorisFeature/FastDetectorisFeature_label4"
set_directive_loop_tripcount -min 5 -max 5 "FastDetectorisFeature/FastDetectorisFeature_label5"
set_directive_loop_tripcount -min 10 -max 13 "FastDetectorisFeature/FastDetectorisFeature_label0"
set_directive_loop_tripcount -min 3 -max 6 "FastDetectorisFeature/FastDetectorisFeature_label1"
set_directive_loop_tripcount -min 4 -max 4 "FastDetectorisFeature/FastDetectorisFeature_label2"
set_directive_array_reshape -type complete -dim 1 "min/label1/label1" inArr
set_directive_loop_tripcount -min 0 -max 16 "insertion_sort/L2"
set_directive_dataflow "fastCornerHW"
set_directive_resource -core RAM_2P_LUTRAM "rwSAE" saeHW
