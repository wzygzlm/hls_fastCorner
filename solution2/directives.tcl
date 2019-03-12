############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
set_directive_loop_tripcount -min 16 -max 16 "FastDetector::isFeature/isFeatureOutterLoop"
set_directive_interface -mode ap_ctrl_hs "FastDetectorisFeature"
set_directive_interface -mode ap_fifo -depth 500 "parseEvents" data
set_directive_interface -mode s_axilite -register "parseEvents"
set_directive_interface -mode axis -register -register_mode both -depth 500 "parseEvents" eventSlice
set_directive_array_reshape -type block -factor 1 -dim 2 "FastDetectorisFeature" sae_
set_directive_loop_tripcount -min 16 -max 16 "FastDetectorisFeature/FastDetectorisFeature_label0"
set_directive_loop_tripcount -min 3 -max 6 "FastDetectorisFeature/FastDetectorisFeature_label1"
set_directive_loop_tripcount -min 4 -max 4 "FastDetectorisFeature/FastDetectorisFeature_label2"
set_directive_loop_tripcount -min 16 -max 16 "FastDetectorisFeature/isFeatureOutterLoop"
set_directive_loop_tripcount -min 12 -max 16 "FastDetectorisFeature/FastDetectorisFeature_label3"
set_directive_loop_tripcount -min 3 -max 7 "FastDetectorisFeature/FastDetectorisFeature_label4"
set_directive_loop_tripcount -min 5 -max 5 "FastDetectorisFeature/FastDetectorisFeature_label5"
set_directive_loop_tripcount -min 20 -max 20 "FastDetectorisFeature/FastDetectorisFeature_label6"
set_directive_pipeline -rewind "FastDetectorisFeature/FastDetectorisFeature_label0"
set_directive_pipeline -rewind "FastDetectorisFeature/FastDetectorisFeature_label3"
set_directive_pipeline -rewind "FastDetectorisFeature/FastDetectorisFeature_label4"
set_directive_pipeline -rewind "FastDetectorisFeature/FastDetectorisFeature_label1"
