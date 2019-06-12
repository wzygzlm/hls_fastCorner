############################################################
## This file is generated automatically by Vivado HLS.
## Please DO NOT edit it.
## Copyright (C) 1986-2018 Xilinx, Inc. All Rights Reserved.
############################################################
set_directive_interface -mode ap_fifo -depth 5000 "parseEventsHW" dataStream
set_directive_interface -mode s_axilite -register "parseEventsHW"
set_directive_interface -mode m_axi -depth 5000 -offset slave -bundle gmem -num_read_outstanding 0 -max_read_burst_length 2 -max_write_burst_length 256 "parseEventsHW" eventSlice
