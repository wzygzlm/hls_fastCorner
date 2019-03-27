#ifndef INSERTION_CELL_SORT
#define INSERTION_CELL_SORT

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define DTYPE ap_uint<32>
#define SIZE 8
#define STAGES 4
#define SYMBOL_BITS 5

// template<int NPC>
void insertion_cell_sort(hls::stream<DTYPE> &in, hls::stream<DTYPE> &out);

#endif
