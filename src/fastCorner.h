#ifndef FASTCORNER
#define FASTCORNER

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define DVS_WIDTH  256
#define DVS_HEIGHT 18*8

#define RESHAPE_FACTOR 18

#define X_TYPE ap_uint<8>
#define Y_TYPE ap_uint<8>

// These two should change together
#define TS_TYPE ap_uint<32>
#define TS_TYPE_BIT_WIDTH 32

#define col_pix_t ap_uint<RESHAPE_FACTOR * TS_TYPE_BIT_WIDTH>


#define INNER_SIZE 16
#define OUTER_SIZE 20

#define TEST_SORT_DATA_SIZE 32

#define MERGE_STAGES 5

#define SYMBOL_BITS 5

#endif
