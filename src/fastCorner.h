#ifndef FASTCORNER
#define FASTCORNER

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"


#define DVS_WIDTH  240
// Change these two together
#define RESHAPE_FACTOR 16
#define DVS_HEIGHT RESHAPE_FACTOR*8


#define X_TYPE ap_uint<8>
#define Y_TYPE ap_uint<8>

// Change them together
#define TS_TYPE_BIT_WIDTH 32
#define LOG_TS_TYPE_BIT_WIDTH 5   // Log(TS_TYPE_BIT_WIDTH), used in pix read and pix write

#define col_pix_t ap_uint<RESHAPE_FACTOR * TS_TYPE_BIT_WIDTH>


#define INNER_SIZE 16
#define OUTER_SIZE 20

#define TEST_SORT_DATA_SIZE 20

#define MERGE_STAGES 5

#define SYMBOL_BITS 5

#endif
