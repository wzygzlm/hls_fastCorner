#ifndef FASTCORNER
#define FASTCORNER

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define DVS_WIDTH  240
#define DVS_HEIGHT 5


#define X_TYPE ap_uint<8>
#define Y_TYPE ap_uint<8>
#define TS_TYPE ap_uint<180>

#define INNER_SIZE 16
#define OUTER_SIZE 20

#define TEST_SORT_DATA_SIZE 32

#define MERGE_STAGES 5

#define SYMBOL_BITS 5

#endif
