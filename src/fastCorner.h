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

#define TEST_SORT_DATA_SIZE 10

#define MERGE_STAGES 5

#define SYMBOL_BITS 5

//#define SIZE 16
void mergeSortParallelWithSize(ap_uint<TS_TYPE_BIT_WIDTH> A[OUTER_SIZE], ap_uint<8> num_symbols,  ap_uint<TS_TYPE_BIT_WIDTH> B[OUTER_SIZE]);
void insertionCellSort(ap_uint<TS_TYPE_BIT_WIDTH> inData[20], ap_uint<TS_TYPE_BIT_WIDTH> outputData[20]);

template<int DATA_SIZE, int NPC>
void insertionSortParallel(ap_uint<TS_TYPE_BIT_WIDTH> A[DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> B[DATA_SIZE]);

void testSortHW(ap_uint<TS_TYPE_BIT_WIDTH> inputA[TEST_SORT_DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outputB[TEST_SORT_DATA_SIZE]);
void testSortedIdxData(ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> newIdx[OUTER_SIZE]);

void testRwSAEHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<5> *size);
void testFromTsDataCheckInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);
void testFromTsDataToIdxDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> idxData[OUTER_SIZE]);
void testFromTsDataToIdxInnerBoolDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<4> idxBoolData[OUTER_SIZE]);
void testFromTsDataToInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);

#endif
