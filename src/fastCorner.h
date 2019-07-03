#ifndef FASTCORNER
#define FASTCORNER

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define DEBUG 0

#define DVS_WIDTH  204

#define POLARITY_SHIFT 1
#define POLARITY_MASK 0x00000001
#define POLARITY_Y_ADDR_SHIFT 2
#define POLARITY_Y_ADDR_MASK 0x000001FF      //  Reduce mask bit width to reduce LUTs
#define POLARITY_X_ADDR_SHIFT 17
#define POLARITY_X_ADDR_MASK 0x000001FF      //  Reduce mask bit width to reduce LUTs

typedef ap_uint<49> apUint49_t;
typedef ap_uint<17> apUint17_t;
typedef ap_uint<15> apUint15_t;
typedef ap_uint<6> apUint6_t;
typedef ap_uint<1> apUint1_t;

// Change these two together
#define RESHAPE_FACTOR 16
#define DVS_HEIGHT RESHAPE_FACTOR*10

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

#define LOOPS_PER_EVENT 2   // To process one event, how many loops are required

//#define SIZE 16
void mergeSortParallelWithSize(ap_uint<TS_TYPE_BIT_WIDTH> A[OUTER_SIZE], ap_uint<8> num_symbols,  ap_uint<TS_TYPE_BIT_WIDTH> B[OUTER_SIZE]);
void insertionCellSort(ap_uint<TS_TYPE_BIT_WIDTH> inData[20], ap_uint<TS_TYPE_BIT_WIDTH> outputData[20]);

template<int DATA_SIZE, int NPC>
void insertionSortParallel(ap_uint<TS_TYPE_BIT_WIDTH> A[DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> B[DATA_SIZE]);

void testSortHW(ap_uint<TS_TYPE_BIT_WIDTH> inputA[TEST_SORT_DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outputB[TEST_SORT_DATA_SIZE]);
void testSortedIdxData(ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> newIdx[OUTER_SIZE]);

void testRwSAEHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<5> *size);
void testIdxDataToIdxInnerBoolDataHW(ap_uint<5> newIdx[OUTER_SIZE], ap_uint<5> size, ap_uint<4> condFlg[INNER_SIZE]);
void testFromTsDataCheckInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);
void testFromTsDataCheckOuterCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);
void testFromTsDataToIdxDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> idxData[OUTER_SIZE]);
void testFromTsDataToIdxInnerBoolDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<4> idxBoolData[OUTER_SIZE]);
void testFromTsDataToInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);
void testFromTsDataToOuterCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner);
void fastCornerInnerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<1> *isCorner);
void fastCornerOuterHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<1> *isCorner);
void fastCornerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<1> *isCorner);
void parseEventsHW(uint64_t * dataStream, int32_t eventsArraySize, uint64_t *eventSlice);

#endif
