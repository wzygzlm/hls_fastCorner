#ifndef FASTCORNER
#define FASTCORNER

#include<stdint.h>
#include "ap_int.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define DEBUG 0

// This two macros will not influence the slice memory size,
// only used to indicate the real chip size.
// Used for removing boarder events.
#define DVS_REAL_WIDTH 346
#define DVS_REAL_HEIGHT 260

// If want to save more sources, this value could be chosen as some value which is easy to be calcuated
// with sub and shift register. For example we could choose 448 here because 448 = 512-64.
// So if some number x multiple by 448, if could be calucated by x*512 - x*64, and it means that
// we can use sub and shift register to get the multiple result.
// The reason we use 454 here is for memory utilization, because the maximum allowed number here is 455.
#define DVS_WIDTH  454

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
#define RESHAPE_FACTOR 32
#define DVS_HEIGHT RESHAPE_FACTOR*9

#define X_TYPE ap_uint<10>
#define Y_TYPE ap_uint<10>

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

typedef struct
{
	uint64_t inEventsNum;
	uint64_t outEventsNum;
	uint64_t cornerEventsNum;
} status_t;

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
void EVFastCornerStreamNoAxiLite(hls::stream< ap_uint<16> > &xStreamIn, hls::stream< ap_uint<16> > &yStreamIn, hls::stream< ap_uint<64> > &tsStreamIn, hls::stream< ap_uint<1> > &polStreamIn,
		hls::stream< ap_uint<16> > &xStreamOut, hls::stream< ap_uint<16> > &yStreamOut, hls::stream< ap_uint<64> > &tsStreamOut, hls::stream< ap_uint<1> > &polStreamOut,
		hls::stream< ap_uint<10> > &custDataStreamOut,
		ap_uint<32> config, status_t *status);
#endif
