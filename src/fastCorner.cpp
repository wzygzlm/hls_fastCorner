#include "fastCorner.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

// SAE (Surface of Active Event)
static col_pix_t saeHW[1][DVS_HEIGHT/RESHAPE_FACTOR][DVS_WIDTH];

//const int innerCircleOffset[INNER_SIZE][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
//      {3, 0}, {3, -1}, {2, -2}, {1, -3},
//      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
//      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
//const int outerCircleOffset[OUTER_SIZE][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
//      {4, 1}, {4, 0}, {4, -1}, {3, -2},
//      {2, -3}, {1, -4}, {0, -4}, {-1, -4},
//      {-2, -3}, {-3, -2}, {-4, -1}, {-4, 0},
//      {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}};
const int innerCircleOffset[INNER_SIZE * 2] = {0, 3, 1, 3, 2, 2, 3, 1,
      3, 0, 3, -1, 2, -2, 1, -3,
      0, -3, -1, -3, -2, -2, -3, -1,
      -3, 0, -3, 1, -2, 2, -1, 3};
const int outerCircleOffset[OUTER_SIZE * 2] = {0, 4, 1, 4, 2, 3, 3, 2,
      4, 1, 4, 0, 4, -1, 3, -2,
      2, -3, 1, -4, 0, -4, -1, -4,
      -2, -3, -3, -2, -4, -1, -4, 0,
      -4, 1, -3, 2, -2, 3, -1, 4};

const ap_int<128> innerTest =  ap_int<128>("031322303f2e1d0cfcee1d03132213223130312213", 16);
const ap_int<160> outerTest = ap_int<160>("041433241404132231404142332414041322314", 16);


ap_uint<TS_TYPE_BIT_WIDTH> readOneDataFromCol(col_pix_t colData, ap_uint<8> idx)
{
#pragma HLS INLINE
	ap_uint<TS_TYPE_BIT_WIDTH> retData;
	// Use bit selection plus for-loop to read multi-bits from a wider bit width value
	// rather than use range selection directly. The reason is that the latter will use
	// a lot of shift-register which will increase a lot of LUTs consumed.
	readWiderBitsLoop: for(int8_t yIndex = 0; yIndex < TS_TYPE_BIT_WIDTH; yIndex++)
	{
#pragma HLS UNROLL
		const int bitOffset = LOG_TS_TYPE_BIT_WIDTH;   // This value should be equal to log(TS_TYPE_BIT_WIDTH)
		ap_uint<8 + bitOffset> colIdx;
		// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
		colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<8 + bitOffset>(idx * TS_TYPE_BIT_WIDTH).range(8 + bitOffset - 1, bitOffset);
		colIdx.range(bitOffset - 1, 0) = ap_uint<bitOffset>(yIndex);

		retData[yIndex] = colData[colIdx];
//		retData[yIndex] = colData[ap_uint<TS_TYPE_BIT_WIDTH>_BIT_WIDTH*idx + yIndex];
	}
	return retData;
}

void writeOneDataToCol(col_pix_t *colData, ap_uint<8> idx, ap_uint<TS_TYPE_BIT_WIDTH> toWriteData)
{
#pragma HLS INLINE
	writeWiderBitsLoop: for(int8_t yIndex = 0; yIndex < TS_TYPE_BIT_WIDTH; yIndex++)
	{
#pragma HLS UNROLL
		const int bitOffset = LOG_TS_TYPE_BIT_WIDTH;   // This value should be equal to log(TS_TYPE_BIT_WIDTH)
		ap_uint<8 + bitOffset> colIdx;
		// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
		colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<8 + bitOffset>(idx * TS_TYPE_BIT_WIDTH).range(8 + bitOffset - 1, bitOffset);
		colIdx.range(bitOffset - 1, 0) = ap_uint<bitOffset>(yIndex);

		(*colData)[colIdx] = toWriteData[yIndex];
	}
}


// create a function with an II=3
// - there are 3 registers explicitly used in the loop
// - we limit the multiplier instances allowed to 1
// -> so the tool can't schedule anything in parallel, so operations have to execute in serial fashion, II=3 is at least needed (or more depending on other clock constraints)
void my_func(int b[4], int &r) {
#pragma HLS inline off
#pragma HLS ALLOCATION instances=mul limit=1 operation
    int t=b[0],i;
mul_loop:
    for(i=1;i<4;i++) {
        t=t*b[i];
    }
    r=t;
}


// this is the top level of this short example
void top( hls::stream<int> &stream_input, hls::stream<int> &stream_output) {

    int i,buff[4];

    ap_uint<2> buff_index=0;
loop_a:
    for (i=0; i<16; i++) {
#pragma HLS pipeline II=1
        buff[buff_index]=stream_input.read();
// *** this is the place whereple uses the directive in TCL ***
		if (buff_index==3) {
			// this is executed every 4 cycles.
			my_occurrence_region:
			{
#pragma HLS OCCURRENCE cycle=4
				buff_index=0;
				int tmp;
				my_func(buff,tmp);
				stream_output.write(tmp);
			} // my_occurrence_region

		} else {
			buff_index++;
		}

    } // for loop
} // top function

void updateSAE(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts)
{
#pragma HLS INLINE
	col_pix_t tmpData;
	Y_TYPE yNewIdx = y%RESHAPE_FACTOR;

	tmpData = saeHW[0][y/RESHAPE_FACTOR][x];

	writeOneDataToCol(&tmpData, yNewIdx, ts);

	saeHW[0][y/RESHAPE_FACTOR][x] = tmpData;
}

void rwSAE(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<TS_TYPE_BIT_WIDTH> innerCircle[INNER_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outerCircle[OUTER_SIZE])
{
	readInnerCircleFromSAE:for(int8_t i = 0; i < INNER_SIZE + 1; i++)
	{
#pragma HLS DEPENDENCE variable=saeHW inter false
#pragma HLS PIPELINE rewind
		if (i == INNER_SIZE)
		{
			updateSAE(x, y, ts);
		}
		else
		{
			ap_uint<8> xInnerTest, xOuterTest;
			rwSAEReadOffsetBitsLoop:
			for (ap_uint<8> j = 0; j < 8; j++)
			{
#pragma HLS UNROLL
				xInnerTest[j] = innerTest[( (i << 3) , j(3,0) )];
				xOuterTest[j] = outerTest[( (i << 3) , j(3,0) )];
			}

			X_TYPE xInner = x + xInnerTest(3, 0);
			Y_TYPE yInner = y + xInnerTest(7, 4);
			Y_TYPE yInnerNewIdx = yInner%RESHAPE_FACTOR;

			X_TYPE xOuter = x + xOuterTest(3, 0);
			Y_TYPE yOuter = y + xOuterTest(7, 4);
			Y_TYPE yOuterNewIdx = yOuter%RESHAPE_FACTOR;

			innerCircle[i] = readOneDataFromCol(saeHW[0][yInner/RESHAPE_FACTOR][xInner], yInnerNewIdx);


			outerCircle[i] = readOneDataFromCol(saeHW[0][yOuter/RESHAPE_FACTOR][xOuter], yOuterNewIdx);
		}
	}
}


template<int DATA_SIZE, int NPC>
void insertionSortParallel(ap_uint<TS_TYPE_BIT_WIDTH> A[DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> B[DATA_SIZE]) {
    #pragma HLS array_partition variable=B complete
    L1:  for(int i = 0; i < DATA_SIZE; i++)
    {
        #pragma HLS pipeline II=1
    	ap_uint<TS_TYPE_BIT_WIDTH> item = A[i];
        L2:
        for (int j = DATA_SIZE - 1; j >= 0; j--)
        {
        	int t;
        	if(j > i)
        	{
        		t = B[j];
        	}
        	else if(j > 0 && B[j - 1] > item)
        	{
        		t = B[j - 1];
        	}
        	else
        	{
        		t = item;
        		if (j > 0)
        		{
        			item = B[j - 1];
        		}
        	}
        	B[j] = t;
        }
    }
}


template<int DATA_SIZE>
void mergeArrays(ap_uint<TS_TYPE_BIT_WIDTH> in[DATA_SIZE], int width, ap_uint<TS_TYPE_BIT_WIDTH> out[DATA_SIZE]) {
  int f1 = 0;
  int f2 = width;
  int i2 = width;
  int i3 = 2*width;
  if(i2 >= DATA_SIZE) i2 = DATA_SIZE;
  if(i3 >= DATA_SIZE) i3 = DATA_SIZE;
 merge_arrays:
  for (int i = 0; i < DATA_SIZE; i++) {
#pragma HLS PIPELINE II=1
      ap_uint<TS_TYPE_BIT_WIDTH> t1 = in[f1];
      ap_uint<TS_TYPE_BIT_WIDTH> t2 = in[f2];
    if((f1 < i2 && t1 <= t2) || f2 == i3) {
      out[i] = t1;
      f1++;
    } else {
      assert(f2 < i3);
      out[i] = t2;
      f2++;
    }
    if(f1 == i2 && f2 == i3) {
      f1 = i3;
      i2 += 2*width;
      i3 += 2*width;
      if(i2 >= DATA_SIZE) i2 = DATA_SIZE;
      if(i3 >= DATA_SIZE) i3 = DATA_SIZE;
      f2 = i2;
     }
  }
}

template<int DATA_SIZE, int STAGES>
void mergeSortParallel(ap_uint<TS_TYPE_BIT_WIDTH> A[DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> B[DATA_SIZE]) {
#pragma HLS DATAFLOW

    ap_uint<TS_TYPE_BIT_WIDTH> temp[STAGES-1][DATA_SIZE];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    int width = 1;

    mergeArrays<DATA_SIZE>(A, width, temp[0]);
    width *= 2;

    stage:
	for (int stage = 1; stage < STAGES-1; stage++) {
#pragma HLS UNROLL
		mergeArrays<DATA_SIZE>(temp[stage-1], width, temp[stage]);
        width *= 2;
    }

	mergeArrays<DATA_SIZE>(temp[STAGES-2], width, B);
}

const unsigned int RADIX = 16;
const unsigned int BITS_PER_LOOP = 4; // should be log2(RADIX)
typedef ap_uint<BITS_PER_LOOP> Digit;
template<int DATA_SIZE, int NPC>
void radixSort(
    /* input */ ap_uint<TS_TYPE_BIT_WIDTH> in[DATA_SIZE],
    /* output */ ap_uint<TS_TYPE_BIT_WIDTH> out[DATA_SIZE]) {
	ap_uint<TS_TYPE_BIT_WIDTH> previous_sorting[DATA_SIZE], sorting[DATA_SIZE];
    ap_uint<SYMBOL_BITS> digit_histogram[RADIX], digit_location[RADIX];
#pragma HLS ARRAY_PARTITION variable=digit_location complete dim=1
#pragma HLS ARRAY_PARTITION variable=digit_histogram complete dim=1
    Digit current_digit[DATA_SIZE];

 copy_in_to_sorting:
    for(int j = 0; j < DATA_SIZE; j++) {
#pragma HLS PIPELINE II=1
        sorting[j] = in[j];
    }

 radix_sort:
    for(int shift = 0; shift < 20; shift += BITS_PER_LOOP) {
    init_histogram:
        for(int i = 0; i < RADIX; i++) {
#pragma HLS pipeline II=1
            digit_histogram[i] = 0;
        }

    compute_histogram:
        for(int j = 0; j < DATA_SIZE; j++) {
#pragma HLS PIPELINE II=1
            Digit digit = (sorting[j] >> shift) & (RADIX - 1); // Extract a digit
            current_digit[j] = digit;  // Store the current digit for each symbol
            digit_histogram[digit]++;
            previous_sorting[j] = sorting[j]; // Save the current sorted order of symbols
        }

        digit_location[0] = 0;
    find_digit_location:
        for(int i = 1; i < RADIX; i++)
#pragma HLS PIPELINE II=1
            digit_location[i] = digit_location[i-1] + digit_histogram[i-1];

    re_sort:
        for(int j = 0; j < DATA_SIZE; j++) {
#pragma HLS PIPELINE II=1
            Digit digit = current_digit[j];
            sorting[digit_location[digit]] = previous_sorting[j]; // Move symbol to new sorted location
            out[digit_location[digit]] = previous_sorting[j]; // Also copy to output
            digit_location[digit]++; // Update digit_location
        }
    }
}


void testSortHW(ap_uint<TS_TYPE_BIT_WIDTH> inputA[TEST_SORT_DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outputB[TEST_SORT_DATA_SIZE])
{
	 mergeSortParallel<TEST_SORT_DATA_SIZE, MERGE_STAGES> (inputA, outputB);
//	insertionSortParallel<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
//	radixSort<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
}

void fastCornerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<TS_TYPE_BIT_WIDTH> B[INNER_SIZE])
{
    ap_uint<TS_TYPE_BIT_WIDTH> inner[INNER_SIZE];
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    rwSAE(x, y, ts, inner, outer);

	 mergeSortParallel<INNER_SIZE, MERGE_STAGES>(inner, B);
//	insertionSortParallel<INNER_SIZE, 1>(A, B);
}
