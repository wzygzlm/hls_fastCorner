#include "fastCorner.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

// SAE (Surface of Active Event)
static TS_TYPE saeHW[1][DVS_WIDTH][DVS_HEIGHT];

const int innerCircleOffset[INNER_SIZE][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
      {3, 0}, {3, -1}, {2, -2}, {1, -3},
      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
const int outerCircleOffset[OUTER_SIZE][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
      {4, 1}, {4, 0}, {4, -1}, {3, -2},
      {2, -3}, {1, -4}, {0, -4}, {-1, -4},
      {-2, -3}, {-3, -2}, {-4, -1}, {-4, 0},
      {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}};

void readCircleFromSAE(X_TYPE x, Y_TYPE y, TS_TYPE ts, TS_TYPE innerCircle[INNER_SIZE])
{
	saeHW[0][x][y] = ts;
	for(int i = 0; i < INNER_SIZE; i++)
	{
		X_TYPE xOffset = innerCircleOffset[i][0];
		Y_TYPE yOffset = innerCircleOffset[i][1];
		innerCircle[i] = saeHW[0][x + xOffset][y + yOffset];
	}
}


template<int DATA_SIZE, int NPC>
void insertionSortParallel(TS_TYPE A[DATA_SIZE], TS_TYPE B[DATA_SIZE]) {
    #pragma HLS array_partition variable=B complete
    L1:  for(int i = 0; i < DATA_SIZE; i++)
    {
        #pragma HLS pipeline II=1
    	TS_TYPE item = A[i];
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
void mergeArrays(TS_TYPE in[DATA_SIZE], int width, TS_TYPE out[DATA_SIZE]) {
  int f1 = 0;
  int f2 = width;
  int i2 = width;
  int i3 = 2*width;
  if(i2 >= DATA_SIZE) i2 = DATA_SIZE;
  if(i3 >= DATA_SIZE) i3 = DATA_SIZE;
 merge_arrays:
  for (int i = 0; i < DATA_SIZE; i++) {
#pragma HLS PIPELINE II=1
      TS_TYPE t1 = in[f1];
      TS_TYPE t2 = in[f2];
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
void mergeSortParallel(TS_TYPE A[DATA_SIZE], TS_TYPE B[DATA_SIZE]) {
#pragma HLS DATAFLOW

    TS_TYPE temp[STAGES-1][DATA_SIZE];
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
    /* input */ TS_TYPE in[DATA_SIZE],
    /* output */ TS_TYPE out[DATA_SIZE]) {
	TS_TYPE previous_sorting[DATA_SIZE], sorting[DATA_SIZE];
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


void testSortHW(TS_TYPE inputA[TEST_SORT_DATA_SIZE], TS_TYPE outputB[TEST_SORT_DATA_SIZE])
{
//	 mergeSortParallel<TEST_SORT_DATA_SIZE, MERGE_STAGES> (inputA, outputB);
	insertionSortParallel<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
//	radixSort<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
}

void fastCornerHW(X_TYPE x, Y_TYPE y, TS_TYPE ts, TS_TYPE B[INNER_SIZE])
{
    TS_TYPE A[INNER_SIZE];
	readCircleFromSAE(x, y, ts, A);

	// mergeSortParallel<INNER_SIZE, MERGE_STAGES>(A, B);
	insertionSortParallel<INNER_SIZE, 1>(A, B);
}
