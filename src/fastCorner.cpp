#include "fastCorner.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

// SAE (Surface of Active Event)
static TS_TYPE saeHW[1][DVS_WIDTH][DVS_HEIGHT];

const int innerCircleOffset[16][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
      {3, 0}, {3, -1}, {2, -2}, {1, -3},
      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
const int outerCircleOffset[20][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
      {4, 1}, {4, 0}, {4, -1}, {3, -2},
      {2, -3}, {1, -4}, {0, -4}, {-1, -4},
      {-2, -3}, {-3, -2}, {-4, -1}, {-4, 0},
      {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}};

void readCircleFromSAE(X_TYPE x, Y_TYPE y, TS_TYPE ts, TS_TYPE innerCircle[16])
{
	saeHW[0][x][y] = ts;
	for(int i = 0; i < 16; i++)
	{
		X_TYPE xOffset = innerCircleOffset[i][0];
		Y_TYPE yOffset = innerCircleOffset[i][1];
		innerCircle[i] = saeHW[0][x + xOffset][y + yOffset];
	}
}


template<int NPC>
void insertionSortParallel(TS_TYPE A[SIZE], TS_TYPE B[SIZE]) {
    #pragma HLS array_partition variable=B complete
    L1:  for(int i = 0; i < SIZE; i++)
    {
        #pragma HLS pipeline II=1
    	TS_TYPE item = A[i];
        L2:
        for (int j = SIZE - 1; j >= 0; j--)
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



void mergeArrays(TS_TYPE in[SIZE], int width, TS_TYPE out[SIZE]) {
  int f1 = 0;
  int f2 = width;
  int i2 = width;
  int i3 = 2*width;
  if(i2 >= SIZE) i2 = SIZE;
  if(i3 >= SIZE) i3 = SIZE;
 merge_arrays:
  for (int i = 0; i < SIZE; i++) {
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
      if(i2 >= SIZE) i2 = SIZE;
      if(i3 >= SIZE) i3 = SIZE;
      f2 = i2;
     }
  }
}

template<int DATA_SIZE>
void mergeSortParallel(TS_TYPE A[DATA_SIZE], TS_TYPE B[DATA_SIZE]) {
#pragma HLS DATAFLOW

    TS_TYPE temp[STAGES-1][DATA_SIZE];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    int width = 1;

    mergeArrays(A, width, temp[0]);
    width *= 2;

    stage:
	for (int stage = 1; stage < STAGES-1; stage++) {
#pragma HLS UNROLL
		mergeArrays(temp[stage-1], width, temp[stage]);
        width *= 2;
    }

	mergeArrays(temp[STAGES-2], width, B);
}


void fastCornerHW(X_TYPE x, Y_TYPE y, TS_TYPE ts, TS_TYPE B[SIZE])
{
    TS_TYPE A[SIZE];
	readCircleFromSAE(x, y, ts, A);

	mergeSortParallel<SIZE>(A, B);
	// insertionSortParallel<1>(A, B);
}
