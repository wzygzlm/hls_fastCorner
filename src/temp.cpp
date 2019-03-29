#include "fastCorner.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

void tempSorted(ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> newIdx[OUTER_SIZE])
{
#pragma HLS ARRAY_PARTITION variable=newIdx complete dim=0
#pragma HLS ARRAY_PARTITION variable=inData complete dim=0
	outerLoop:
	for(uint8_t i = 0; i < OUTER_SIZE; i++)
	{
#pragma HLS PIPELINE

		ap_uint<TS_TYPE_BIT_WIDTH> tmpData = inData[i];
		ap_uint<5> tempIdx = 0;

		innerLoop:
		for(uint8_t j = 0; j < OUTER_SIZE; j++)
		{
			if(inData[j] < tmpData)
			{
				tempIdx += 1;
			}

			if (size == INNER_SIZE && j == INNER_SIZE - 1)
			{
				newIdx[i] = tempIdx;
				return;
			}
		}

		newIdx[i] = tempIdx;
	}
}
