#include "fastCorner.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

// SAE (Surface of Active Event)
// Some important things about the hardware memory overflow.
// Here for example RESHAPE_FACTOR = 16 is the power of 2.
// Then in this memory saeHW, the address will be 11bits. How can we know that it is 11bits?
// If we represent it as saeHW[0][y][x].
// The address is calculated from y/RESHAPE_FACTOR * (256 - 16) + x.
// y >> 4 will have 3 bits and x have 8 bits. The maximum of y is assumed 128.
// So the real capacity of the memory are 128/16 * 240 = 1920 words.
// Howvere, we have 11bits memory, so it can represent a memory capacity of 2^11 = 2048 words.
// That means, if we try to read words range from 1920 to 2048, it will return XXX;
// But if we try to read words range which is bigger than 2048, it will go back to 2048.
// For example, if y = 128, x = 127, then 128/16 * 240 +  127 = 2047, it will return XXX
// But if y = 128, x = 128, then 128/16 * 128 + 128 = 2048, it will return of value of the 0th word.
// It has observed from the simulation.
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

//const ap_int<128> innerTest =  ap_int<128>("03132231303f2e1d0dfdeedfd0d1e2f3", 16);
//const ap_int<160> outerTest = ap_int<160>("0414233241404f3e2d1c0cfceddecfc0c1d2e3f4", 16);

const ap_int<128> innerTest =  ap_int<128>("3f2e1d0dfdeedfd0d1e2f30313223130", 16);
const ap_int<160> outerTest = ap_int<160>("4f3e2d1c0cfceddecfc0c1d2e3f4041423324140", 16);

// Function Description: return the minimum value of an array.
template<typename DATA_TYPE, int DATA_SIZE>
DATA_TYPE min(DATA_TYPE inArr[DATA_SIZE], int8_t *index)
{
//#pragma HLS PIPELINE
//#pragma HLS ARRAY_RESHAPE variable=inArr complete dim=1
#pragma HLS INLINE off
	DATA_TYPE tmp = inArr[0];
	int8_t tmpIdx = 0;
	minLoop: for(int8_t i = 0; i < DATA_SIZE; i++)
	{
		// Here is a bug. Use the if-else statement,
		// cannot use the question mark statement.
		// Otherwise a lot of muxs will be generated,
		// DON'T KNOW WHY. SHOULD BE A BUG.
		if(inArr[i] < tmp) tmpIdx = i;
		if(inArr[i] < tmp) tmp = inArr[i];
//		tmp = (inArr[i] < tmp) ? inArr[i] : tmp;
	}
	*index = tmpIdx;
	return tmp;
}

// Function Description: convert ap_memory to several ap_none ports
template<int NPC>
void convertInterface(ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> size, hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > &inStream)
{
#pragma HLS ARRAY_PARTITION variable=inData cyclic factor=NPC dim=0
#pragma HLS INLINE off

	ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> tmpData = 0;

	for(uint8_t i = 0; i < OUTER_SIZE; i = i + NPC)
	{
#pragma HLS LOOP_TRIPCOUNT min=0 max=20/NPC
#pragma HLS PIPELINE rewind
		if (i >= size)
		{
			break;
		}
		for(uint8_t j = 0; j < NPC; j++)
		{
			for (uint8_t yIndex = 0; yIndex < TS_TYPE_BIT_WIDTH; yIndex++)
			{
#pragma HLS UNROLL
				const int bitOffset = LOG_TS_TYPE_BIT_WIDTH;   // This value should be equal to log(TS_TYPE_BIT_WIDTH)
				ap_uint<8 + bitOffset> colIdx;
				// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
				colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<8 + bitOffset>((i + j) * TS_TYPE_BIT_WIDTH).range(8 + bitOffset - 1, bitOffset);
				colIdx.range(bitOffset - 1, 0) = ap_uint<bitOffset>(yIndex);

				tmpData[colIdx] = inData[i + j][yIndex];
			}
		}
	}
	inStream.write(tmpData);
}

// Function Description: return the idx of current data in the sorted data array.
void idxSorted(ap_uint<TS_TYPE_BIT_WIDTH> oriData, ap_uint<TS_TYPE_BIT_WIDTH> tsData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> *newIdx)
{
#pragma HLS ARRAY_PARTITION variable=tsData complete dim=0
#pragma HLS PIPELINE
#pragma HLS INLINE
	ap_uint<5> temp = 0;
	for(uint8_t i = 0; i < OUTER_SIZE; i++ )
	{
		ap_uint<1> cond1 = (tsData[i] < oriData);  // Notice the difference between < and <= here.
		temp += cond1;
//		if (size == INNER_SIZE && i == INNER_SIZE - 1)
//		{
//			*newIdx = temp;
//			return;
//		}
	}
	*newIdx = temp;
}

// Function Description: return the idx of current data in the sorted data array.
//void idxSorted(ap_uint<TS_TYPE_BIT_WIDTH> oriData, ap_uint<TS_TYPE_BIT_WIDTH> tsData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> *newIdx)
//{
//#pragma HLS ARRAY_PARTITION variable=tsData cyclic factor=8 dim=0
//#pragma HLS INLINE
//	assert(size==16||size==20||size==0);
//	if(size == 0)
//	{
//		*newIdx = 0;
//		return;
//	}
//
//	ap_uint<5> temp = 0;
//	for(uint8_t i = 0; i < size; i = i + 16 )
//	{
//#pragma HLS PIPELINE
//#pragma HLS LOOP_TRIPCOUNT min=0 max=2
//		for(uint8_t j = 0; j < 16; j++)
//		{
//			ap_uint<1> cond1 = (tsData[i + j] < oriData);  // Notice the difference between < and <= here.
//			temp += cond1;
//		}
//	}
//	*newIdx = temp;
//}


void sortedTest(ap_uint<5> oriIdx, ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> *newIdx)
{
#pragma HLS ARRAY_PARTITION variable=inData complete dim=0
#pragma HLS PIPELINE
#pragma HLS INLINE
	ap_uint<5> temp = 0;
	for(uint8_t i = 0; i < OUTER_SIZE; i++ )
	{
		ap_uint<1> cond1 = (inData[i] < inData[oriIdx]);  // Notice the difference between < and <= here.
		temp += cond1;
		if (size == INNER_SIZE && i == INNER_SIZE - 1)
		{
			*newIdx = temp;
			return;
		}
	}
	*newIdx = temp;
}



// Function Description: convert the current data array to sorted idx array.
template<int NPC>
void sortedIdxData(ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> newIdx[OUTER_SIZE])
{
#pragma HLS INLINE
//#pragma HLS ARRAY_PARTITION variable=newIdx complete dim=0
	for(uint8_t i = 0; i < size; i = i + NPC)
	{
#pragma HLS LOOP_TRIPCOUNT min=0 max=20/NPC
#pragma HLS PIPELINE
		for(uint8_t j = 0; j < NPC; j++)
		{
//			ap_uint<5> tmpIdx;
//			idxSorted(inData[i + j], inData, size, &newIdx[i + j]);
			sortedTest(i + j, inData, size, &newIdx[i + j]);
//			newIdx[i + 0] = tmpIdx;
		}
	}
}


// Function Description: convert the current data array to sorted idx array.
template<int NPC>
void sortedIdxStream(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > &tsStream, ap_uint<5> size, ap_uint<5> newIdx[OUTER_SIZE])
{
assert(size <= OUTER_SIZE);
#pragma HLS INLINE off
	ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> tmpData = tsStream.read();
	ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE];

	for(uint8_t j = 0; j < OUTER_SIZE; j++)
	{
#pragma HLS UNROLL
		for (uint8_t yIndex = 0; yIndex < TS_TYPE_BIT_WIDTH; yIndex++)
		{
#pragma HLS UNROLL
			const int bitOffset = LOG_TS_TYPE_BIT_WIDTH;   // This value should be equal to log(TS_TYPE_BIT_WIDTH)
			ap_uint<8 + bitOffset> colIdx;
			// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
			colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<8 + bitOffset>(j * TS_TYPE_BIT_WIDTH).range(8 + bitOffset - 1, bitOffset);
			colIdx.range(bitOffset - 1, 0) = ap_uint<bitOffset>(yIndex);

			inData[j][yIndex] = tmpData[colIdx];
		}
	}

	for(uint8_t i = 0; i < OUTER_SIZE/NPC; i = i + 1)
	{
// #pragma HLS LOOP_TRIPCOUNT min=0 max=20/NPC
#pragma HLS PIPELINE rewind
		if (i * NPC >= size)
		{
			break;
		}
		for(uint8_t j = 0; j < NPC; j++)
		{
//			ap_uint<5> tmpIdx;
			idxSorted(inData[i * NPC + j], inData, size, &newIdx[i * NPC + j]);
//			newIdx[i + 0] = tmpIdx;
		}
	}
}

// Convert index data to the bool version data. It has two types, one for INNER circle and the other for OUTER circle.
template<int NPC>
void idxDataToIdxInnerBoolData(ap_uint<5> newIdx[OUTER_SIZE], ap_uint<5> size, ap_uint<4> condFlg[OUTER_SIZE])
{
assert(size<=OUTER_SIZE);
#pragma HLS ARRAY_PARTITION variable=condFlg cyclic factor=NPC dim=0
//#pragma HLS ARRAY_PARTITION variable=condFlg complete dim=0
#pragma HLS ARRAY_PARTITION variable=newIdx cyclic factor=NPC/2 dim=0
// #pragma HLS ARRAY_PARTITION variable=newIdx complete dim=0

	for(uint8_t i = 0; i <= OUTER_SIZE/NPC; i = i + 1)
	{
#pragma HLS PIPELINE
		InitRegion:
		{
//#pragma HLS LATENCY min=1
			if (i * NPC >= size)
			{
				break;
			}
		}
		for(uint8_t j = 0; j < NPC; j++)
		{
			uint8_t tmpIndex = i * NPC + j;
			ap_uint<5> tmpNewIdx = newIdx[tmpIndex];
			ap_uint<4> tmpTmp;
			// The condition should be the idxData > (INNER_SIZE -3).
			// However, in order to make the idxSorted could be shared by inner circle and outer circle together.
			// We use a method that compare "size" values to all the input data which has OUTER_SIZE values in total.
			// On the other hand, if the valid input data number is less than OUTER_SIZE, the other input data will be filled with 0.
			// Thus, all the idxData for inner circle value will be added 4 (OUTER_SIZE - INNER_SIZE = 20 - 16 =4)
			// When we check the innner idx data, we need to remove it.
			tmpTmp[0] = (tmpNewIdx  >= INNER_SIZE - 3 + OUTER_SIZE - INNER_SIZE);
			tmpTmp[1] = (tmpNewIdx  >= INNER_SIZE - 4 + OUTER_SIZE - INNER_SIZE);
			tmpTmp[2] = (tmpNewIdx  >= INNER_SIZE - 5 + OUTER_SIZE - INNER_SIZE);
			tmpTmp[3] = (tmpNewIdx  >= INNER_SIZE - 6 + OUTER_SIZE - INNER_SIZE);
			condFlg[tmpIndex] = tmpTmp;
		}
	}
}

template<int NPC>
void idxInnerBoolDataToCorner(ap_uint<4> condFlg[INNER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS ARRAY_PARTITION variable=condFlg cyclic factor=NPC dim=0

	ap_uint<1> isCornerTemp = 0;
	ap_uint<4> tempCond[NPC];
#pragma HLS ARRAY_PARTITION variable=tempCond complete dim=0

	ap_uint<4> cond[INNER_SIZE];
#pragma HLS ARRAY_PARTITION variable=cond complete dim=0

	for (uint8_t i = 0; i < INNER_SIZE; i = i + 2 * NPC)
	{
#pragma HLS PIPELINE
		for (uint8_t k = 0; k < 2 * NPC; k++)
		{
			cond[i + k] = condFlg[i + k];
		}
	}

	for(uint8_t i = 0; i <= OUTER_SIZE/NPC; i = i + 1)
	{
#pragma HLS PIPELINE
		InitRegion:
		{
//#pragma HLS LATENCY min=1
			if (i * NPC >= size)
			{
				break;
			}
		}

		for (uint8_t k = 0; k < NPC; k++)
		{
			tempCond[k] = ap_uint<4>(15);
			for (uint8_t n = 0; n < 4; n++)
			{
				for (uint8_t j = 0; j < 3 + n; j++)
				{
					ap_uint<5> tmpIdx = i * NPC + j + k;
					if (tmpIdx >= INNER_SIZE) tmpIdx = tmpIdx - INNER_SIZE;

					ap_uint<1> tmpTmp = tempCond[k][n];
					tmpTmp = tmpTmp & cond[tmpIdx][n];
					tempCond[k][n] = tmpTmp;
				}
				isCornerTemp |= tempCond[k][n];
			}
		}
	}
	*isCorner = isCornerTemp;
}

void testSortedIdxData(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > &tsStream, ap_uint<5> size, ap_uint<5> newIdx[OUTER_SIZE])
{
#pragma HLS ARRAY_PARTITION variable=newIdx complete dim=0

	ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> tmpData = tsStream.read();
	ap_uint<TS_TYPE_BIT_WIDTH> inData[OUTER_SIZE];


	for(uint8_t j = 0; j < OUTER_SIZE; j++)
	{
#pragma HLS UNROLL
		for (uint8_t yIndex = 0; yIndex < TS_TYPE_BIT_WIDTH; yIndex++)
		{
#pragma HLS UNROLL
			const int bitOffset = LOG_TS_TYPE_BIT_WIDTH;   // This value should be equal to log(TS_TYPE_BIT_WIDTH)
			ap_uint<8 + bitOffset> colIdx;
			// Concatenate and bit shift rather than multiple and accumulation (MAC) can save area.
			colIdx.range(8 + bitOffset - 1, bitOffset) = ap_uint<8 + bitOffset>(j * TS_TYPE_BIT_WIDTH).range(8 + bitOffset - 1, bitOffset);
			colIdx.range(bitOffset - 1, 0) = ap_uint<bitOffset>(yIndex);

			inData[j][yIndex] = tmpData[colIdx];
		}
	}

	sortedIdxData<2>(inData, size, newIdx);

//	sortedIdxStream<2>(tsStream, size, newIdx);
}

template<int NPC>
void checkInnerIdx(ap_uint<5> idxData[INNER_SIZE + 6 - 1], ap_uint<5> size, ap_uint<1> *isCorner)
{
	/* This is a good example to show the LUTs as a function of the NPC
	 * Decreasing factor doesn't mean decreasing the performance.
	 * In this example, change the partition from factor to completely, the LUTs will increase a lot.
	 * The start index is fixed (always the multiple of NPC), so the partition could increase a lot.
	 * Another very interesting thing here is that: if NPC equals to 1, then i become arbitrary again
	 * which is not in a fixed pattern. In this case, multiplxer will be generated again.
	 *
	 * It's also a good example to compare the parameter cyclic and block here.
	 *
	 * The final expression LUTs# is about:  (NPC+2)*icmp + (NPC*2)*and + (NPC+1)*or + 2*adder for only streak length=3.
	 * If multiple streaks are used, the formulation should be derived again.
	 *
	 * The max number of data read from the memory M should be less than 2*factor.
	 * At the same time, NPC should be >= factor. Otherwise select will be synthesed.
	 *
	 * For example, if we the max streak lenght is 6, then the max number of data read from the memory is M = 6 + NPC - 1 = NPC + 5
	 * If we set NPC=2, then factor can only be 1 or 2. Say we choose factor = 2, then M = 7 > 2 * 2. So this is not a good combination.
	 * Set NPC = 4 and factor = 4, M = 9 still bigger than 2 * factor = 8.
	 * So we should use NPC = 5 and factor = 5. The factor should be divided by the size (if size = 20).
	 *
	 * TODO: compare these cases: 1. the loop_count is the multiple of NPC 2. the loop count is not the multiple of NPC
	 * 						      3. decrease M to make II = 1 under NPC = 4 and NPC = 2 to check the resource usage reducing.
	 * */
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=idxData cyclic factor=NPC dim=0
	ap_uint<1> isCornerTemp = 0;
//	if (size == 0)
//	{
//		*isCorner = isCornerTemp;
//		return;
//	}
	for(uint8_t i = 0; i <= OUTER_SIZE/NPC; i = i + 1)   // This multipler form is easier to understand but consume more resources.
	{
#pragma HLS LOOP_TRIPCOUNT min=0 max=16/NPC
#pragma HLS PIPELINE
		InitRegion:
		{
//#pragma HLS LATENCY min=1
			if (i * NPC >= size)
			{
				break;
			}
		}
		ap_uint<1> cond[4][6 + NPC - 1];
		for (uint8_t m = 0; m < 3 + NPC - 1; m++)
		{
			// The condition should be the idxData > (INNER_SIZE -3).
			// However, in order to make the idxSorted could be shared by inner circle and outer circle together.
			// We use a method that compare "size" values to all the input data which has OUTER_SIZE values in total.
			// On the other hand, if the valid input data number is less than OUTER_SIZE, the other input data will be filled with 0.
			// Thus, all the idxData for inner circle value will be added 4 (OUTER_SIZE - INNER_SIZE = 20 - 16 = 4)
			// When we check the innner idx data, we need to remove it.
			cond[0][m] = (idxData[(i * NPC + m)%16] >= INNER_SIZE - 3 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond2[4 + NPC - 1];
		for (uint8_t m = 0; m < 4 + NPC - 1; m++)
		{
			cond[1][m] = (idxData[(i * NPC + m)%16] >= INNER_SIZE - 4 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond3[5 + NPC - 1];
		for (uint8_t m = 0; m < 5 + NPC - 1; m++)
		{
			cond[2][m] = (idxData[(i * NPC + m)%16] >= INNER_SIZE - 5 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond4[6 + NPC - 1];
		for (uint8_t m = 0; m < 6 + NPC - 1; m++)
		{
			cond[3][m] = (idxData[(i * NPC + m)%16] >= INNER_SIZE - 6 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> tempCond[4][NPC];

		for (uint8_t k = 0; k < NPC; k++)
		{
			for (uint8_t n = 0; n < 4; n++)
			{
				tempCond[n][k] = 1;
				for (uint8_t j = 0; j < 3 + n; j++)
				{
					tempCond[n][k] &= cond[n][j + k];
				}
				isCornerTemp |= tempCond[n][k];

//				if (isCornerTemp == 1)
//				{
//					*isCorner = isCornerTemp ;
//					std::cout << "HW: Position is :" << (int)(i + k) << " and streak size is: " << (int)(n + 3) << std::endl;
//					return;
//				}
			}
		}
		*isCorner = isCornerTemp ;
	}
}


template<int NPC>
void checkInnerIdxV2(ap_uint<5> idxData[INNER_SIZE + 6 - 1], ap_uint<5> size, ap_uint<1> *isCorner)
{
	/* This is a good example to show the LUTs as a function of the NPC
	 * Decreasing factor doesn't mean decreasing the performance.
	 * In this example, change the partition from factor to completely, the LUTs will increase a lot.
	 * The start index is fixed (always the multiple of NPC), so the partition could increase a lot.
	 * Another very interesting thing here is that: if NPC equals to 1, then i become arbitrary again
	 * which is not in a fixed pattern. In this case, multiplxer will be generated again.
	 *
	 * It's also a good example to compare the parameter cyclic and block here.
	 *
	 * The final expression LUTs# is about:  (NPC+2)*icmp + (NPC*2)*and + (NPC+1)*or + 2*adder for only streak length=3.
	 * If multiple streaks are used, the formulation should be derived again.
	 *
	 * The max number of data read from the memory M should be less than 2*factor.
	 * At the same time, NPC should be >= factor. Otherwise select will be synthesed.
	 *
	 * For example, if we the max streak lenght is 6, then the max number of data read from the memory is M = 6 + NPC - 1 = NPC + 5
	 * If we set NPC=2, then factor can only be 1 or 2. Say we choose factor = 2, then M = 7 > 2 * 2. So this is not a good combination.
	 * Set NPC = 4 and factor = 4, M = 9 still bigger than 2 * factor = 8.
	 * So we should use NPC = 5 and factor = 5. The factor should be divided by the size (if size = 20).
	 *
	 * TODO: compare these cases: 1. the loop_count is the multiple of NPC 2. the loop count is not the multiple of NPC
	 * 						      3. decrease M to make II = 1 under NPC = 4 and NPC = 2 to check the resource usage reducing.
	 * */
#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=idxData cyclic factor=NPC dim=0
	ap_uint<1> isCornerTemp = 0;
//	if (size == 0)
//	{
//		*isCorner = isCornerTemp;
//		return;
//	}
	for(uint8_t i = 0; i < OUTER_SIZE; i = i + NPC)
	{
#pragma HLS LOOP_TRIPCOUNT min=0 max=16/NPC
#pragma HLS PIPELINE
		InitRegion:
		{
//#pragma HLS LATENCY min=1
			if (i >= size)
			{
				break;
			}
		}
		ap_uint<1> cond[4][6 + NPC - 1];
		for (uint8_t m = 0; m < 3 + NPC - 1; m++)
		{
			// The condition should be the idxData > (INNER_SIZE -3).
			// However, in order to make the idxSorted could be shared by inner circle and outer circle together.
			// We use a method that compare "size" values to all the input data which has OUTER_SIZE values in total.
			// On the other hand, if the valid input data number is less than OUTER_SIZE, the other input data will be filled with 0.
			// Thus, all the idxData for inner circle value will be added 4 (OUTER_SIZE - INNER_SIZE = 20 - 16 =4)
			// When we check the innner idx data, we need to remove it.
			cond[0][m] = (idxData[(i + m)%16] >= INNER_SIZE - 3 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond2[4 + NPC - 1];
		for (uint8_t m = 0; m < 4 + NPC - 1; m++)
		{
			cond[1][m] = (idxData[(i + m)%16] >= INNER_SIZE - 4 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond3[5 + NPC - 1];
		for (uint8_t m = 0; m < 5 + NPC - 1; m++)
		{
			cond[2][m] = (idxData[(i + m)%16] >= INNER_SIZE - 5 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> cond4[6 + NPC - 1];
		for (uint8_t m = 0; m < 6 + NPC - 1; m++)
		{
			cond[3][m] = (idxData[(i + m)%16] >= INNER_SIZE - 6 + OUTER_SIZE - INNER_SIZE);
		}

		ap_uint<1> tempCond[4][NPC];

		for (uint8_t k = 0; k < NPC; k++)
		{
			for (uint8_t n = 0; n < 4; n++)
			{
				tempCond[n][k] = 1;
				for (uint8_t j = 0; j < 3 + n; j++)
				{
					tempCond[n][k] &= cond[n][j + k];
				}
				isCornerTemp |= tempCond[n][k];

//				if (isCornerTemp == 1)
//				{
//					*isCorner = isCornerTemp ;
//					std::cout << "HW: Position is :" << (int)(i + k) << " and streak size is: " << (int)(n + 3) << std::endl;
//					return;
//				}
			}
		}
		*isCorner = isCornerTemp ;
	}
}


template<int NPC>
void checkOuterIdx(ap_uint<5> idxData[OUTER_SIZE + 8 - 1], ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=idxData cyclic factor=NPC dim=0

	ap_uint<1> isCornerTemp = 0;
	for(uint8_t i = 0; i < OUTER_SIZE; i = i + NPC)
	{
#pragma HLS LOOP_TRIPCOUNT min=0 max=20/NPC
#pragma HLS PIPELINE
		InitRegion:
		{
//#pragma HLS LATENCY min=1
			if (i >= size)
			{
				break;
			}
		}
		ap_uint<1> cond[5][8 + NPC - 1];
		for (uint8_t m = 0; m < 4 + NPC - 1; m++)
		{
			ap_uint<5> tmpIdx = i + m;
			if (tmpIdx >= OUTER_SIZE) tmpIdx = tmpIdx - OUTER_SIZE;

			cond[0][m] = (idxData[tmpIdx] >= OUTER_SIZE - 4);
		}

		ap_uint<1> cond2[5 + NPC - 1];
		for (uint8_t m = 0; m < 5 + NPC - 1; m++)
		{
			ap_uint<5> tmpIdx = i + m;
			if (tmpIdx >= OUTER_SIZE) tmpIdx = tmpIdx - OUTER_SIZE;

			cond[1][m] = (idxData[tmpIdx] >= OUTER_SIZE - 5);
		}

		ap_uint<1> cond3[6 + NPC - 1];
		for (uint8_t m = 0; m < 6 + NPC - 1; m++)
		{
			ap_uint<5> tmpIdx = i + m;
			if (tmpIdx >= OUTER_SIZE) tmpIdx = tmpIdx - OUTER_SIZE;

			cond[2][m] = (idxData[tmpIdx] >= OUTER_SIZE - 6);
		}

		ap_uint<1> cond4[7 + NPC - 1];
		for (uint8_t m = 0; m < 7 + NPC - 1; m++)
		{
			ap_uint<5> tmpIdx = i + m;
			if (tmpIdx >= OUTER_SIZE) tmpIdx = tmpIdx - OUTER_SIZE;

			cond[3][m] = (idxData[tmpIdx] >= OUTER_SIZE - 7);
		}

		ap_uint<1> cond5[8 + NPC - 1];
		for (uint8_t m = 0; m < 8 + NPC - 1; m++)
		{
			ap_uint<5> tmpIdx = i + m;
			if (tmpIdx >= OUTER_SIZE) tmpIdx = tmpIdx - OUTER_SIZE;

			cond[4][m] = (idxData[tmpIdx] >= OUTER_SIZE - 8);
		}

		ap_uint<1> tempCond[5][NPC];

		for (uint8_t k = 0; k < NPC; k++)
		{
			for (uint8_t n = 0; n < 5; n++)
			{
				tempCond[n][k] = 1;
				for (uint8_t j = 0; j < 4 + n; j++)
				{
					tempCond[n][k] &= cond[n][j + k];
				}
				isCornerTemp |= tempCond[n][k];

//				if (isCornerTemp == 1)
//				{
//					*isCorner = isCornerTemp ;
//					std::cout << "HW: Position is :" << (int)(i + k) << " and streak size is: " << (int)(n + 4) << std::endl;
//					return;
//				}

			}
		}

		*isCorner = isCornerTemp ;
	}
}

void testCheckInnerIdx(ap_uint<5> idxData[INNER_SIZE + 6 - 1], ap_uint<5> size, ap_uint<1> *isCorner)
{
//	checkInnerIdx<5>(idxData, size, isCorner);   // If resource is not enough, decrease this number to increase II a little.
	checkInnerIdxV2<4>(idxData, size, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void testCheckOuterIdx(ap_uint<5> idxData[OUTER_SIZE + 8 - 1], ap_uint<5> size, ap_uint<1> *isCorner)
{
	checkOuterIdx<5>(idxData, size, isCorner);    // NPC = 7 could make II = 1 but we might not need so fast.
}

void checkIdx(ap_uint<5> inData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner)
{
	if(size == INNER_SIZE)
	{
		ap_uint<5> idxData[INNER_SIZE + 6 - 1];
		for (uint8_t i = INNER_SIZE; i < INNER_SIZE + 6 - 1; i++)
		{
			idxData[i] = inData[i - INNER_SIZE];
		}
		checkInnerIdx<5>(idxData, size, isCorner);
	}
	else if(size == 20)
	{
		ap_uint<5> idxData[OUTER_SIZE + 8 - 1];
		for (uint8_t i = OUTER_SIZE; i < OUTER_SIZE + 8 - 1; i++)
		{
			idxData[i] = inData[i - OUTER_SIZE];
		}
		checkOuterIdx<5>(idxData, size, isCorner);
	}
	else
	{
		*isCorner = 0;
	}
}

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

void getXandY(const uint64_t * data, X_TYPE *x, Y_TYPE *y, ap_uint<TS_TYPE_BIT_WIDTH> *ts, hls::stream<apUint17_t> &packetEventDataStream)
{
#pragma HLS PIPELINE
	uint64_t tmp = *data;
	*x = ((tmp) >> POLARITY_X_ADDR_SHIFT) & POLARITY_X_ADDR_MASK;
	*y = ((tmp) >> POLARITY_Y_ADDR_SHIFT) & POLARITY_Y_ADDR_MASK;
	bool pol  = ((tmp) >> POLARITY_SHIFT) & POLARITY_MASK;
	*ts = tmp >> 32;
}


template<int READ_NPC>   //  Due to the memory has 2 ports at most for arbitrary reading, here this number could be only 1 or 2.
void rwSAE(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<5> *size)
{
//#pragma HLS RESOURCE variable=saeHW core=RAM_T2P_BRAM
#pragma HLS INLINE off
	if(stage == 0)
	{
//		updateSAE(x, y, ts);

		readInnerCircleFromSAE:for(ap_uint<8> i = 0; i < INNER_SIZE + 1; i = i + READ_NPC)
		{
	#pragma HLS DEPENDENCE variable=saeHW inter false
	#pragma HLS PIPELINE rewind
			if (i >= INNER_SIZE)
			{
				updateSAE(x, y, ts);
			}
			else
			{
	            ap_uint<8 * READ_NPC> xInnerTest, xOuterTest;
	            X_TYPE xInner[READ_NPC];
	            Y_TYPE yInner[READ_NPC], yInnerNewIdx[READ_NPC];
	        	rwSAEReadOffsetBitsLoop:
	            for (ap_uint<8> j = 0; j < 8 * READ_NPC; j++)
	            {
	#pragma HLS UNROLL
	            	ap_uint<8> tmpIndex;   // In  order to save the resource, we use bit operation to get the index of the inner/outer offset.
	            	tmpIndex.range(7, 2 + READ_NPC) = ap_uint<8>(i * 8).range(7, 2 + READ_NPC);
	            	tmpIndex.range(1 + READ_NPC, 0) = j.range(1 + READ_NPC, 0);
	            	xInnerTest[j] = innerTest[tmpIndex];
	            	xOuterTest[j] = outerTest[tmpIndex];
	            }

	            readNPCLoop:
	            for (ap_uint<8> k = 0; k < READ_NPC; k++)
	            {
	                xInner[k] = x + ap_int<4>(xInnerTest(8 * k + 3, 8  * k));  // Change back from unsigned to signed.
	                yInner[k] = y + ap_int<4>(xInnerTest(8 * k + 7, 8 * k + 4));          // Change back from unsigned to signed.
	                yInnerNewIdx[k] = yInner[k]%RESHAPE_FACTOR;

	                outputData[i + k] = readOneDataFromCol(saeHW[0][yInner[k]/RESHAPE_FACTOR][xInner[k]], yInnerNewIdx[k]);
	            }

	//			X_TYPE xOuter = x + xOuterTest(3, 0);
	//			Y_TYPE yOuter = y + xOuterTest(7, 4);
	//			Y_TYPE yOuterNewIdx = yOuter%RESHAPE_FACTOR;
	//
	//			outerCircle[i] = readOneDataFromCol(saeHW[0][yOuter/RESHAPE_FACTOR][xOuter], yOuterNewIdx);
			}
		}

		*size = INNER_SIZE;
	}
	else if(stage == 1)
	{
		readOuterCircleFromSAE:for(ap_uint<8> i = 0; i < OUTER_SIZE + 1; i = i + READ_NPC)
		{
	#pragma HLS DEPENDENCE variable=saeHW inter false
	#pragma HLS PIPELINE rewind
			if (i >= OUTER_SIZE)
			{
				updateSAE(x, y, ts);
			}
			else
			{
				ap_uint<8 * READ_NPC> xOuterTest;
				X_TYPE xOuter[READ_NPC];
				Y_TYPE yOuter[READ_NPC], yOuterNewIdx[READ_NPC];
				rwSAEReadOuterOffsetBitsLoop:
				for (ap_uint<8> j = 0; j < 8 * READ_NPC; j++)
				{
		#pragma HLS UNROLL
					ap_uint<8> tmpIndex; // In  order to save the resource, we use bit operation to get the index of the inner/outer offset.
					tmpIndex.range(7, 2 + READ_NPC) = ap_uint<8>(i * 8).range(7, 2 + READ_NPC);
					tmpIndex.range(1 + READ_NPC, 0) = j.range(1 + READ_NPC, 0);
					xOuterTest[j] = outerTest[tmpIndex];
				}

				readOuterNPCLoop:
				for (ap_uint<8> k = 0; k < READ_NPC; k++)
				{
					xOuter[k] = x + ap_int<4>(xOuterTest(8 * k + 3, 8  * k));       // Change back from unsigned to signed.
					yOuter[k] = y + ap_int<4>(xOuterTest(8 * k + 7, 8 * k + 4));    // Change back from unsigned to signed.
					yOuterNewIdx[k] = yOuter[k]%RESHAPE_FACTOR;

					outputData[i + k] = readOneDataFromCol(saeHW[0][yOuter[k]/RESHAPE_FACTOR][xOuter[k]], yOuterNewIdx[k]);
				}
			}
		}
		*size = OUTER_SIZE;
	}
	else
	{
		*size = 0;
	}
}


template<int READ_NPC>
void readOutterCircle(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> outerCircle[OUTER_SIZE])
{
	readOuterCircleFromSAE:for(ap_uint<8> i = 0; i < OUTER_SIZE; i = i + READ_NPC)
	{
#pragma HLS DEPENDENCE variable=saeHW inter false
#pragma HLS PIPELINE rewind

        ap_uint<8 * READ_NPC> xOuterTest;
        X_TYPE xOuter[READ_NPC];
        Y_TYPE yOuter[READ_NPC], yOuterNewIdx[READ_NPC];
    	rwSAEReadOuterOffsetBitsLoop:
        for (ap_uint<8> j = 0; j < 8 * READ_NPC; j++)
        {
#pragma HLS UNROLL
        	ap_uint<8> tmpIndex;
        	tmpIndex.range(7, 2 + READ_NPC) = ap_uint<8>(i * 8).range(7, 2 + READ_NPC);
        	tmpIndex.range(1 + READ_NPC, 0) = j.range(1 + READ_NPC, 0);
        	xOuterTest[j] = outerTest[tmpIndex];
        }

        readOuterNPCLoop:
        for (ap_uint<8> k = 0; k < READ_NPC; k++)
        {
            xOuter[k] = x + xOuterTest(8 * k + 3, 8  * k);
            yOuter[k] = y + xOuterTest(8 * k + 7, 8 * k + 4);
            yOuterNewIdx[k] = yOuter[k]%RESHAPE_FACTOR;

            outerCircle[i + k] = readOneDataFromCol(saeHW[0][yOuter[k]/RESHAPE_FACTOR][xOuter[k]], yOuterNewIdx[k]);
        }
	}
}

template<int DATA_SIZE, int NPC>
void insertionSortParallel(ap_uint<TS_TYPE_BIT_WIDTH> A[DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> B[DATA_SIZE])
{
    #pragma HLS array_partition variable=B complete
    L1:  for(int i = 0; i < DATA_SIZE; i++)
    {
        #pragma HLS pipeline II=1 rewind
    	ap_uint<TS_TYPE_BIT_WIDTH> item = A[i];
        L2:
        for (int j = DATA_SIZE - 1; j >= 0; j--)
        {
        	ap_uint<TS_TYPE_BIT_WIDTH> t;

        	ap_uint<1> cond1 = (j <= i) ? 1 : 0;
        	ap_uint<1> cond2 = (j > 0) ? 1 : 0;
        	ap_uint<1> cond3 = (cond2 && (B[j - 1] > item)) ? 1 : 0;

//        	B[j] = (cond1 == 1 && cond3 == 1) ? B[j-1] : B[j];
//        	B[j] = (cond1 == 1 && (cond3 == 0)) ? item : B[j];
//        	item = (cond1 == 1 && cond2 == 1  && cond3 == 0) ? B[j-1]: item;

        	if(cond1 == 1 && cond3 == 1) B[j] = B[j - 1];
        	else if((cond1 == 1 && (cond3 == 0))) B[j] = item;
        	if(cond1 == 1 && cond2 == 1  && cond3 == 0) item = B[j-1];

//        	if(j <= i)
//        	{
//				if(j > 0 && B[j - 1] > item)
//				{
//					B[j] = B[j - 1];
//				}
//				else
//				{
//					B[j] = item;
//					if (j > 0)
//					{
//						item = B[j - 1];
//					}
//				}
//        	}
        }
    }
}



static uint8_t sortedIndex[INNER_SIZE] = {0};
template<int CELL_ID>
void cellExt(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & in, hls::stream<ap_uint< TS_TYPE_BIT_WIDTH> > & out,
		ap_uint<TS_TYPE_BIT_WIDTH> initVal,
		hls::stream<uint8_t> & indexStream)
{
    const static uint8_t cellId = CELL_ID;    // The cell's id.
    static ap_uint<TS_TYPE_BIT_WIDTH> local = 0;
    static uint8_t dataIndex = 0;
    static uint8_t counter = 0;
    ap_uint<TS_TYPE_BIT_WIDTH> in_copy = in.read();

    if (counter == 0)
    {
    	local = initVal;
    }

    if(in_copy < initVal)
    {
    	dataIndex++;
    }

    counter++;

    if (counter >= INNER_SIZE)
    {
    	counter = 0;
    	indexStream.write(dataIndex);
    }

    out.write(in_copy);
}

void testCellExt(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & in, hls::stream<ap_uint< TS_TYPE_BIT_WIDTH> > & out,
		ap_uint<TS_TYPE_BIT_WIDTH> initVal[8],
		hls::stream<uint8_t>  indexStream[8])
{
#pragma HLS ARRAY_PARTITION variable=initVal complete dim=0
#pragma HLS DATAFLOW
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > tmpOut[7];

	cellExt<0>(in, tmpOut[0], initVal[0], indexStream[0]);
	cellExt<1>(tmpOut[0], tmpOut[1], initVal[1], indexStream[1]);
	cellExt<2>(tmpOut[1], tmpOut[2], initVal[2], indexStream[2]);
	cellExt<3>(tmpOut[2], tmpOut[3], initVal[3], indexStream[3]);
	cellExt<4>(tmpOut[3], tmpOut[4], initVal[4], indexStream[4]);
	cellExt<5>(tmpOut[4], tmpOut[5], initVal[5], indexStream[5]);
	cellExt<6>(tmpOut[5], tmpOut[6], initVal[6], indexStream[6]);
	cellExt<7>(tmpOut[6], out, initVal[7], indexStream[7]);
}

void insertionCellExtSort(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & in, hls::stream<ap_uint< TS_TYPE_BIT_WIDTH> > & out,
		ap_uint<TS_TYPE_BIT_WIDTH> initVal[8], ap_uint<TS_TYPE_BIT_WIDTH> outputIndex[8])
{
#pragma HLS ARRAY_PARTITION variable=initVal complete dim=0
#pragma HLS DATAFLOW
#pragma HLS ARRAY_PARTITION variable=outputIndex complete dim=0

    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > tmpOut[7];
    hls::stream<uint8_t>  indexStream[8];

	cellExt<0>(in, tmpOut[0], initVal[0], indexStream[0]);
	cellExt<1>(tmpOut[0], tmpOut[1], initVal[1], indexStream[1]);
	cellExt<2>(tmpOut[1], tmpOut[2], initVal[2], indexStream[2]);
	cellExt<3>(tmpOut[2], tmpOut[3], initVal[3], indexStream[3]);
	cellExt<4>(tmpOut[3], tmpOut[4], initVal[4], indexStream[4]);
	cellExt<5>(tmpOut[4], tmpOut[5], initVal[5], indexStream[5]);
	cellExt<6>(tmpOut[5], tmpOut[6], initVal[6], indexStream[6]);
	cellExt<7>(tmpOut[6], out, initVal[7], indexStream[7]);

	for(uint8_t i = 0; i < 8; i ++)
	{
#pragma HLS UNROLL
		outputIndex[i] = indexStream[i].read();
	}
}

template<int UNUSE_PARAMETER>
void cell(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & in, hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & out,
		hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > & localStream)
{
    static ap_uint< TS_TYPE_BIT_WIDTH > local = 0;
    static uint8_t counter = 0;
    ap_uint<TS_TYPE_BIT_WIDTH> in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
    counter++;
    if(counter == 8)
    {
    	counter = 0;
        localStream.write(local);
    }
}


void insertionCells(hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > &in, hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > &out,
		hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > localStream[8])
{
#pragma HLS DATAFLOW
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > tmpOut[7];

	cell<0>(in, tmpOut[0], localStream[0]);
	cell<1>(tmpOut[0], tmpOut[1], localStream[1]);
	cell<2>(tmpOut[1], tmpOut[2], localStream[2]);
	cell<3>(tmpOut[2], tmpOut[3], localStream[3]);
	cell<4>(tmpOut[3], tmpOut[4], localStream[4]);
	cell<5>(tmpOut[4], tmpOut[5], localStream[5]);
	cell<6>(tmpOut[5], tmpOut[6], localStream[6]);
	cell<7>(tmpOut[6], out, localStream[7]);
}

void insertionCellSort(ap_uint<TS_TYPE_BIT_WIDTH> inData[20], ap_uint<TS_TYPE_BIT_WIDTH> outputData[20])
{
#pragma HLS DATAFLOW
#pragma HLS ARRAY_PARTITION variable=outputData complete dim=0

    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > inStream, outStream;
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH> > localStream[8];

    for (uint8_t i = 0; i < 8; i++)
    {
#pragma HLS pipeline rewind
    	inStream.write(inData[i]);
    }

	insertionCells(inStream, outStream, localStream);

    for (uint8_t i = 0; i < 8; i++)
    {
#pragma HLS pipeline rewind
    	ap_uint<TS_TYPE_BIT_WIDTH> tmp = outStream.read();
    	outputData[i] = localStream[i].read();
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
	/* input */ ap_uint<8> num_symbols,
    /* output */ ap_uint<TS_TYPE_BIT_WIDTH> out[DATA_SIZE]) {
	ap_uint<TS_TYPE_BIT_WIDTH> previous_sorting[DATA_SIZE], sorting[DATA_SIZE];
    ap_uint<SYMBOL_BITS> digit_histogram[RADIX], digit_location[RADIX];
#pragma HLS ARRAY_PARTITION variable=digit_location complete dim=1
#pragma HLS ARRAY_PARTITION variable=digit_histogram complete dim=1

    Digit current_digit[DATA_SIZE];

 copy_in_to_sorting:
    for(int j = 0; j < num_symbols; j++) {
#pragma HLS PIPELINE II=1
        sorting[j] = in[j];
    }

 radix_sort_step1:
    for(int shift = 0; shift < TS_TYPE_BIT_WIDTH; shift += BITS_PER_LOOP) {
    init_histogram:
        for(int i = 0; i < RADIX; i++) {
#pragma HLS pipeline rewind
            digit_histogram[i] = 0;
        }
    }

radix_sort_step2:
   for(int shift = 0; shift < TS_TYPE_BIT_WIDTH; shift += BITS_PER_LOOP) {
    compute_histogram:
        for(int j = 0; j < num_symbols; j++) {
#pragma HLS pipeline rewind
            Digit digit = (sorting[j] >> shift) & (RADIX - 1); // Extract a digit
            current_digit[j] = digit;  // Store the current digit for each symbol
            digit_histogram[digit]++;
            previous_sorting[j] = sorting[j]; // Save the current sorted order of symbols
        }
     }

radix_sort_step3:
	  for(int shift = 0; shift < TS_TYPE_BIT_WIDTH; shift += BITS_PER_LOOP) {
        digit_location[0] = 0;
    find_digit_location:
        for(int i = 1; i < RADIX; i++)
#pragma HLS pipeline rewind
            digit_location[i] = digit_location[i-1] + digit_histogram[i-1];

	  }

radix_sort_step4:
  for(int shift = 0; shift < TS_TYPE_BIT_WIDTH; shift += BITS_PER_LOOP) {
    re_sort:
        for(int j = 0; j < num_symbols; j++) {
#pragma HLS pipeline rewind
            Digit digit = current_digit[j];
            sorting[digit_location[digit]] = previous_sorting[j]; // Move symbol to new sorted location
            out[digit_location[digit]] = previous_sorting[j]; // Also copy to output
            digit_location[digit]++; // Update digit_location
        }
    }
}


void mergeArraysWithSize(ap_uint<TS_TYPE_BIT_WIDTH> in[OUTER_SIZE], ap_uint<6> width,  ap_uint<5> size, ap_uint<TS_TYPE_BIT_WIDTH> out[OUTER_SIZE])
{
#pragma HLS FUNCTION_INSTANTIATE variable=width
#pragma HLS FUNCTION_INSTANTIATE variable=size

  assert(size <= 20);
  if(size > width)
  {
	  ap_uint<6> f1 = 0;
	  ap_uint<6> f2 = width;
	  ap_uint<6> i2 = width;
	  ap_uint<6> i3 = 2*width;
	  if(i2 >= size) i2 = size;
	  if(i3 >= size) i3 = size;
	 merge_arrays:
	  for (int i = 0; i < size; i++) {
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
	      if(i2 >= size) i2 = size;
	      if(i3 >= size) i3 = size;
	      f2 = i2;
	     }
	  }
  }
  else
  {
	  for (int i = 0; i < size; i++) {
	#pragma HLS LOOP_TRIPCOUNT min=0 max=20
	#pragma HLS PIPELINE II=1
		  out[i] = in[i];
	  }
  }
}

void mergeSortParallelWithSize(ap_uint<TS_TYPE_BIT_WIDTH> A[OUTER_SIZE], ap_uint<8> num_symbols,  ap_uint<TS_TYPE_BIT_WIDTH> B[OUTER_SIZE])
{
#pragma HLS DATAFLOW

    ap_uint<TS_TYPE_BIT_WIDTH> temp[5-1][20];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1
    int width = 1;

    mergeArraysWithSize(A, width, num_symbols, temp[0]);
    width *= 2;

    stage:
	for (int stage = 1; stage < 5-1; stage++) {
#pragma HLS UNROLL
		mergeArraysWithSize(temp[stage-1], width, num_symbols, temp[stage]);
        width *= 2;
    }

	mergeArraysWithSize(temp[5-2], width, num_symbols, B);
}


void testMergeArrays(ap_uint<TS_TYPE_BIT_WIDTH> in[OUTER_SIZE], ap_uint<6> width,  ap_uint<5> size, ap_uint<TS_TYPE_BIT_WIDTH> out[OUTER_SIZE])
{
#pragma HLS DATAFLOW

    ap_uint<TS_TYPE_BIT_WIDTH> temp[5-1][20];
#pragma HLS ARRAY_PARTITION variable=temp complete dim=1

	mergeArraysWithSize(in, 5, size, temp[0]);
	mergeArraysWithSize(temp[0], 10, size, out);
}

void testSortHW(ap_uint<TS_TYPE_BIT_WIDTH> inputA[TEST_SORT_DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outputB[TEST_SORT_DATA_SIZE])
{
//	 mergeSortParallel<TEST_SORT_DATA_SIZE, MERGE_STAGES> (inputA, outputB);
	insertionSortParallel<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
//	radixSort<TEST_SORT_DATA_SIZE, 1> (inputA, 16, outputB);
}

void testRwSAEHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<5> *size)
{
    rwSAE<2>(x, y, ts, stage, outputData, size);
}

void testIdxDataToIdxInnerBoolDataHW(ap_uint<5> newIdx[OUTER_SIZE], ap_uint<5> size, ap_uint<4> condFlg[OUTER_SIZE])
{
	idxDataToIdxInnerBoolData<6>(newIdx, size, condFlg);
}

void testFromTsDataToIdxDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> idxData[OUTER_SIZE])
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
#pragma HLS ARRAY_PARTITION variable=idxData cyclic factor=5 dim=0

    convertInterface<4>(inputRawData, INNER_SIZE, inStream);
    // The NPC value of sortedIdxStream should be equal to the value of idxData factor.
	sortedIdxStream<5>(inStream, INNER_SIZE, idxData);
//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)idxData[i]<< "\t";
//	}
//	std::cout << std::endl;
}

void testFromTsDataToIdxInnerBoolDataHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<4> idxBoolData[INNER_SIZE])
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> idxData[OUTER_SIZE];
//#pragma HLS ARRAY_PARTITION variable=idxData cyclic factor=5 dim=0

    convertInterface<2>(inputRawData, size, inStream);
    // The NPC value of sortedIdxStream should be equal to the value of idxData factor.
	sortedIdxStream<2>(inStream, size, idxData);
	idxDataToIdxInnerBoolData<4>(idxData, size, idxBoolData);
//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)idxData[i]<< "\t";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Idx Bool Data HW is: " << std::endl;
//	for (int i = 0; i < INNER_SIZE; i++)
//	{
//		std::cout << idxBoolData[i][3] << idxBoolData[i][2] << idxBoolData[i][1] << idxBoolData[i][0] << "\t";
//	}
//	std::cout << std::dec << std::endl;
}

void testFromTsDataToInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> idxData[OUTER_SIZE];
    ap_uint<4> idxBoolData[INNER_SIZE];
#pragma HLS RESOURCE variable=idxBoolData core=RAM_2P_LUTRAM

    convertInterface<2>(inputRawData, size, inStream);
    // The NPC value of sortedIdxStream should be equal to the value of idxData factor.
	sortedIdxStream<2>(inStream, size, idxData);
	idxDataToIdxInnerBoolData<4>(idxData, size, idxBoolData);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)idxData[i]<< "\t";
//	}
//	std::cout << std::endl;
//
//	std::cout << "Idx Bool Data HW is: " << std::endl;
//	for (int i = 0; i < INNER_SIZE; i++)
//	{
//		std::cout << idxBoolData[i][3] << idxBoolData[i][2] << idxBoolData[i][1] << idxBoolData[i][0] << "\t";
//	}
//	std::cout << std::dec << std::endl;

	idxInnerBoolDataToCorner<4>(idxBoolData, size, isCorner);

}


void testFromTsDataCheckInnerCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> idxData[OUTER_SIZE];

    convertInterface<2>(inputRawData, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)idxData[i]<< "\t";
//	}
//	std::cout << std::endl;

	checkInnerIdx<4>(idxData, INNER_SIZE, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void testFromTsDataCheckOuterCornerHW(ap_uint<TS_TYPE_BIT_WIDTH> inputRawData[OUTER_SIZE], ap_uint<5> size, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> idxData[OUTER_SIZE];

    convertInterface<2>(inputRawData, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)idxData[i]<< "\t";
//	}
//	std::cout << std::endl;

	checkOuterIdx<4>(idxData, OUTER_SIZE, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void fastCornerInnerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> size;
    ap_uint<5> idxData[OUTER_SIZE];

    rwSAE<2>(x, y, ts, stage, outer, &size);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)outer[i]<< "\t";
//	}
//	std::cout << std::endl;

//    sortedIdxData<2>(outer, size, idxData);
    convertInterface<4>(outer, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);
	checkInnerIdx<4>(idxData, size, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void fastCornerOuterHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> size;
    ap_uint<5> idxData[OUTER_SIZE];

    rwSAE<2>(x, y, ts, stage, outer, &size);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)outer[i]<< "\t";
//	}
//	std::cout << std::endl;

//    sortedIdxData<2>(outer, size, idxData);
    convertInterface<4>(outer, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);
	checkOuterIdx<4>(idxData, size, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void fastCornerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<1> *isCorner)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL
    ap_uint<5> size;
    ap_uint<5> idxData[OUTER_SIZE];

    rwSAE<2>(x, y, ts, stage, outer, &size);

//	std::cout << "Idx Data HW is: " << std::endl;
//	for (int i = 0; i < size; i++)
//	{
//		std::cout << (int)outer[i]<< "\t";
//	}
//	std::cout << std::endl;

//    sortedIdxData<2>(outer, size, idxData);
    convertInterface<4>(outer, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);
	checkInnerIdx<4>(idxData, size, isCorner);   // If resource is not enough, decrease this number to increase II a little.
}

void outputResult(ap_uint<1> isCorner,  hls::stream<apUint17_t> &packetEventDataStream, int32_t *eventSlice)
{
	apUint17_t tmp1 = packetEventDataStream.read();
//	apUint15_t miniSumRet = 0;
//	ap_int<9> tmp2 = miniSumRet.range(8, 0);
//	apUint6_t tmpOF = isCorner;

	ap_uint<32> output = tmp1;
	output[31] = isCorner;
//		std :: cout << "tmp1 is "  << std::hex << tmp1 << std :: endl;
//		std :: cout << "tmp2 is "  << std::hex << tmp2 << std :: endl;
//		std :: cout << "output is "  << std::hex << output << std :: endl;
//		std :: cout << "eventSlice is "  << std::hex << output.to_int() << std :: endl;
	*eventSlice++ = output.to_int();
}

void parseEventsHW(uint64_t * dataStream, int32_t eventsArraySize, int32_t *eventSlice)
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    hls::stream< ap_uint<TS_TYPE_BIT_WIDTH * OUTER_SIZE> > inStream("dataStream");
#pragma HLS STREAM variable=inStream depth=2 dim=1
#pragma HLS RESOURCE variable=inStream core=FIFO_SRL

	hls::stream<apUint17_t> pktEventDataStream("EventStream");
#pragma HLS STREAM variable=pktEventDataStream depth=2 dim=1
#pragma HLS RESOURCE variable=pktEventDataStream core=FIFO_SRL
    ap_uint<5> size;
    ap_uint<5> idxData[OUTER_SIZE];
	X_TYPE x;
	Y_TYPE y;
	ap_uint<TS_TYPE_BIT_WIDTH> ts;
	ap_uint<1> isCorner;

	getXandY(dataStream, &x, &y, &ts, pktEventDataStream);
    rwSAE<2>(x, y, ts, 0, outer, &size);
//    sortedIdxData<2>(outer, size, idxData);
    convertInterface<4>(outer, size, inStream);
	sortedIdxStream<4>(inStream, size, idxData);
	checkInnerIdx<4>(idxData, size, &isCorner);   // If resource is not enough, decrease this number to increase II a little.
	outputResult(isCorner, pktEventDataStream, eventSlice++);
}
