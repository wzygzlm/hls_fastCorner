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

const ap_int<128> innerTest =  ap_int<128>("03132231303f2e1d0dfdeedfd0d1e2f3", 16);
const ap_int<160> outerTest = ap_int<160>("0414233241404f3e2d1c0cfceddecfc0c1d2e3f4", 16);


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


template<int READ_NPC>   //  Due to the memory has 2 ports at most for arbitrary reading, here this number could be only 1 or 2.
void rwSAE(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<8> *size)
{
	if(stage == 0)
	{
		updateSAE(x, y, ts);

		readInnerCircleFromSAE:for(ap_uint<8> i = 0; i < INNER_SIZE; i = i + READ_NPC)
		{
	#pragma HLS DEPENDENCE variable=saeHW inter false
	#pragma HLS PIPELINE rewind
//			if (i >= INNER_SIZE)
//			{
//				updateSAE(x, y, ts);
//			}
//			else
//			{
	            ap_uint<8 * READ_NPC> xInnerTest, xOuterTest;
	            X_TYPE xInner[READ_NPC];
	            Y_TYPE yInner[READ_NPC], yInnerNewIdx[READ_NPC];
	        	rwSAEReadOffsetBitsLoop:
	            for (ap_uint<8> j = 0; j < 8 * READ_NPC; j++)
	            {
	#pragma HLS UNROLL
	            	ap_uint<8> tmpIndex;
	            	tmpIndex.range(7, 2 + READ_NPC) = ap_uint<8>(i * 8).range(7, 2 + READ_NPC);
	            	tmpIndex.range(1 + READ_NPC, 0) = j.range(1 + READ_NPC, 0);
	            	xInnerTest[j] = innerTest[tmpIndex];
	            	xOuterTest[j] = outerTest[tmpIndex];
	            }

	            readNPCLoop:
	            for (ap_uint<8> k = 0; k < READ_NPC; k++)
	            {
	                xInner[k] = x + xInnerTest(8 * k + 3, 8  * k);
	                yInner[k] = y + xInnerTest(8 * k + 7, 8 * k + 4);
	                yInnerNewIdx[k] = yInner[k]%RESHAPE_FACTOR;

	                outputData[i + k] = readOneDataFromCol(saeHW[0][yInner[k]/RESHAPE_FACTOR][xInner[k]], yInnerNewIdx[k]);
	            }

	//			X_TYPE xOuter = x + xOuterTest(3, 0);
	//			Y_TYPE yOuter = y + xOuterTest(7, 4);
	//			Y_TYPE yOuterNewIdx = yOuter%RESHAPE_FACTOR;
	//
	//			outerCircle[i] = readOneDataFromCol(saeHW[0][yOuter/RESHAPE_FACTOR][xOuter], yOuterNewIdx);
//			}
		}

		*size = INNER_SIZE;
	}
	else if(stage == 1)
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

	            outputData[i + k] = readOneDataFromCol(saeHW[0][yOuter[k]/RESHAPE_FACTOR][xOuter[k]], yOuterNewIdx[k]);
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


void mergeArraysWithSize(ap_uint<TS_TYPE_BIT_WIDTH> in[OUTER_SIZE], ap_uint<8> width,  ap_uint<8> size, ap_uint<TS_TYPE_BIT_WIDTH> out[OUTER_SIZE])
{
  if(size != 0)
  {
		int f1 = 0;
	  int f2 = width;
	  int i2 = width;
	  int i3 = 2*width;
	  if(i2 >= size) i2 = size;
	  if(i3 >= size) i3 = size;
	 merge_arrays:
	  for (int i = 0; i < size; i++) {
	#pragma HLS LOOP_TRIPCOUNT min=0 max=20
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




void testSortHW(ap_uint<TS_TYPE_BIT_WIDTH> inputA[TEST_SORT_DATA_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outputB[TEST_SORT_DATA_SIZE])
{
//	 mergeSortParallel<TEST_SORT_DATA_SIZE, MERGE_STAGES> (inputA, outputB);
//	insertionSortParallel<TEST_SORT_DATA_SIZE, 1> (inputA, outputB);
	radixSort<TEST_SORT_DATA_SIZE, 1> (inputA, 16, outputB);
}

void fastCornerInnerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage,
		ap_uint<TS_TYPE_BIT_WIDTH> innerSort[INNER_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> outerSort[OUTER_SIZE])
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> inner[INNER_SIZE];
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
    ap_uint<8> size;
    int8_t index;
    rwSAE<2>(x, y, ts, stage, outer, &size);
    mergeSortParallelWithSize(outer, size, outerSort);
//    return min< ap_uint<TS_TYPE_BIT_WIDTH>, INNER_SIZE >(inner, &index);

//    insertionSortParallel<OUTER_SIZE, 1>(outer, outerSort);
//	mergeSortParallel<INNER_SIZE, 4>(inner, innerSort);
}

ap_uint<TS_TYPE_BIT_WIDTH> fastCornerOuterHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> Outer_B[OUTER_SIZE])
{
#pragma HLS DATAFLOW
    ap_uint<TS_TYPE_BIT_WIDTH> outer[OUTER_SIZE];
	readOutterCircle<2>(x, y, outer);
    int8_t index;
    return min< ap_uint<TS_TYPE_BIT_WIDTH>, OUTER_SIZE >(outer, &index);

//	mergeSortParallel<OUTER_SIZE, MERGE_STAGES>(outer, Outer_B);
//		insertionSortParallel<OUTER_SIZE, 1>(outer, Outer_B);
}

void fastCornerHW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<1>  skip,
		ap_uint<TS_TYPE_BIT_WIDTH> Inner_B[INNER_SIZE], ap_uint<TS_TYPE_BIT_WIDTH> Outer_B[OUTER_SIZE],
		ap_uint<TS_TYPE_BIT_WIDTH> *minimum)
{
//	if(skip == 0)
//	{
//		*minimum = fastCornerInnerHW(x, y, ts, Inner_B);
//	}
//	else
//	{
//		*minimum = fastCornerOuterHW(x, y, Outer_B);
//	}
}
