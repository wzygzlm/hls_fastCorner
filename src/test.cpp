#include "fastCorner.h"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;
#include "ap_fixed.h"
#include "time.h"
//#include "insertion_cell_sort.h"

const static int DEBUG=1;
const static int MAX_NUMBER=1000;
#define DTYPE ap_uint<32>
#define TEST_TIMES 20

// SAE (Surface of Active Event)
static uint32_t saeSW[1][DVS_HEIGHT][DVS_WIDTH];

const int innerCircleOffset[INNER_SIZE][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
      {3, 0}, {3, -1}, {2, -2}, {1, -3},
      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
const int outerCircleOffset[OUTER_SIZE][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
      {4, 1}, {4, 0}, {4, -1}, {3, -2},
      {2, -3}, {1, -4}, {0, -4}, {-1, -4},
      {-2, -3}, {-3, -2}, {-4, -1}, {-4, 0},
      {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}};

void insertion_sortSW(ap_uint<TS_TYPE_BIT_WIDTH> A[20], int size, ap_uint<TS_TYPE_BIT_WIDTH> sortOut[20]) {
    L1:
    for(int i = 1; i < size; i++) {
        DTYPE item = A[i];
        int j = i;
    DTYPE t = A[j-1];
        L2:
        while(j > 0 && A[j - 1] > item ) {
            #pragma HLS pipeline II=1
            A[j] = A[j - 1];
            j--;
        }
        A[j] = item;
    }
}

// Simple binary search algorithm
ap_uint<5> binarySearch(ap_uint<TS_TYPE_BIT_WIDTH> arr[], int l, int r, int x)
{
    if (r>=l)
    {
    	ap_uint<5> mid = l + (r - l)/2;
        if (arr[mid] == x)
            return mid;
        if (arr[mid] > x)
            return binarySearch(arr, l, mid-1, x);
        return binarySearch(arr, mid+1, r, x);
    }
    return -1;
}

// function takes an infinite size array and a key to be
//  searched and returns its position if found else -1.
// We don't know size of arr[] and we can assume size to be
// infinite in this function.
// NOTE THAT THIS FUNCTION ASSUMES arr[] TO BE OF INFINITE SIZE
// THEREFORE, THERE IS NO INDEX OUT OF BOUND CHECKING
ap_uint<5> findPos(ap_uint<TS_TYPE_BIT_WIDTH> arr[], ap_uint<TS_TYPE_BIT_WIDTH> key)
{
    int l = 0, h = 1;
    ap_uint<TS_TYPE_BIT_WIDTH> val = arr[0];

    // Find h to do binary search
    while (val < key)
    {
        l = h;        // store previous high
        h = 2*h;      // double high index
        val = arr[h]; // update new val
    }

    // at this point we have updated low and
    // high indices, Thus use binary search
    // between them
    return binarySearch(arr, l, h, key);
}

void sortedIndexSW(ap_uint<TS_TYPE_BIT_WIDTH> A[20], int size, ap_uint<5> sortOut[20])
{
	for(int i = 0; i < size; i++)
	{
		int temp = 0;
		for(int j = 0; j < size; j++)
		{
			if(A[j] < A[i]) temp += 1;
		}
		sortOut[i] = temp;
	}
}


void rwSAESW(X_TYPE x, Y_TYPE y, ap_uint<TS_TYPE_BIT_WIDTH> ts, ap_uint<2>  stage, ap_uint<TS_TYPE_BIT_WIDTH> outputData[OUTER_SIZE], ap_uint<5> *size)
{

	if(stage == 0)
	{
		saeSW[0][y][x] = ts;
		for(ap_uint<8> i = 0; i < INNER_SIZE; i = i + 1)
		{
			outputData[i] = saeSW[0][y + innerCircleOffset[i][1]][x + innerCircleOffset[i][0]];
		}
		*size = INNER_SIZE;
	}
	else if(stage == 1)
	{
		for(ap_uint<8> i = 0; i < OUTER_SIZE; i = i + 1)
		{
			outputData[i] = saeSW[0][y + outerCircleOffset[i][1]][x + outerCircleOffset[i][0]];
		}
		*size = OUTER_SIZE;
	}
	else
	{
		*size = 0;
	}
}


int main ()
 {
	int testTimes = TEST_TIMES;

    int total_err_cnt = 0;
	int retval=0;

	/******************* Test rwSAE module from random value**************************/
	srand((unsigned)time(NULL));
	int16_t eventCnt = 500;
	ap_uint<TS_TYPE_BIT_WIDTH> outputDataSW[OUTER_SIZE], outputDataHW[OUTER_SIZE];
	ap_uint<5> sizeSW, sizeHW;

	uint32_t x, y;
	uint32_t ts[eventCnt];

	for(int k = 0; k < TEST_TIMES; k++)
	{
		cout << "Test " << k << ":" << endl;

		int err_cnt = 0;

 		for (int i = 0; i < eventCnt; i++)
		{
 			ts[i]  = rand();
		}
 	    sort(ts, ts+eventCnt);

// 	    cout << "\nArray after sorting using "
// 	         "default sort is : \n";
// 	    for (int i = 0; i < eventCnt; ++i)
// 	        cout << ts[i] << " ";

 		for (int i = 0; i < eventCnt; i++)
		{
			x = rand()%50 + 20;
			y = rand()%50 + 20;
//			idx = rand()%3;
	//		x = 255;
	//		y = 240;
//			cout << "x : " << x << endl;
//			cout << "y : " << y << endl;
//			cout << "idx : " << idx << endl;

//			data[i] = (uint64_t)(x << 17) + (uint64_t)(y << 2) + (1 << 1);
//			cout << "data[" << i << "] is: "<< hex << data[i]  << endl;
	 		rwSAESW(x, y, ts[i], 1, outputDataSW, &sizeSW);
	 		testRwSAEHW(x, y, ts[i], 1, outputDataHW, &sizeHW);

	 		for (int  j = 0; j < OUTER_SIZE; j++)
	 		{
	 			if (outputDataSW[j] != outputDataHW[j])
	 			{
 	 				err_cnt++;
	 			}
	 		}
		}


 		if(err_cnt == 0)
		{
			cout << "Test " << k << " passed." << endl;
		}
		total_err_cnt += err_cnt;
		cout << endl;
	}

	/******************* Test SortedIdxData module from random value**************************/
//	ap_uint<TS_TYPE_BIT_WIDTH> input[OUTER_SIZE];
//	ap_uint<5> outputSortedIdxHW[OUTER_SIZE], outputSortedIdxSW[OUTER_SIZE];
//	for(int k = 0; k < TEST_TIMES; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//		//generate random data to sort
//		if(DEBUG) std::cout << "Random Input Data\n";
//		int size = rand()%20;
//		size = 20;
//		for(int i = 0; i < size; i++) {
//			input[i] = rand() % MAX_NUMBER + 1;
//			if(DEBUG) std::cout << input[i] << "\t";
//		}
//
//		testSortedIdxData(input, outputSortedIdxHW);
//		sortedIndex(input, size, outputSortedIdxSW);
//
//		//compare the results of insertion_sort to insertion_cell_sort; fail if they differ
//		if(DEBUG) std::cout << "\nSorted Output\n";
//		for(int i = 0; i < size; i++) {
//			if(DEBUG) std::cout << outputSortedIdxHW[i] << "\t";
//		}
//		for(int i = 0; i < size; i++) {
//			if(outputSortedIdxSW[i] != outputSortedIdxHW[i]) {
//				std::cout << "\n";
//				err_cnt = 1;
//				std::cout << "golden= " << outputSortedIdxSW[i] << " hw=" << outputSortedIdxHW[i] << "\n";
//			}
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}
	/******************* Test sort module from random value**************************/
//    srand((unsigned)time(NULL)); //change me if you want different numbers
//	srand(20);
//	int32_t eventCnt = 500;
//	uint64_t data[eventCnt];
//	int32_t eventSlice[eventCnt], eventSliceSW[eventCnt];
//
//    int fail = 0;
//    DTYPE input[20];
//    DTYPE insertion_output[20] = {0}, insertion_cell_output[20] = {0};
//    DTYPE merge_sort_input[20]= {0}, merge_sort_output[20] = {0};
//    hls::stream<DTYPE> in, out;
//
//	for(int k = 0; k < TEST_TIMES; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//	    int err_cnt = 0;
//
//	    //generate random data to sort
//	    if(DEBUG) std::cout << "Random Input Data\n";
// 		int size = rand()%20;
// 		size = 8;
//	    for(int i = 0; i < size; i++) {
//	        input[i] = rand() % MAX_NUMBER + 1;
//	        merge_sort_input[i] = input[i];
//	        if(DEBUG) std::cout << input[i] << "\t";
//	    }
//
//	    //process the data through the insertion_cell_sort function
////	    for(int i = 0; i < size*2; i++) {
////	        if(i < size) {
////	            //feed in the SIZE elements to be sorted
////	            in.write(input[i]);
////	            insertion_cell_sort(in, out);
////	            insertion_cell_output[i] = out.read();
////	        } else {
////	            //then send in dummy data to flush pipeline
////	            in.write(MAX_NUMBER);
////	            insertion_cell_sort(in, out);
////	            insertion_cell_output[i-size] = out.read();
////	        }
////	    }
//
////	    insertionCellSort(input, insertion_cell_output);
//	    testSortHW(input, insertion_cell_output);
//	    //sort the data using the insertion_sort function
//	    insertion_sortSW(input, size, insertion_output);
//	    mergeSortParallelWithSize(merge_sort_input, size, merge_sort_output);
//
//	    //compare the results of insertion_sort to insertion_cell_sort; fail if they differ
//	    if(DEBUG) std::cout << "\nSorted Output\n";
//	    for(int i = 0; i < size; i++) {
//	        if(DEBUG) std::cout << insertion_cell_output[i] << "\t";
//	    }
//	    for(int i = 0; i < size; i++) {
//	        if(input[i] != insertion_cell_output[i]) {
//	        	err_cnt = 1;
//	            std::cout << "golden= " << input[i] << " hw=" << insertion_cell_output[i] << "\n";
//	        }
//	    }
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}


	if (total_err_cnt == 0)
	{
			cout<<"*** TEST PASSED ***" << endl;
			retval = 0;
	} else
	{
			cout<<"!!! TEST FAILED - " << total_err_cnt << " mismatches detected !!!";
			cout<< endl;
			retval = -1;
	}

	// Return 0 if the test passes
	return retval;
}
