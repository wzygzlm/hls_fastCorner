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

const static int DEBUG=1;
const static int MAX_NUMBER=1000;
#define DTYPE ap_uint<TS_TYPE_BIT_WIDTH>
#define TEST_TIMES 20

void insertion_sortSW(DTYPE A[20], int size, DTYPE sortOut[20]) {
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

int main ()
{
	int testTimes = TEST_TIMES;

    int total_err_cnt = 0;
	int retval=0;

	/******************* Test parseEvents module from random value**************************/
    srand((unsigned)time(NULL)); //change me if you want different numbers
	int32_t eventCnt = 500;
	uint64_t data[eventCnt];
	int32_t eventSlice[eventCnt], eventSliceSW[eventCnt];

    int fail = 0;
    DTYPE input[20];
    DTYPE cell_output[20] = {0};
    DTYPE merge_sort_input[20]= {0}, merge_sort_output[20] = {0};
    hls::stream<DTYPE> in, out;

	for(int k = 0; k < TEST_TIMES; k++)
	{
		cout << "Test " << k << ":" << endl;

	    int err_cnt = 0;

	    //generate random data to sort
	    if(DEBUG) std::cout << "Random Input Data\n";
 		int size = rand()%20;
// 		size = 3;
	    for(int i = 0; i < size; i++) {
	        input[i] = rand() % MAX_NUMBER + 1;
	        merge_sort_input[i] = input[i];
	        if(DEBUG) std::cout << input[i] << "\t";
	    }

	    //process the data through the insertion_cell_sort function
	//    for(int i = 0; i < SIZE*2; i++) {
	//        if(i < SIZE) {
	//            //feed in the SIZE elements to be sorted
	//            in.write(input[i]);
	//            insertion_cell_sort(in, out);
	//            cell_output[i] = out.read();
	//        } else {
	//            //then send in dummy data to flush pipeline
	//            in.write(MAX_NUMBER);
	//            insertion_cell_sort(in, out);
	//            cell_output[i-SIZE] = out.read();
	//        }
	//    }

	    //sort the data using the insertion_sort function
	    insertion_sortSW(input, size, cell_output);
	    mergeSortParallelWithSize(merge_sort_input, size, merge_sort_output);

	    //compare the results of insertion_sort to insertion_cell_sort; fail if they differ
	    if(DEBUG) std::cout << "\nSorted Output\n";
	    for(int i = 0; i < size; i++) {
	        if(DEBUG) std::cout << merge_sort_output[i] << "\t";
	    }
	    for(int i = 0; i < size; i++) {
	        if(input[i] != merge_sort_output[i]) {
	        	err_cnt = 1;
	            std::cout << "golden= " << input[i] << " hw=" << merge_sort_output[i] << "\n";
	        }
	    }

		if(err_cnt == 0)
		{
			cout << "Test " << k << " passed." << endl;
		}
		total_err_cnt += err_cnt;
		cout << endl;
	}


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
