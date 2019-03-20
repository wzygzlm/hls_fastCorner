#include "insertion_cell_sort.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"
#include "assert.h"

void insertion_sort(DTYPE A[SIZE])
{
    L1:
    for (int i = 1; i < SIZE; i++)
    {
    	DTYPE item = A[i];
        int j = i;
        DTYPE t = A[j-1];
        L2:
		while (j > 0 && t > item)
		{
			#pragma HLS pipeline II=1
			A[j] = t;
			t = A[j - 2];
			j--;
		}
        A[j] = item;
    }
}
template<int NPC>
void insertion_sort_parallel(DTYPE A[SIZE], DTYPE B[SIZE]) {
    #pragma HLS array_partition variable=B complete
    L1:  for(int i = 0; i < SIZE; i++)
    {
        #pragma HLS pipeline II=1
    	DTYPE item = A[i];
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

void cell0(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell1(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell2(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell3(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell4(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell5(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell6(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}
void cell7(hls::stream<DTYPE> & in, hls::stream<DTYPE> & out)
{
    static DTYPE local = 0;
    DTYPE in_copy = in.read();
    if(in_copy > local) {
        out.write(local);
        local = in_copy;
    }
    else
    {
        out.write(in_copy);
    }
}

// template<int NPC>
void insertion_cell_sort(hls::stream<DTYPE> &in, hls::stream<DTYPE> &out)
{
    #pragma HLS DATAFLOW
    hls::stream<DTYPE> out0("out0 stream");
    hls::stream<DTYPE> out1("out1 stream");
    hls::stream<DTYPE> out2("out2 stream");
    hls::stream<DTYPE> out3("out3 stream");
    hls::stream<DTYPE> out4("out4 stream");
    hls::stream<DTYPE> out5("out5 stream");
    hls::stream<DTYPE> out6("out6 stream");

    // Function calls;
    cell0(in, out0);
    cell1(out0, out1);
    cell2(out1, out2);
    cell3(out2, out3);
    cell4(out3, out4);
    cell5(out4, out5);
    cell6(out5, out6);
    cell7(out6, out);
}

void merge_sort(DTYPE A[SIZE]) {
    DTYPE temp[SIZE];
 stage:
    for (int width = 1; width < SIZE; width = 2 * width) {
        int f1 = 0;
        int f2 = width;
        int i2 = width;
        int i3 = 2*width;
        if(i2 >= SIZE) i2 = SIZE;
        if(i3 >= SIZE) i3 = SIZE;
    merge_arrays:
        for (int i = 0; i < SIZE; i++) {
#pragma HLS pipeline II=1
            DTYPE t1 = A[f1];
            DTYPE t2 = A[f2];
            if((f1 < i2 && t1 <= t2) || f2 == i3) {
                temp[i] = t1;
                f1++;
            } else {
                assert(f2 < i3);
                temp[i] = t2;
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

    copy:
        for(int i = 0; i < SIZE; i++) {
#pragma HLS pipeline II=1
            A[i] = temp[i];
        }
    }
}

void merge_arrays(DTYPE in[SIZE], int width, DTYPE out[SIZE]) {
  int f1 = 0;
  int f2 = width;
  int i2 = width;
  int i3 = 2*width;
  if(i2 >= SIZE) i2 = SIZE;
  if(i3 >= SIZE) i3 = SIZE;
 merge_arrays:
  for (int i = 0; i < SIZE; i++) {
#pragma HLS pipeline II=1
      DTYPE t1 = in[f1];
      DTYPE t2 = in[f2];
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

void merge_sort_parallel(DTYPE A[SIZE], DTYPE B[SIZE]) {
#pragma HLS dataflow

    DTYPE temp[STAGES-1][SIZE];
#pragma HLS array_partition variable=temp complete dim=1
    int width = 1;

    merge_arrays(A, width, temp[0]);
    width *= 2;

    for (int stage = 1; stage < STAGES-1; stage++) {
#pragma HLS unroll
        merge_arrays(temp[stage-1], width, temp[stage]);
        width *= 2;
    }

    merge_arrays(temp[STAGES-2], width, B);
}


const unsigned int RADIX = 16;
const unsigned int BITS_PER_LOOP = 4; // should be log2(RADIX)
typedef ap_uint<BITS_PER_LOOP> Digit;

void radix_sort(
    /* input */ DTYPE in[SIZE],
    /* output */ DTYPE out[SIZE]) {
	DTYPE previous_sorting[SIZE], sorting[SIZE];
    ap_uint<SYMBOL_BITS> digit_histogram[RADIX], digit_location[RADIX];
#pragma HLS ARRAY_PARTITION variable=digit_location complete dim=1
#pragma HLS ARRAY_PARTITION variable=digit_histogram complete dim=1
    Digit current_digit[SIZE];

 copy_in_to_sorting:
    for(int j = 0; j < SIZE; j++) {
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
        for(int j = 0; j < SIZE; j++) {
#pragma HLS PIPELINE II=1
            Digit digit = (sorting[j] >> shift) & (RADIX - 1); // Extrract a digit
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
        for(int j = 0; j < SIZE; j++) {
#pragma HLS PIPELINE II=1
            Digit digit = current_digit[j];
            sorting[digit_location[digit]] = previous_sorting[j]; // Move symbol to new sorted location
            out[digit_location[digit]] = previous_sorting[j]; // Also copy to output
            digit_location[digit]++; // Update digit_location
        }
    }
}


void testSort(ap_uint<20> inputA[16], ap_uint<20> outputB[16])
{
	insertion_sort_parallel<1> (inputA, outputB);
}
