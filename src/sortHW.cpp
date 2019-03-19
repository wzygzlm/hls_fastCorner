#include "insertion_cell_sort.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"

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


void testSort(ap_uint<20> inputA[16], ap_uint<20> outputB[16])
{
	insertion_sort_parallel<1> (inputA, outputB);
}
