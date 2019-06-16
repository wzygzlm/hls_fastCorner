#include "fastCorner.h"
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

const static int MAX_NUMBER=1000;
#define DTYPE ap_uint<32>
#define TEST_TIMES 200

static const int sensor_width_ = 240;
static const int sensor_height_ = 180;

// SAE (Surface of Active Event)
static int sae_[2][DVS_HEIGHT][DVS_WIDTH];

const int circle3_[INNER_SIZE][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
      {3, 0}, {3, -1}, {2, -2}, {1, -3},
      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
const int circle4_[OUTER_SIZE][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
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
		sae_[0][y][x] = ts;
		for(ap_uint<8> i = 0; i < INNER_SIZE; i = i + 1)
		{
			outputData[i] = sae_[0][y + circle3_[i][1]][x + circle3_[i][0]];
		}
		*size = INNER_SIZE;
	}
	else if(stage == 1)
	{
		for(ap_uint<8> i = 0; i < OUTER_SIZE; i = i + 1)
		{
			outputData[i] = sae_[0][y + circle4_[i][1]][x + circle4_[i][0]];
		}
		*size = OUTER_SIZE;
	}
	else
	{
		*size = 0;
	}
}


// Compares two intervals according to staring times.
bool compareInterval(pair<uint32_t, int> i1, pair<uint32_t, int> i2)
{
    return (i1.first < i2.first);
}


void sortArr(uint32_t arr[], int n, uint8_t outputIdx[])
{

    // Vector to store element
    // with respective present index
    vector<pair<uint32_t, int> > vp;

    // Inserting element in pair vector
    // to keep track of previous indexes
    for (int i = 0; i < n; ++i) {
        vp.push_back(make_pair(arr[i], i));
    }

    // Sorting pair vector
    sort(vp.begin(), vp.end(), compareInterval);

    // Displaying sorted element
    // with previous indexes
    // corresponding to each element
//    cout << "Element\t"
//         << "index" << endl;
    for (int i = 0; i < vp.size(); i++) {
    	outputIdx[vp[i].second] = i;
        cout << vp[i].first << "\t"
             << vp[i].second << "\t" << i << endl;
    }
}

void testConvertandSortedIdxSW(uint32_t rawData[OUTER_SIZE], uint8_t size, uint8_t outputIdxData[OUTER_SIZE])
{
	assert(size <= OUTER_SIZE);
	cout << "Raw Data is: " << hex << endl;
	for (int i = 0; i < size; i++)
	{
		cout << rawData[i] << "\t";
	}
	cout << dec << endl;

	sortArr(rawData, size, outputIdxData);

	cout << "Idx Data SW is: " << endl;
	for (int i = 0; i < size; i++)
	{
		cout << (int)outputIdxData[i]<< "\t";
	}
	cout << endl;
}

void testConvertandSortedInnerIdxSW(uint32_t rawData[OUTER_SIZE], ap_uint<4> condIdxData[INNER_SIZE])
{
	uint8_t outputIdxData[OUTER_SIZE];
	testConvertandSortedIdxSW(rawData, INNER_SIZE, outputIdxData);

	for (int i = 0; i < INNER_SIZE; i++)
	{
		for (int streak_size = 3; streak_size<=6; streak_size++)
		{
			condIdxData[i][streak_size - 3] =  (outputIdxData[i] >= (INNER_SIZE - streak_size));
		}
	}

	cout << "Idx Bool Data SW is: " << endl;
	for (int i = 0; i < INNER_SIZE; i++)
	{
		cout << condIdxData[i][3] << condIdxData[i][2] << condIdxData[i][1] << condIdxData[i][0] << "\t";
	}
	cout << dec << endl;
}

void testIdxDataToIdxInnerBoolDataSW(int idxData[OUTER_SIZE], ap_uint<5> size, ap_uint<5> condFlg[OUTER_SIZE])
{
	if(size == INNER_SIZE || size == 18)
	{
		for(int i = 0; i < size; i++)
		{
			condFlg[i][0] = (idxData[i] >= INNER_SIZE - 3 + OUTER_SIZE - INNER_SIZE);
			condFlg[i][1] = (idxData[i]  >= INNER_SIZE - 4 + OUTER_SIZE - INNER_SIZE);
			condFlg[i][2] = (idxData[i]  >= INNER_SIZE - 5 + OUTER_SIZE - INNER_SIZE);
			condFlg[i][3] = (idxData[i]  >= INNER_SIZE - 6 + OUTER_SIZE - INNER_SIZE);
		}
	}
	else if(size == OUTER_SIZE)
	{
		for(int i = 0; i < size; i++)
		{
			condFlg[i][0] = (idxData[i] >= OUTER_SIZE - 4);
			condFlg[i][1] = (idxData[i] >= OUTER_SIZE - 5);
			condFlg[i][2] = (idxData[i] >= OUTER_SIZE - 6);
			condFlg[i][3] = (idxData[i] >= OUTER_SIZE - 7);
			condFlg[i][4] = (idxData[i] >= OUTER_SIZE - 8);
		}
	}

}

void testFromTsDataCheckInnerCornerSW(uint32_t rawData[OUTER_SIZE], uint8_t size, ap_uint<1> *isCorner)
{
	uint8_t idxData[INNER_SIZE];

	testConvertandSortedIdxSW(rawData, INNER_SIZE, idxData);

	*isCorner = 0;

	for (int streak_size = 3; streak_size<=6; streak_size++)
	{
		for (uint8_t i = 0; i < INNER_SIZE; i++)
		{
			ap_uint<1> tempCond = 1;
			uint8_t j =  0;
			for (j =  0; j < streak_size; j++)
			{
				uint8_t tmpData = ((i + j) >= INNER_SIZE) ? idxData[i + j - INNER_SIZE] : idxData[i + j];
				tempCond &= (tmpData >= (INNER_SIZE - streak_size));
			}
			if (tempCond == 1)
			{
				*isCorner = 1;
				cout << "Position is :" << (int)i << " and streak size is: " << (int)j << endl;
				return;
			}
		}
	}

}

void testFromTsDataCheckOuterCornerSW(uint32_t rawData[OUTER_SIZE], uint8_t size, ap_uint<1> *isCorner)
{
	uint8_t idxData[OUTER_SIZE];

	testConvertandSortedIdxSW(rawData, OUTER_SIZE, idxData);

	*isCorner = 0;

	for (int streak_size = 4; streak_size<=8; streak_size++)
	{
		for (uint8_t i = 0; i < OUTER_SIZE; i++)
		{
			ap_uint<1> tempCond = 1;
			uint8_t j =  0;
			for (j =  0; j < streak_size; j++)
			{
				uint8_t tmpData = ((i + j) >= OUTER_SIZE) ? idxData[i + j - OUTER_SIZE] : idxData[i + j];
				tempCond &= (tmpData >= (OUTER_SIZE - streak_size));
			}
			if (tempCond == 1)
			{
				*isCorner = 1;
				cout << "Position is :" << (int)i << " and streak size is: " << (int)j << endl;
				return;
			}
		}
	}

}

void FastDetectorisInnerFeature(int pix_x, int pix_y, int timesmp, bool polarity, bool *found_streak)
{

  // update SAE
  //const int pol = polarity ? 1 : 0; //conver plo to 1 or 0
  const int pol = 0;
  sae_[pol][pix_x][pix_y] = timesmp;//

  const int max_scale = 1;

  // only check if not too close to border
  const int cs = max_scale*4;
  if (pix_x < cs || pix_x >= sensor_width_-cs ||
      pix_y < cs || pix_y >= sensor_height_-cs)
  {
    *found_streak = false;
  }

  *found_streak = false;

  std::cout << "Idx Inner Data SW is: " << std::endl;
  for (int i=0; i<16; i++)
  {
		cout << sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]] << "\t";
  }
	std::cout << std::endl;


  isFeatureOutterLoop:for (int i=0; i<16; i++)
  {
    FastDetectorisFeature_label2:for (int streak_size = 3; streak_size<=6; streak_size++)
    {
      // check that streak event is larger than neighbor
      if ((sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]]) <  (sae_[pol][pix_x+circle3_[(i-1+16)%16][0]][pix_y+circle3_[(i-1+16)%16][1]]))
        continue;

      // check that streak event is larger than neighbor
      if (sae_[pol][pix_x+circle3_[(i+streak_size-1)%16][0]][pix_y+circle3_[(i+streak_size-1)%16][1]] < sae_[pol][pix_x+circle3_[(i+streak_size)%16][0]][pix_y+circle3_[(i+streak_size)%16][1]])
        continue;

      // find the smallest timestamp in corner min_t
      double min_t = sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]];
      FastDetectorisFeature_label1:for (int j=1; j<streak_size; j++)
      {
        const double tj = sae_[pol][pix_x+circle3_[(i+j)%16][0]][pix_y+circle3_[(i+j)%16][1]];
        if (tj < min_t)
          min_t = tj;
      }

      //check if corner timestamp is higher than corner
      bool did_break = false;
      FastDetectorisFeature_label0:for (int j=streak_size; j<16; j++)
      {
        const double tj = sae_[pol][pix_x+circle3_[(i+j)%16][0]][pix_y+circle3_[(i+j)%16][1]];

        if (tj >= min_t)
        {
          did_break = true;
          break;
        }
      }

      if (!did_break)
      {
        *found_streak = true;
		cout << "Inner Corner position is : " << i << " and streak size is: " << streak_size << endl;
        break;
      }
    }

    if (*found_streak)
    {
      break;
    }
  }
}

void FastDetectorisOuterFeature(int pix_x, int pix_y, int timesmp, bool polarity, bool *found_streak)
{
	// update SAE
	//const int pol = polarity ? 1 : 0; //conver plo to 1 or 0
	const int pol = 0;
	sae_[pol][pix_x][pix_y] = timesmp;//

	const int max_scale = 1;

	// only check if not too close to border
	const int cs = max_scale*4;
	if (pix_x < cs || pix_x >= sensor_width_-cs ||
	  pix_y < cs || pix_y >= sensor_height_-cs)
	{
	*found_streak = false;
	}

	*found_streak = false;

	std::cout << "Idx Outer Data SW is: " << std::endl;
	for (int i=0; i<20; i++)
	{
	cout << sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]] << "\t";
	}
	std::cout << std::endl;

	FastDetectorisFeature_label6:for (int i=0; i<20; i++)
	{

		FastDetectorisFeature_label5:for (int streak_size = 4; streak_size<=8; streak_size++)
		{
			// check that first event is larger than neighbor
			if (sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]] <  sae_[pol][pix_x+circle4_[(i-1+20)%20][0]][pix_y+circle4_[(i-1+20)%20][1]])
			  continue;

			// check that streak event is larger than neighbor
			if (sae_[pol][pix_x+circle4_[(i+streak_size-1)%20][0]][pix_y+circle4_[(i+streak_size-1)%20][1]] < sae_[pol][pix_x+circle4_[(i+streak_size)%20][0]][pix_y+circle4_[(i+streak_size)%20][1]])
			  continue;

			double min_t = sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]];
			FastDetectorisFeature_label4:for (int j=1; j<streak_size; j++)
			{
				  const double tj = sae_[pol][pix_x+circle4_[(i+j)%20][0]][pix_y+circle4_[(i+j)%20][1]];
				  if (tj < min_t)
					min_t = tj;
			}

			bool did_break = false;
			FastDetectorisFeature_label3:for (int j=streak_size; j<20; j++)
			{
			  const double tj = sae_[pol][pix_x+circle4_[(i+j)%20][0]][pix_y+circle4_[(i+j)%20][1]];
			  if (tj >= min_t)
			  {
				did_break = true;
				break;
			  }
			}

			if (!did_break)
			{
			  *found_streak = true;
			  cout << "Outer Corner position is : " << i << " and streak size is: " << streak_size << endl;
			  break;
			}
		}
		if (*found_streak)
		{
			break;
		}
	}
}


void FastDetectorisFeature(int pix_x, int pix_y, int timesmp, bool polarity, bool *found_streak)
{
	if(polarity==0)
	{
		*found_streak = false;
		return;
	}

	if(pix_y == 159)
	{
		int tmp = 0;
	}

	const int max_scale = 1;
	// only check if not too close to border
	const int cs = max_scale*20;
	if (pix_x < cs || pix_x >= sensor_width_-cs - 4 ||
			pix_y < cs || pix_y >= sensor_height_-cs -4)
	{
		*found_streak = false;
		return;
	}

	const int pol = 0;
	// update SAE
	sae_[pol][pix_x][pix_y] = timesmp;

    *found_streak = false;

#if DEBUG
  std::cout << "Idx Inner Data SW is: " << std::endl;
  for (int i=0; i<16; i++)
  {
		cout << sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]] << "\t";
  }
	std::cout << std::endl;
#endif

  isFeatureOutterLoop:for (int i=0; i<16; i++)
  {
    FastDetectorisFeature_label2:for (int streak_size = 3; streak_size<=6; streak_size++)
    {
      // check that streak event is larger than neighbor
      if ((sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]]) <  (sae_[pol][pix_x+circle3_[(i-1+16)%16][0]][pix_y+circle3_[(i-1+16)%16][1]]))
        continue;

      // check that streak event is larger than neighbor
      if (sae_[pol][pix_x+circle3_[(i+streak_size-1)%16][0]][pix_y+circle3_[(i+streak_size-1)%16][1]] < sae_[pol][pix_x+circle3_[(i+streak_size)%16][0]][pix_y+circle3_[(i+streak_size)%16][1]])
        continue;

      // find the smallest timestamp in corner min_t
      double min_t = sae_[pol][pix_x+circle3_[i][0]][pix_y+circle3_[i][1]];
      FastDetectorisFeature_label1:for (int j=1; j<streak_size; j++)
      {
        const double tj = sae_[pol][pix_x+circle3_[(i+j)%16][0]][pix_y+circle3_[(i+j)%16][1]];
        if (tj < min_t)
          min_t = tj;
      }

      //check if corner timestamp is higher than corner
      bool did_break = false;
      FastDetectorisFeature_label0:for (int j=streak_size; j<16; j++)
      {
        const double tj = sae_[pol][pix_x+circle3_[(i+j)%16][0]][pix_y+circle3_[(i+j)%16][1]];

        if (tj >= min_t)
        {
          did_break = true;
          break;
        }
      }

      if (!did_break)
      {
        *found_streak = true;
#if DEBUG
		cout << "Inner Corner position is : " << i << " and streak size is: " << streak_size << endl;
#endif
        break;
      }
    }

	  if (*found_streak)
	  {
		  break;
	  }
  }

  if (*found_streak)
  {
    *found_streak = false;

#if DEBUG
	std::cout << "Idx Outer Data SW is: " << std::endl;
	for (int i=0; i<20; i++)
	{
	cout << sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]] << "\t";
	}
	std::cout << std::endl;
#endif

    FastDetectorisFeature_label6:for (int i=0; i<20; i++)
    {
      FastDetectorisFeature_label5:for (int streak_size = 4; streak_size<=8; streak_size++)
      {
        // check that first event is larger than neighbor
        if (sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]] <  sae_[pol][pix_x+circle4_[(i-1+20)%20][0]][pix_y+circle4_[(i-1+20)%20][1]])
          continue;

        // check that streak event is larger than neighbor
        if (sae_[pol][pix_x+circle4_[(i+streak_size-1)%20][0]][pix_y+circle4_[(i+streak_size-1)%20][1]] < sae_[pol][pix_x+circle4_[(i+streak_size)%20][0]][pix_y+circle4_[(i+streak_size)%20][1]])
          continue;

        double min_t = sae_[pol][pix_x+circle4_[i][0]][pix_y+circle4_[i][1]];
        FastDetectorisFeature_label4:for (int j=1; j<streak_size; j++)
        {
          const double tj = sae_[pol][pix_x+circle4_[(i+j)%20][0]][pix_y+circle4_[(i+j)%20][1]];
          if (tj < min_t)
            min_t = tj;
        }

        bool did_break = false;
        FastDetectorisFeature_label3:for (int j=streak_size; j<20; j++)
        {
          const double tj = sae_[pol][pix_x+circle4_[(i+j)%20][0]][pix_y+circle4_[(i+j)%20][1]];
          if (tj >= min_t)
          {
            did_break = true;
            break;
          }
        }

        if (!did_break)
        {
          *found_streak = true;
		  cout << "Outer Corner position is : " << i << " and streak size is: " << streak_size << endl;
		  break;
        }
      }
      if (*found_streak)
      {
        break;
      }
    }

  }

  //return *found_streak;
}

void parseEventsSW(uint64_t * dataStream, int32_t eventsArraySize, uint32_t *eventSlice)
{
	for (int i = 0; i < eventsArraySize; i++)
	{
		uint64_t tmp = *dataStream++;
		uint32_t x = ((tmp) >> POLARITY_X_ADDR_SHIFT) & POLARITY_X_ADDR_MASK;
		uint32_t y = ((tmp) >> POLARITY_Y_ADDR_SHIFT) & POLARITY_Y_ADDR_MASK;
		bool pol  = ((tmp) >> POLARITY_SHIFT) & POLARITY_MASK;
		int ts = tmp >> 32;
#if DEBUG
		cout << "x : " << x << endl;
		cout << "y : " << y << endl;
		cout << "ts : " << ts << endl;
#endif

		bool isCorner = 0;

		FastDetectorisFeature(x, y, ts, pol, &isCorner);

		x = 239 - x;
		y = 179 - y;

		ap_uint<32> tmpOutput = (0 << 31) + (y << 22) + (x << 12)  + (pol << 11) + isCorner;

		ap_uint<32> output;

		// Changed to small endian mode to send it to jAER
		output.range(7,0) = tmpOutput.range(31,24);
		output.range(15,8) = tmpOutput.range(23,16);
		output.range(23,16) = tmpOutput.range(15,8);
		output.range(31,24) = tmpOutput.range(7,0);

		*eventSlice++ = output.to_uint();
	}
}

int main ()
 {
	int testTimes = TEST_TIMES;

    int total_err_cnt = 0;
	int retval=0;
	/******************* Test parseEvents module from random value**************************/
	int32_t eventCnt = 5000;
	uint8_t x[eventCnt], y[eventCnt];
	uint64_t ts[eventCnt];
	bool pol[eventCnt];
	uint64_t data[eventCnt];
	uint32_t eventSlice[eventCnt], eventSliceSW[eventCnt];

	testTimes = 10;

	for(int k = 0; k < testTimes; k++)
	{
		cout << "Test " << k << ":" << endl;

	    int err_cnt = 0;

		for (int i = 0; i < eventCnt; i++)
		{
			ts[i]  = rand();
		}
		sort(ts, ts+eventCnt);

		for (int i = 0; i < eventCnt; i++)
		{
			x[i] = rand()%240;
			y[i] = rand()%180;
			pol[i] = 1;
//			idx = rand()%3;
	//		x = 255;
	//		y = 240;
//			cout << "x : " << x << endl;
//			cout << "y : " << y << endl;
//			cout << "idx : " << idx << endl;

			data[i] = (uint64_t)(ts[i] << 32) + (uint64_t)(x[i] << 17) + (uint64_t)(y[i] << 2) + (pol[i] << 1);
//			cout << "data[" << i << "] is: "<< hex << data[i]  << endl;
		}

		if (k == 17)
		{
			int tmp = 0;
		}
		parseEventsSW(data, eventCnt, eventSliceSW);
		parseEventsHW(data, eventCnt, eventSlice);

		for (int j = 0; j < eventCnt; j++)
		{
			if (eventSlice[j] != eventSliceSW[j])
			{
				std::cout << "eventSliceSW is: " << eventSliceSW[j] << std::endl;
				std::cout << "eventSlice is: " << eventSlice[j] << std::endl;

				cout << "j : " << j << endl;
				cout << "x : " << int(x[j]) << endl;
				cout << "y : " << int(y[j]) << endl;
				cout << "ts : " << ts[j] << endl;
				cout << "pol : " << pol[j] << endl;

				err_cnt++;
				cout << "Mismatch detected on TEST " << k << " and the mismatch index is: " << j << endl;
			}
		}

		if(err_cnt == 0)
		{
			cout << "Test " << k << " passed." << endl;
		}
		else
		{
			cout << "Test " << k << " failed!!!" << endl;
		}
		total_err_cnt += err_cnt;
		cout << endl;
	}

//	/******************* Test FastCheckOuterCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 15000;
//
//	// The raw data for SW and HW are exactly the same, except the data type.
//	uint8_t x, y;
//	uint32_t ts[testTimes];
//	bool pol;
//	ap_uint<2> stage[2];
//	stage[0] = ap_uint<2>(0);
//	stage[1] = ap_uint<2>(1);
//
////	uint8_t outputIdxSW[OUTER_SIZE];
////	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	bool isOuterCornerSW = 0;
//	ap_uint<1>  isOuterCornerHW = 0;
//
//	uint8_t size = OUTER_SIZE;
//
//	for (int i = 0; i < testTimes; i++)
//	{
//		ts[i]  = rand();
//	}
//	sort(ts, ts+testTimes);
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//// 	    cout << "\nArray after sorting using "
//// 	         "default sort is : \n";
//// 	    for (int i = 0; i < eventCnt; ++i)
//// 	        cout << ts[i] << " ";
//
//			x = rand()%220 + 10;
//			y = rand()%110 + 10;
////			idx = rand()%3;
//	//		x = 255;
//	//		y = 240;
//			cout << "x : " << (int)x << endl;
//			cout << "y : " << (int)y << endl;
//			cout << "ts : " << ts[k] << endl;
//
//		if (k == 1235)
//		{
//			int tmp = 0;
//		}
//
//		fastCornerHW(x, y, ts[k], &isOuterCornerHW);
//		FastDetectorisFeature(x, y, ts[k], pol, &isOuterCornerSW);
//
//		cout << "isCornerSW is: " << isOuterCornerSW << endl;
//		cout << "isCornerHW is: " << isOuterCornerHW << endl;
//
//		if (isOuterCornerSW != isOuterCornerHW.to_bool())
//		{
//			err_cnt++;
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		else
//		{
//			cout << "Test " << k << " failed!!!" << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test FastCheckOuterCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 10000;
//
//	// The raw data for SW and HW are exactly the same, except the data type.
//	uint8_t x, y;
//	uint32_t ts[testTimes];
//	bool pol;
////	uint8_t outputIdxSW[OUTER_SIZE];
////	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	bool isOuterCornerSW = 0;
//	ap_uint<1>  isOuterCornerHW = 0;
//
//	uint8_t size = OUTER_SIZE;
//
//	for (int i = 0; i < testTimes; i++)
//	{
//		ts[i]  = rand();
//	}
//	sort(ts, ts+testTimes);
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//// 	    cout << "\nArray after sorting using "
//// 	         "default sort is : \n";
//// 	    for (int i = 0; i < eventCnt; ++i)
//// 	        cout << ts[i] << " ";
//
//			x = rand()%220 + 10;
//			y = rand()%110 + 10;
////			idx = rand()%3;
//	//		x = 255;
//	//		y = 240;
//			cout << "x : " << (int)x << endl;
//			cout << "y : " << (int)y << endl;
//			cout << "ts : " << ts[k] << endl;
//
//	    FastDetectorisOuterFeature(x, y, ts[k], pol, &isOuterCornerSW);
//		fastCornerOuterHW(x, y, ts[k], 1, &isOuterCornerHW);
//
//		cout << "isCornerSW is: " << isOuterCornerSW << endl;
//		cout << "isCornerHW is: " << isOuterCornerHW << endl;
//
//		if (isOuterCornerSW != isOuterCornerHW.to_bool())
//		{
//			err_cnt++;
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		else
//		{
//			cout << "Test " << k << " failed!!!" << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test FastCheckInnerCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 10000;
//
//	// The raw data for SW and HW are exactly the same, except the data type.
//	uint32_t x, y;
//	uint32_t ts[testTimes];
//	bool pol;
////	uint8_t outputIdxSW[OUTER_SIZE];
////	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	bool isInnerCornerSW = 0;
//	ap_uint<1>  isInnerCornerHW = 0;
//
//	uint8_t size = INNER_SIZE;
//
//	for (int i = 0; i < testTimes; i++)
//	{
//		ts[i]  = rand();
//	}
//	sort(ts, ts+testTimes);
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//// 	    cout << "\nArray after sorting using "
//// 	         "default sort is : \n";
//// 	    for (int i = 0; i < eventCnt; ++i)
//// 	        cout << ts[i] << " ";
//
//			x = rand()%220 + 10;
//			y = rand()%110 + 10;
////			idx = rand()%3;
//	//		x = 255;
//	//		y = 240;
//			cout << "x : " << x << endl;
//			cout << "y : " << y << endl;
//			cout << "ts : " << ts[k] << endl;
//
//		FastDetectorisInnerFeature(x, y, ts[k], pol, &isInnerCornerSW);
//		fastCornerInnerHW(x, y, ts[k], 0, &isInnerCornerHW);
//
//		cout << "isCornerSW is: " << isInnerCornerSW << endl;
//		cout << "isCornerHW is: " << isInnerCornerHW << endl;
//
//		if (isInnerCornerSW != isInnerCornerHW.to_bool())
//		{
//			err_cnt++;
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		else
//		{
//			cout << "Test " << k << " failed!!!" << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test testFromTsDataCheckOuterCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 1000;
//
//	// The raw data for SW and HW are exactly the same, except the data type.
//	uint32_t testRawDataSW[OUTER_SIZE];
//	ap_uint<32> testRawDataHW[OUTER_SIZE];
////	uint8_t outputIdxSW[OUTER_SIZE];
////	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	ap_uint<1> isCornerSW = 0, isCornerHW = 0;
//
//	uint8_t size = OUTER_SIZE;
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//		for (int i = 0; i < size; i++)
//		{
//			testRawDataSW[i]  = rand();
//			testRawDataHW[i] = testRawDataSW[i];
//		}
//
//		// The following pattern is only used to test the boundary behavior.
//		// On normal test condition, uncomment them.
//// 		testRawDataSW[0] = 100000;
//// 		testRawDataSW[15] = 100000;
//// 		testRawDataSW[14] = 100000;
//// 		testRawDataHW[0] = 100000;
//// 		testRawDataHW[15] = 100000;
//// 		testRawDataHW[14] = 100000;
//
//		testFromTsDataCheckOuterCornerHW(testRawDataHW, size, &isCornerHW);
//		testFromTsDataCheckOuterCornerSW(testRawDataSW, size, &isCornerSW);
//
//		cout << "isCornerSW is: " << isCornerSW << endl;
//		cout << "isCornerHW is: " << isCornerHW << endl;
//
//		if (isCornerSW != isCornerHW)
//		{
//			err_cnt++;
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		else
//		{
//			cout << "Test " << k << " failed!!!" << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test testFromTsDataCheckInnerCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 1000;
//
//    // The raw data for SW and HW are exactly the same, except the data type.
//	uint32_t testRawDataSW[OUTER_SIZE];
//	ap_uint<32> testRawDataHW[OUTER_SIZE];
////	uint8_t outputIdxSW[OUTER_SIZE];
////	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	ap_uint<1> isCornerSW = 0, isCornerHW = 0;
//
//	uint8_t size = INNER_SIZE;
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
// 		for (int i = 0; i < size; i++)
//		{
// 			testRawDataSW[i]  = rand();
// 			testRawDataHW[i] = testRawDataSW[i];
//		}
//
// 		// The following pattern is only used to test the boundary behavior.
// 		// On normal test condition, uncomment them.
//// 		testRawDataSW[0] = 100000;
//// 		testRawDataSW[15] = 100000;
//// 		testRawDataSW[14] = 100000;
//// 		testRawDataHW[0] = 100000;
//// 		testRawDataHW[15] = 100000;
//// 		testRawDataHW[14] = 100000;
//
// 		testFromTsDataCheckInnerCornerSW(testRawDataSW, size, &isCornerSW);
// 		testFromTsDataToInnerCornerHW(testRawDataHW, size, &isCornerHW);
// 		// 		testFromTsDataCheckInnerCornerHW(testRawDataHW, size, &isCornerHW);
//
//		cout << "isCornerSW is: " << isCornerSW << endl;
//		cout << "isCornerHW is: " << isCornerHW << endl;
//
//		if (isCornerSW != isCornerHW)
//		{
//			err_cnt++;
//		}
//
// 		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
// 		else
// 		{
//			cout << "Test " << k << " failed!!!" << endl;
// 		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test testConvertandSortedInnerIdxSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 1000;
//
//	// The raw data for SW and HW are exactly the same, except the data type.
//	uint32_t testRawDataSW[OUTER_SIZE];
//	ap_uint<32> testRawDataHW[OUTER_SIZE];
//	ap_uint<4> outputIdxBoolSW[INNER_SIZE];
//	ap_uint<4> outputIdxBoolHW[INNER_SIZE];
//
//	uint8_t size = INNER_SIZE;
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
//		for (int i = 0; i < size; i++)
//		{
//			testRawDataSW[i]  = rand();
//
//		}
//		for (int i = 0; i < size; i++)
//		{
//			for(int j = i + 1; j < size; j++)
//			if(testRawDataSW[i] == testRawDataSW[j])  // If the same, generate again.
//				testRawDataSW[j]  = rand();
//
//			testRawDataHW[i] = testRawDataSW[i];
//		}
//
//		testFromTsDataToIdxInnerBoolDataHW(testRawDataHW, size, outputIdxBoolHW);
//		testConvertandSortedInnerIdxSW(testRawDataSW, outputIdxBoolSW);
//
//		for (int  j = 0; j < INNER_SIZE; j++)
//		{
//			for (int i = 0; i < 4; i++ )
//			{
//				if (outputIdxBoolSW[j][i] != outputIdxBoolHW[j][i])
//				{
//					err_cnt++;
//				}
//			}
//		}
//
//		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		else
//		{
//			cout << "Test " << k << " failed!!!" << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test testFromTsDataToIdxData module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 1000;
//
//    // The raw data for SW and HW are exactly the same, except the data type.
//	uint32_t testRawDataSW[OUTER_SIZE];
//	ap_uint<32> testRawDataHW[OUTER_SIZE];
//	uint8_t outputIdxSW[OUTER_SIZE];
//	ap_uint<5> outputIdxHW[OUTER_SIZE];
//
//	uint8_t size = INNER_SIZE;
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
// 		for (int i = 0; i < size; i++)
//		{
// 			testRawDataSW[i]  = rand();
//
//		}
// 		for (int i = 0; i < size; i++)
//		{
//	        for(int j = i + 1; j < size; j++)
//	        if(testRawDataSW[i] == testRawDataSW[j])  // If the same, generate again.
//	        	testRawDataSW[j]  = rand();
//
// 			testRawDataHW[i] = testRawDataSW[i];
//		}
//
// 		testConvertandSortedIdxSW(testRawDataSW, size, outputIdxSW);
// 		testFromTsDataToIdxDataHW(testRawDataHW, size, outputIdxHW);
//
//		for (int  j = 0; j < size; j++)
//		{
//			if (size == INNER_SIZE)
//			{
//				if (outputIdxSW[j] + (OUTER_SIZE - INNER_SIZE) != outputIdxHW[j])
//				{
//					err_cnt++;
//				}
//			}
//			else
//			{
//				if (outputIdxSW[j] != outputIdxHW[j])
//				{
//					err_cnt++;
//				}
//			}
//		}
//
// 		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
// 		else
// 		{
//			cout << "Test " << k << " failed!!!" << endl;
// 		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test testFromTsDataCheckInnerCornerSW module from random value**************************/
////	srand((unsigned)time(NULL));
//	testTimes = 1000;
//
//    // The raw data for SW and HW are exactly the same, except the data type.
//	int idxDataSW[OUTER_SIZE];
//	ap_uint<5> idxDataHW[OUTER_SIZE];
//	ap_uint<5> outputIdxBoolSW[OUTER_SIZE];
//	ap_uint<4> outputIdxBoolHW[OUTER_SIZE];
//
//	ap_uint<1> isCornerSW = 0, isCornerHW = 0;
//
//	uint8_t size = 18;
//
//	for(int k = 0; k < testTimes; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
// 		for (int i = 0; i < size; i++)
//		{
// 			idxDataSW[i]  = rand()%size;
// 			idxDataHW[i] = idxDataSW[i];
//		}
//
// 		testIdxDataToIdxInnerBoolDataSW(idxDataSW, size, outputIdxBoolSW);
// 		testIdxDataToIdxInnerBoolDataHW(idxDataHW, size, outputIdxBoolHW);
//
//		for (int  j = 0; j < 18; j++)
//		{
//			for (int i = 0; i < 4; i++ )
//			{
//				if (outputIdxBoolSW[j][i] != outputIdxBoolHW[j][i])
//				{
//					err_cnt++;
//				}
//			}
//		}
//
// 		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
// 		else
// 		{
//			cout << "Test " << k << " failed!!!" << endl;
// 		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

//	/******************* Test rwSAE module from random value**************************/
//	srand((unsigned)time(NULL));
//	int16_t eventCnt = 500;
//	ap_uint<TS_TYPE_BIT_WIDTH> outputDataSW[OUTER_SIZE], outputDataHW[OUTER_SIZE];
//	ap_uint<5> sizeSW, sizeHW;
//
//	uint32_t x, y;
//	uint32_t ts[eventCnt];
//
//	for(int k = 0; k < TEST_TIMES; k++)
//	{
//		cout << "Test " << k << ":" << endl;
//
//		int err_cnt = 0;
//
// 		for (int i = 0; i < eventCnt; i++)
//		{
// 			ts[i]  = rand();
//		}
// 	    sort(ts, ts+eventCnt);
//
//// 	    cout << "\nArray after sorting using "
//// 	         "default sort is : \n";
//// 	    for (int i = 0; i < eventCnt; ++i)
//// 	        cout << ts[i] << " ";
//
// 		for (int i = 0; i < eventCnt; i++)
//		{
//			x = rand()%50 + 20;
//			y = rand()%50 + 20;
////			idx = rand()%3;
//	//		x = 255;
//	//		y = 240;
////			cout << "x : " << x << endl;
////			cout << "y : " << y << endl;
////			cout << "idx : " << idx << endl;
//
////			data[i] = (uint64_t)(x << 17) + (uint64_t)(y << 2) + (1 << 1);
////			cout << "data[" << i << "] is: "<< hex << data[i]  << endl;
//	 		rwSAESW(x, y, ts[i], 1, outputDataSW, &sizeSW);
//	 		testRwSAEHW(x, y, ts[i], 1, outputDataHW, &sizeHW);
//
//	 		for (int  j = 0; j < OUTER_SIZE; j++)
//	 		{
//	 			if (outputDataSW[j] != outputDataHW[j])
//	 			{
// 	 				err_cnt++;
//	 			}
//	 		}
//		}
//
//
// 		if(err_cnt == 0)
//		{
//			cout << "Test " << k << " passed." << endl;
//		}
//		total_err_cnt += err_cnt;
//		cout << endl;
//	}

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
