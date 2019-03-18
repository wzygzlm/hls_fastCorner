#include "fast_detector.h"
#include "ap_int.h"
#include <iostream>
#include "ap_int.h"
#include "hls_stream.h"

static col_pix_t glPLSlices[SLICES_NUMBER][SLICE_WIDTH][SLICE_HEIGHT/COMBINED_PIXELS];
static sliceIdx_t glPLActiveSliceIdx = 0, glPLTminus1SliceIdx, glPLTminus2SliceIdx;

#define INPUT_COLS 4

static const int sensor_width_ = 240;
static const int sensor_height_ = 180;

static uint16_t eventIterSize = 100;
//static uint16_t areaEventThr = 1000;
//static uint16_t OFRetRegs[8][8];
//uint16_t areaEventThrBak = areaEventThr;

// SAE (Surface of Active Event)
ap_uint<18> sae_[1][sensor_width_][sensor_height_];

const int circle3_[16][2] = {{0, 3}, {1, 3}, {2, 2}, {3, 1},
      {3, 0}, {3, -1}, {2, -2}, {1, -3},
      {0, -3}, {-1, -3}, {-2, -2}, {-3, -1},
      {-3, 0}, {-3, 1}, {-2, 2}, {-1, 3}};
const int circle4_[20][2] = {{0, 4}, {1, 4}, {2, 3}, {3, 2},
      {4, 1}, {4, 0}, {4, -1}, {3, -2},
      {2, -3}, {1, -4}, {0, -4}, {-1, -4},
      {-2, -3}, {-3, -2}, {-4, -1}, {-4, 0},
      {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}};



void FastDetectorisFeature(int pix_x, int pix_y, int timesmp, bool polarity, bool *found_streak)
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
        break;
      }


    if (*found_streak)
    {
      break;
    }

  }

  if (*found_streak)
  {
    *found_streak = false;
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
}


void fastCorner(hls::stream<uint8_t>  &xInStream, hls::stream<uint8_t> &yInStream, hls::stream<uint32_t> &tsInStream,
		hls::stream<apUint15_t> &minStream,  hls::stream<apUint6_t> &OFStream)
{
	ap_uint<8> x, y;
	x = xInStream.read();
	y = yInStream.read();
	uint32_t ts = tsInStream.read();

	bool polarity = true;
	bool foundFlg = false;
	bool *found_streak = &foundFlg;

	FastDetectorisFeature(x, y, ts, polarity, found_streak);

	minStream.write((apUint15_t)found_streak);
	OFStream.write((apUint6_t)found_streak);

}

// Function Description: return the minimum value of an array.
ap_int<16> min(ap_int<16> inArr[2*SEARCH_DISTANCE + 1], int8_t *index)
{
#pragma HLS INLINE
#pragma HLS PIPELINE
#pragma HLS ARRAY_RESHAPE variable=inArr complete dim=1
	ap_int<16> tmp = inArr[0];
	int8_t tmpIdx = 0;
	minLoop: for(int8_t i = 0; i < 2*SEARCH_DISTANCE + 1; i++)
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


void getXandY(const uint64_t * data, hls::stream<uint8_t>  &xStream, hls::stream<uint8_t> &yStream, hls::stream<uint32_t> &tsStream, hls::stream<apUint17_t> &packetEventDataStream)
//void getXandY(const uint64_t * data, int32_t eventsArraySize, ap_uint<8> *xStream, ap_uint<8> *yStream)
{
#pragma HLS INLINE off

	// Every event always consists of 2 int32_t which is 8bytes.
//	getXandYLoop:for(int32_t i = 0; i < eventIterSize; i++)
//	{
		uint64_t tmp = *data;
		ap_uint<8> xWr, yWr;
		xWr = ((tmp) >> POLARITY_X_ADDR_SHIFT) & POLARITY_X_ADDR_MASK;
		yWr = ((tmp) >> POLARITY_Y_ADDR_SHIFT) & POLARITY_Y_ADDR_MASK;
		bool pol  = ((tmp) >> POLARITY_SHIFT) & POLARITY_MASK;
		uint32_t ts = tmp >> 32;

//		writePix(xWr, yWr, glPLActiveSliceIdx);
//		resetPix(xWr, yWr, glPLActiveSliceIdx + 3);

//		shiftCnt = 0;
//		miniRetVal = 0x7fff;
//		for(int8_t i = 0; i <= 2*SEARCH_DISTANCE; i++)
//		{
//				miniSumTmp[i] = 0;
//		}
//		for(int8_t i = 0; i <= 2*SEARCH_DISTANCE; i++)
//		{
//			for(int8_t j = 0; j <= 2*SEARCH_DISTANCE; j++)
//			{
//				localSumReg[i][j] = 0;
//			}
//		}

		xStream << xWr;
		yStream << yWr;
		tsStream << ts;
		packetEventDataStream << apUint17_t(xWr.to_int() + (yWr.to_int() << 8) + (pol << 16));
//		*xStream++ = xWr;
//		*yStream++ = yWr;
//	}
}


void feedbackWrapperAndOutputResult(hls::stream<apUint15_t> &miniSumStream, hls::stream<apUint6_t> &OFRetStream,
						hls::stream<apUint17_t> &packetEventDataStream,
					 uint32_t *eventSlice)
{
	apUint17_t tmp1 = packetEventDataStream.read();
	apUint15_t tmpMiniSumRet = miniSumStream.read();
	ap_int<9> tmp2 = tmpMiniSumRet.range(8, 0);
	apUint6_t tmpOF = OFRetStream.read();

	ap_uint<32> output = (tmp2, (tmpOF, tmp1));
//		std :: cout << "tmp1 is "  << std::hex << tmp1 << std :: endl;
//		std :: cout << "tmp2 is "  << std::hex << tmp2 << std :: endl;
//		std :: cout << "output is "  << std::hex << output << std :: endl;
//		std :: cout << "eventSlice is "  << std::hex << output.to_int() << std :: endl;
	*eventSlice++ = output.to_uint();
}

#pragma SDS data access_pattern(data:SEQUENTIAL, eventSlice:SEQUENTIAL)
// #pragma SDS data data_mover(data:AXIFIFO:1, eventSlice:AXIFIFO:2)
// #pragma SDS data buffer_depth(data:512, eventSlice:1)
#pragma SDS data data_mover(data:AXIDMA_SIMPLE:1, eventSlice:AXIDMA_SIMPLE:2)
#pragma SDS data copy(data[0:eventsArraySize], eventSlice[0:eventsArraySize])
#pragma SDS data mem_attribute(data:PHYSICAL_CONTIGUOUS, eventSlice:PHYSICAL_CONTIGUOUS)
// #pragma SDS data zero_copy(eventSlice[0:DVS_WIDTH * DVS_HEIGHT])
//#pragma SDS data sys_port(data:AFI, eventSlice:AFI)
void parseEvents(const uint64_t * data, int32_t eventsArraySize, uint32_t *eventSlice, ap_uint<1> *outLed)
{
//#pragma HLS INTERFACE axis port=data
//#pragma HLS INTERFACE axis register both port=eventSlice

//#pragma HLS INTERFACE m_axi port=eventSlice offset=slave bundle=gmem max_read_burst_length=2 max_write_burst_length=256

//#pragma HLS INTERFACE s_axilite port=return bundle=control
	hls::stream<uint8_t>  xInStream("xInStream"), yInStream("yInStream");
	hls::stream<uint8_t>  xOutStream("xOutStream"), yOutStream("yOutStream");
	hls::stream<uint32_t>  tsInStream("tsInStream");

	hls::stream<sliceIdx_t> idxStream("idxStream");
	hls::stream<apUint17_t> pktEventDataStream("EventStream");
	hls::stream<apUint6_t> OFRetStream("OFStream");
	hls::stream<apIntBlockCol_t> refStream("refStream"), tagStreamIn("tagStream");
	hls::stream<apUint15_t> miniSumStream("miniSumStream");

	hls::stream<uint16_t> thrStream("thresholdStream");
#pragma HLS STREAM variable=thrStream depth=3 dim=1
	hls::stream<apUint1_t> rotatFlgStream("rotationFlgStream");

	hls::stream<uint8_t>  xWrStream("xWrStream"), yWrStream("yWrStream");
	hls::stream<sliceIdx_t> idxWrStream("idxWrStream");
	hls::stream<col_pix_t> currentColStream("currentColStream");

	hls::stream<apUint112_t> outStream("sumStream");
#pragma HLS STREAM variable=outStream depth=2 dim=1
#pragma HLS RESOURCE variable=outStream core=FIFO_SRL
	hls::stream<int16_t> outSumStream("outSumStream");
	hls::stream<int8_t> OF_yStream("OF_yStream");

	hls::stream<apUint6_t> refZeroCntStream("refZeroCntStream");
#pragma HLS STREAM variable=refZeroCntStream depth=2 dim=1
	hls::stream<uint16_t> refZeroCntSumStream("refZeroCntSumStream");

	hls::stream<apUint42_t> tagColValidCntStream("tagColValidCntStream");
#pragma HLS STREAM variable=tagColValidCntStream depth=2 dim=1
	hls::stream<uint16_t> tagColValidCntSumStream("tagColValidCntSumStream");

	hls::stream<apUint42_t> refTagValidCntStream("tagColValidCntStream");
#pragma HLS STREAM variable=refTagValidCntStream depth=2 dim=1
	hls::stream<uint16_t> refTagValidCntSumStream("tagColValidCntSumStream");

	eventIterSize = eventsArraySize;

	parseEventsLoop:for(int32_t i = 0; i < eventIterSize; i++)
	{
#pragma HLS LOOP_TRIPCOUNT min=1 max=10000
		DFRegion:
		{
#pragma HLS DATAFLOW
#pragma HLS STREAM variable=pktEventDataStream depth=2 dim=1
#pragma HLS RESOURCE variable=pktEventDataStream core=FIFO_SRL
#pragma HLS STREAM variable=miniSumStream depth=2 dim=1
#pragma HLS RESOURCE variable=miniSumStream core=FIFO_SRL
#pragma HLS STREAM variable=tagStreamIn depth=6 dim=1
#pragma HLS RESOURCE variable=tagStreamIn core=FIFO_SRL
#pragma HLS STREAM variable=refStream depth=2 dim=1
#pragma HLS RESOURCE variable=refStream core=FIFO_SRL
			// This one has wrong block sad sum module.
//			getXandY(dataStream++, xInStream, yInStream, pktEventDataStream);
//			rotateSlice(xInStream, yInStream, xOutStream, yOutStream, idxStream);
//			rwSlices(xOutStream, yOutStream, idxStream, refStream, tagStreamIn);
//			miniSADSumWrapper(refStream, tagStreamIn, miniSumStream, OFRetStream);
//			outputResult(miniSumStream, OFRetStream, pktEventDataStream, eventSlice++);

//			getXandY(dataStream++, xInStream, yInStream, pktEventDataStream);
//			rotateSliceNoRotationFlg(xInStream, yInStream, xOutStream, yOutStream, idxStream);
//			rwSlices(xOutStream, yOutStream, idxStream, refStream, tagStreamIn);
//			colStreamToColSum(refStream, tagStreamIn, outStream);
//			accumulateStream(outStream, outSumStream, OF_yStream);
//			findStreamMin(outSumStream, OF_yStream, miniSumStream, OFRetStream);
//			outputResult(miniSumStream, OFRetStream, pktEventDataStream, eventSlice++);

			// With feedback
			getXandY(data++, xInStream, yInStream, tsInStream, pktEventDataStream);
			fastCorner(xInStream, yInStream, tsInStream, miniSumStream, OFRetStream);
			//rotateSlice(xInStream, yInStream, tsInStream, thrStream, xOutStream, yOutStream, idxStream);
			//rwSlices(xOutStream, yOutStream, idxStream, refStream, tagStreamIn);
			//colStreamToColSum(refStream, tagStreamIn, outStream, refZeroCntStream, tagColValidCntStream, refTagValidCntStream);
			//accumulateStream(outStream, outSumStream, OF_yStream, refZeroCntStream, tagColValidCntStream,  refTagValidCntStream);
			//findStreamMin(outSumStream, OF_yStream, miniSumStream, OFRetStream);
			feedbackWrapperAndOutputResult(miniSumStream, OFRetStream, pktEventDataStream, eventSlice++);

			// This is the version combined rwSlices and colStreamToColSum together
			// It consumes less resources but has higher II.
//			getXandY(dataStream++, xInStream, yInStream, pktEventDataStream);
//			rotateSlice(xInStream, yInStream, xOutStream, yOutStream, idxStream);
//			rwSlicesAndColStreams(xOutStream, yOutStream, idxStream, outStream);
//			accumulateStream(outStream, outSumStream, OF_yStream);
//			findStreamMin(outSumStream, OF_yStream, miniSumStream, OFRetStream);
//			outputResult(miniSumStream, OFRetStream, pktEventDataStream, eventSlice++);

// read and write array in separate process function is not supported in dataflow.
//			getXandY(dataStream++, xInStream, yInStream, pktEventDataStream);
//			rotateSlice(xInStream, yInStream, xOutStream, yOutStream, idxStream);
//			readSlices( xOutStream, yOutStream, idxStream, xWrStream, yWrStream, idxWrStream,currentColStream, refStream, tagStreamIn);
//			writeSlices(xWrStream, yWrStream, idxWrStream, currentColStream);
//			miniSADSumWrapper(refStream, tagStreamIn, miniSumStream, OFRetStream);
//			outputResult(miniSumStream, OFRetStream, pktEventDataStream, eventSlice++);

//			getXandY(dataStream, xInStream, yInStream, pktEventDataStream);
//			rotateSlice(xInStream, yInStream, thrStream, xOutStream, yOutStream, idxStream, rotatFlgStream);
//			rwSlices(xOutStream, yOutStream, idxStream, refStream, tagStreamIn);
//			miniSADSumWrapper(refStream, tagStreamIn, miniSumStream, OFRetStream);
//			feedbackWrapperAndOutputResult(miniSumStream, OFRetStream, pktEventDataStream, rotatFlgStream, thrStream, eventSlice);
		}
	}
}
