#include "fast_detector.h"
#include <iostream>
#include <math.h>
#include <sys/socket.h> 
#include <arpa/inet.h>
#include <unistd.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "opencv2/opencv.hpp"

using namespace cv;
using namespace corner_event_detector;
//using namespace std;

int main(int argc, char* argv[])
{
  //ros::init(argc, argv, "corner_event_detector");
  char*    filename;
  bool     detector;

  filename = argv[1];

  std::ifstream file(filename);
  std::string line;

  //FastDetector xiaohung = new FastDetector();

  //creat picture for fram 
  Mat img, img_color, img_resize, img_blank, img_color_blank, image;
  img = Mat::ones(180 , 240, CV_8UC1)*127;
  img_blank = Mat::ones(180 , 240, CV_8UC1)*127;
  int imgSize = img.total() * img.elemSize();
  uchar *iptr = img.data;
  int bytes = 0;
  int key;
  int scalsz = 3;

  cvtColor(img, img_color, COLOR_GRAY2BGR);
  cvtColor(img_blank, img_color_blank, COLOR_GRAY2BGR);
  cv::resize(img_color, img_resize, cv::Size(), scalsz, scalsz);

  //make img continuos
  if ( ! img.isContinuous() ) { 
        img = img.clone();
  
          
    std::cout << line << std::endl;
        //for reading the file
      img_resize.setTo(130);
      bytes = 10000;

      // start reading from the fourth bytes. The first bytes are used to store some debug information from the server.
      for(int bufIndex = 4; bufIndex  < bytes; bufIndex = bufIndex + 4)
      {
            if(!std::getline(file, line))
            {
                file.clear();
                file.seekg(0);
            }
        std::stringstream stream(line);
        uint64_t ts;
        int x;
        int y;
        int pol;
        int OF_x;
        int OF_y;
        bool isfeature;

        stream >> ts;
        stream >> x;
        stream >> y;
        stream >> pol;
        x = 239 - x;
        y = 179 - y;
        //stream >> OF_x;
        //stream >> OF_y;
        detector = FastDetectorisFeature(x, y, ts, pol, &isfeature);
        std::cout<<detector << std::endl;
        std::cout<<line << std::endl;
        std::cout<<pol << std::endl;


        //drow gray frame and event pixel 
        if(pol == 1)
            {
                for(int i = 0; i < scalsz; i++)
                {
                    for(int j = 0; j < scalsz; j++)
                    {
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[0] = 150;
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[1] = 155;
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[2] = 155;
                    }
                }
            }
            else
            {
                for(int i = 0; i < scalsz; i++)
                {
                    for(int j = 0; j < scalsz; j++)
                    {
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[0] = 100;
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[1] = 100;
                        img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[2] = 100;
                    }
                }
            }

        if (detector == 1){
          connernum += 1;
                for(int i = 0; i < scalsz; i++)
                {
                    for(int j = 0; j < scalsz; j++)
                    {
                      img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[0] = 255;
                      img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[1] = 255;
                      img_resize.at<Vec3b>(y*scalsz+i, x*scalsz+j)[2] = 0;
                    }
                }   
        }
      //std::cout<<img_resize << std::endl;
      // cv::imshow("Event slice Client", img_resize); 
      // if (key = cv::waitKey(10) >= 0);
      }
            cv::imshow("Event slice Client", img_resize); 
      if (key = cv::waitKey(1) >= 0);
    }

  // create feature detecotr
  //detector = xiaohung.isFeature(9, 15, 200, true);
  //std::cout<<"corner number = ";
  //std::cout<<connernum;
  //delete detector;

  return 0;
}