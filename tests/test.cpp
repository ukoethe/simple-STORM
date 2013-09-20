#include <iostream>

#include <clocale>
#include "wienerStorm.hxx"

#include "dataparams.cpp"
#include <string>
#include <random>

#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WIN32
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
 #endif



class CliProgressFunctor : public ProgressFunctor
{
public:
    CliProgressFunctor() : ProgressFunctor(), m_stacksize(0), m_frame(0) {};
    ~CliProgressFunctor() {};
    virtual void setStage(WienerStormStage stage)
    {
        m_stage = stage;
        switch (m_stage) {
            case CameraParameters:
                std::cout << std::endl << "Estimating camera gain and offset..." << std::endl;
                break;
            case PSFWidth:
                std::cout << std::endl << "Estimating PSF width..." << std::endl;
                break;
            case ParameterCheck:
                std::cout << std::endl << "Checking parameter..." << std::endl;
                break;
            case Localization:
                std::cout << std::endl << "Localizing molecules..." << std::endl;
                break;
        }
    }

    virtual void setStackSize(int stacksize)
    {
        m_stacksize = stacksize;
        m_frame = 0;
        helper::progress(-1, -1);
    }

    virtual void frameFinished(int frame)
    {
        ++m_frame;
#ifdef OPENMP_FOUND
        if(omp_get_thread_num()==0) { // master thread
            helper::progress(m_frame, m_stacksize);
        }
#else
        helper::progress(m_frame, m_stacksize);
#endif
    }
private:
    WienerStormStage m_stage;
    int m_stacksize;
    std::atomic<int> m_frame;
};

class Counter {
public:
  int suc, fail, tot;
  Counter()
  : suc(0), fail(0), tot(0)
  {}
  void update(bool status, std::string message){
    tot += 1;
    if (status) {
      suc += 1;
      std::cout<<message<<" "<<"succeeded"<<std::endl;
    }
    else {
      fail += 1;
      std::cout<<message<<" "<<"failed"<<std::endl;
    }
  }
};

 #ifdef HDF5_FOUND
void checkHDF5Import(char *argv[], Counter &counter){
    std::string p = argv[0];
    size_t pos = p.find_last_of("/\\");
    const char * c = p.substr(0,pos).append("\test.hdf5").c_str();
    vigra::HDF5File* info = new vigra::HDF5File(c, HDF5File::Open);
    ArrayVector<hsize_t> shape2 = info->getDatasetShape("/data");
    Shape3 shape(shape2[0],shape2[1],shape2[2]);
    vigra::MultiArray<3, float> array(shape);
    HDF5File* info2 = reinterpret_cast<HDF5File*>(info);
    info2->read("/data", array);
    int sum = 0;
    for (auto it = array.begin(); it!=array.end();++it){
      sum += *it;
    }

    std::string str = "HDF5 import test";
    counter.update(sum == 404265, str);
}
#endif

void checkSlopeFitting(Counter &counter){
	DataParams params;
	double y[199]={112.257, 130.122, 127.714, 138.698, 168.373, 210.012, 181.645, 145.412, 216.07, 210.501, 212.668, 212.866, 275.415, 223.914, 304.817, 263.768, 306.686, 404.824, 1765.64, 304.13, 347.583, 1163, 807.224, 717.69, 701.695, 479.984, 880.516, 276.231, 5952, 618.141, 2084.57, 7048.31, 2357.25, 9388.01, 1210.07, 6167.57, 468.566, 595.57, 3234.34, 1167.18, 2732.19, 706.2, 553.393, 638.334, 698.045, 4181.42, 624.212, 945.567, 6384.41, 938.5, 9279.77, 543.834, 851.112, 1725.43, 708.523, 3088.33, 746.009, 7109.13, 7425.94, 1689.41, 1098.08, 1480.27, 4223.14, 9061.06, 1128.88, 781.44, 9881.85, 992.672, 933.362, 2499.68, 1016.99, 1016.99, 883.461, 3078.25, 1330.53, 4828.74, 1016.69, 967.394, 1023.55, 828.5, 1626.89, 882.029, 984.163, 4341.89, 1749.83, 1191.99, 782.203, 1556.45, 1056.17, 1063.65, 1636.7, 1307.19, 1515.03, 1073.64, 1122.99, 2125.2, 1186.86, 1264.2, 1077.52, 1023.45, 1263.21, 1344.87, 1616.35, 1404.44, 1084.35, 1266.46, 2023.47, 2153.64, 1410.55, 2829.72, 2234.01, 1610.17, 1690.91, 1577.3, 3158.41, 1581.02, 2231.57, 2238.11, 1376.21, 3061.87, 1333.25, 1838.35, 1386.34, 1339.25, 2483.4, 2153.71, 1407.03, 2248.55, 2293.17, 1428.32, 1667.49, 2630.26, 1763.13, 2458.82, 2116.54, 1712.12, 1666.49, 1817.84, 2788.45, 1839.55, 1840.85, 2832.7, 1630.17, 2145.92, 1677.52, 1943.53, 2277.66, 2004.42, 4084.13, 1874.97, 2703.86, 3055.17, 1800.39, 3384.29, 3473.15, 2333.21, 2410.02, 2813.38, 2262.14, 2952.61, 1823.17, 4162.82, 2236.18, 2408.4, 2061.7, 2612.79, 2194.09, 2043.17, 3491.08, 3040.04, 2798.53, 2852.31, 2651.45, 2427.94, 2569.53, 2934.04, 2503.19, 3055.83, 2916.36, 2795.86, 4480.15, 2641.22, 3325.14, 3807.39, 3315.28, 2996.29, 4286.11, 3119.74, 5070.81, 4472.28, 5342.66, 3581.01, 5261.45, 4363.62, 5305.16, 4112.49, 5150.07, 4639.34, 6136.55};
	double x[199] = {404.32, 407.585, 409.15, 412.515, 413.445, 415.62, 417.93, 420.225, 422.535, 424.84, 427.145, 429.455, 431.76, 434.07, 436.375, 438.68, 440.99, 443.295, 445.605, 447.91, 450.215, 452.525, 454.83, 457.135, 459.455, 461.75, 464.06, 466.365, 468.705, 470.99, 473.31, 475.615, 477.925, 480.205, 482.53, 484.825, 487.135, 489.455, 491.75, 494.05, 496.46, 498.69, 500.98, 503.305, 505.69, 507.97, 510.255, 512.665, 514.81, 517.26, 519.44, 522.005, 524.535, 526.74, 528.675, 531.19, 533.285, 535.595, 537.895, 541.015, 542.535, 544.875, 547.325, 549.485, 551.995, 554.035, 556.395, 558.885, 561.075, 563.425, 566.095, 567.95, 570.655, 573.17, 575.82, 577.415, 579.6, 581.835, 584.82, 586.855, 588.75, 590.98, 593.28, 595.64, 597.945, 600.39, 602.57, 605.04, 607.16, 609.885, 611.93, 614.965, 619.18, 619.985, 621.65, 623.86, 625.78, 627.89, 630.485, 632.855, 635.135, 641.4, 642.8, 643.755, 644.64, 646.5, 648.86, 654.875, 654.925, 655.62, 658.14, 662.405, 662.585, 665.34, 667.32, 670.485, 671.82, 674.385, 676.435, 680.485, 682.765, 685.435, 685.625, 689.23, 690.455, 695.32, 696.18, 698.705, 702.065, 703.69, 705.345, 706.55, 709.83, 713.78, 713.97, 716.75, 718.95, 722.93, 725.45, 729.565, 730.63, 733.325, 738.35, 747.63, 754.27, 755.015, 760.94, 762.195, 765.015, 765.17, 768.865, 770.865, 773.5, 775.905, 777.7, 781.94, 784.405, 791.565, 794.19, 796.2, 814.7, 819.375, 824.435, 834.045, 836.58, 837.525, 838.375, 844.44, 857.165, 860.425, 864.55, 873.34, 874.695, 889.655, 903.6, 910.54, 911.25, 918.08, 924.205, 931.095, 954.46, 956.19, 967.28, 976.11, 1010.53, 1012.94, 1016.24, 1030.37, 1053.03, 1115.05, 1162.13, 1217.96, 1235.89, 1237.95, 1311.09, 1342.69, 1410.18, 1421.11, 1557.55};
	doRansac(params,x,y, 199);
	//std::cout<<"gain: "<<params.getSlope()<<" offset: "<<params.getIntercept()<<std::endl;
	std::string str = "Skellam fit test";
	counter.update(std::abs(params.getSlope()-4.0)<0.5 && std::abs(params.getIntercept()-370)<10, str);
	}

void checkPSFFitting(Counter &counter) {
  DataParams params;
  params.setSlope(1);
  params.setIntercept(0);
  vigra::MultiArray<2,double> img(vigra::Shape2(30,30));
  for (int i =0; i< 30; ++i){
    for (int j = 0; j < 30; ++j){
      img(i,j) = 10*std::exp(-0.5*(std::pow((i-15.)/2,2)+std::pow((j-15.)/2,2)));
    }
  }
  std::vector<double> BGVar(1);
  fitPSF(params, img);
  std::string str = "PSF fitting test";
  //std::cout<<"sigma: "<<params.getSigma()<<std::endl;
  counter.update(std::abs(params.getSigma() - 2)< 0.2, str);

}

void checkParameterCheck(Counter &counter) {
  DataParams params;
  params.setSlope(1);
  params.setIntercept(0);
  vigra::MultiArray<2,double> img(vigra::Shape2(100,100));
  std::default_random_engine generator (100);
  std::normal_distribution<double> distribution(0,10);
  for (int i =0; i< 100; ++i){
    for (int j = 0; j < 100; ++j){
      img(i,j) = distribution(generator);
    }
  }
  std::vector<double> BGVar;
  getBGVariance(params, img, BGVar, 0);
  std::string str = "Check Parameter test";
  counter.update(std::abs(BGVar[0]-10)<0.5,str);
}

void checkGetMask(Counter &counter){
  MultiArray<2,float> mask(Shape2(40,40));
  DataParams params;
  params.setAlpha(0.001f);
  params.setSigma(1.0);
  BasicImage<float> array(40,40);
  for (int i =0; i< 40; ++i){
    for (int j = 0; j < 40; ++j){
      array(i,j) = 10*std::exp(-0.5*(std::pow((i-15.)/1.0,2)+std::pow((j-15.)/1.0,2)));
    }
  }
  getMask(params, array, 0, mask);
  float sum = 0;
  for (int i = 0; i< 40; ++i){
    for (int j = 0; j < 40; ++j){
       sum += mask(i,j);
    }
  }
  std::string str = "Mask test";
  counter.update(int(sum) == 9, str);
}



template <class C>
void drawCoordsToImage2(const std::vector<std::set<C> >& coords) {
    typename std::vector<std::set<C> >::const_iterator it;
    //  loop over the images
    for(it = coords.begin(); it != coords.end(); ++it) {
        drawCoordsToImage2( *it);
    }
}

/**
 *  Draw coordinates detected in one frame into the resulting image
 */
template <class C>
void drawCoordsToImage2(const std::set<C>& coords) {
    //  loop over the coordinates
    typename std::set<C>::iterator it2;

    for(it2 = coords.begin(); it2 != coords.end(); it2++) {
        const C& c = *it2;
        std::cout<<c.x<<" "<<c.y<<""<<c.val<<std::endl;
    }
}


void checkWholeProgram(Counter &counter){
  DataParams params;
  char c[FILENAME_MAX];

  GetCurrentDir(c, sizeof(c));
  sprintf(c,"%s%s", c,"\\test2.tif");
  params.setInFile(c, false);
  params.doSanityChecks();
  int stacksize = params.shape(2);
  std::vector<std::set<Coord<float> > > res_coords(stacksize);
  CliProgressFunctor func;
  std::ofstream cf;
  if(!params.getCoordsFile().empty()) {
      cf.open(params.getCoordsFile());
      vigra_precondition(cf.is_open(), "Could not open coordinate-file for writing.");
  }
  wienerStorm(params, res_coords, func);
//   typename std::vector<std::set<float> >::const_iterator it;
//   typename std::set<float>::iterator it2;
  float sumx=0, sumy=0, sumval=0;
  int i =0;
  for(auto it = res_coords.begin(); it != res_coords.end(); ++it) {
      for(auto it2 = it->begin(); it2 != it->end(); it2++, ++i) {
          const Coord<float>& c = *it2;
          sumx+= c.x;
          sumy+= c.y;
          sumval += c.val;

      }
  }
  std::string str = "Whole Test";
  std::cout<<"sumx: "<<sumx<<" sumy: "<<sumy<<"  "<<"sumval: "<<sumval<<std::abs(sumx -1336300)<<"  "<<std::abs(sumy -1134510)<<std::endl;
  counter.update(std::abs(sumx -1336300)<100 && std::abs(sumy -1134510)<100 && std::abs(sumval -2363.74)<10, str);

}



int main(int argc, char* argv[]) {
  Counter counter;
  int successCounter = 0, failCounter = 0, testCounter = 0;

  #ifdef HDF5_FOUND
    checkHDF5Import(argv, counter);
  #endif
  checkSlopeFitting(counter);
  checkParameterCheck(counter);
  checkPSFFitting(counter);
  checkGetMask(counter);
  checkWholeProgram(counter);
  printf("successful tests: (%d/%d)\n", counter.suc, counter.tot);
  printf("failed     tests: (%d/%d)", counter.fail, counter.tot);
  std::cin.get();
  return 0;
}

