#include "wienerStorm.hxx"
#define CSTACK_DEFNS
#define HAVE_UINTPTR_T
#include <cstdint>



void fitPSF(DataParams &params, MultiArray<2, double> &ps) {
    double minval=9999999, maxval = 0;
    int size = 2*ps.shape(0)*ps.shape(1);
    std::vector<double> data2(2*ps.shape(0)*ps.shape(1));
    for(int i=0, counter = 0;i<ps.shape(0);++i) {
        for(int j=0;j<ps.shape(1);++j,++counter) {
            //std::cout<<i<<" "<<j<" "<<counter <<std::endl;
            data2[2*counter] = std::sqrt(std::abs(std::pow(i-ps.shape(0)/2,2)+std::pow(j-ps.shape(1)/2,2)));
            data2[2*counter+1] = ps(i,j);
            if (ps(i,j)< minval){minval = ps(i,j);}
            if (ps(i,j)> maxval){maxval = ps(i,j);}
        }
    }
    double sigma = 2.0, scale = maxval - minval, offset = minval;
    std::cout<<sigma<<" "<<scale<<" "<<offset<<std::endl;
    fitGaussian(&data2[0], size/2, sigma, scale, offset);
    params.setSigma(sigma);

}


