#include "wienerStorm.hxx"
#define CSTACK_DEFNS
#define HAVE_UINTPTR_T
#include <cstdint>
#include <iostream>



void fitPSF(DataParams &params, MultiArray<2, double> &ps) {
    double minval=9999999, maxval = 0;
    int size = 2*ps.shape(0)*ps.shape(1);
    std::vector<double> data2(2*ps.shape(0)*ps.shape(1));
    for(int i=0, counter = 0;i<ps.shape(0);++i) {
        for(int j=0;j<ps.shape(1);++j,++counter) {
            //std::cout<<i<<" "<<j<<" "<<counter<<std::endl;
            data2[2*counter] = std::sqrt(std::abs(std::pow(i-ps.shape(0)/2,2)+std::pow(j-ps.shape(1)/2,2)));
            data2[2*counter+1] = ps(i,j);
            if (ps(i,j)< minval){minval = ps(i,j);}
            if (ps(i,j)> maxval){maxval = ps(i,j);}
        }
    }
    double sigma = 2.0, scale = maxval - minval, offset = minval;
    //std::cout<<sigma<<" "<<scale<<" "<<offset<<std::endl;
    fitGaussian(&data2[0], size/2, sigma, scale, offset);
    params.setSigma(sigma);

}



double fitPSF2D(MultiArray<2, double> &ps, double &sigma) {
    double minval=9999999, maxval = 0;
	int w = ps.shape(0), h = ps.shape(1);
    int size = ps.shape(0)*ps.shape(1);
    std::vector<double> data2(3*ps.shape(0)*ps.shape(1));
    //std::ofstream outputRoi;
	//outputRoi.open("c:\\tmp\\outputRoi.txt");
	for(int i=0, counter = 0;i<ps.shape(0);++i) {
        for(int j=0;j<ps.shape(1);++j,++counter) {
            //std::cout<<i<<" "<<j<<" "<<counter<<std::endl;
		//	outputRoi <<i<<" "<<j<<" "<<ps(i,j)<<std::endl;
            data2[3*counter] = i-ps.shape(0)/2;
			data2[3*counter+1] = j-ps.shape(1)/2;
            data2[3*counter+2] = ps(i,j);
            if (ps(i,j)< minval){minval = ps(i,j);}
            if (ps(i,j)> maxval){maxval = ps(i,j);}
        }
    }
	//outputRoi.close();
    double scale = maxval - minval, offset = minval, x0=0, y0=0;
	sigma = 2.0;
    //std::cout<<"sigma: "<<sigma<<" scale: "<<scale<<" offset:"<<offset<<std::endl;
    fitGaussian2D(&data2[0], size, sigma, scale, offset, x0, y0);

	double SSE = 0, SSD = 0, sumY = 0, sumY2 = 0;
	for(int i= w/2-2, counter = 0;i<w/2+2;++i) {
        for(int j=h/2-2;j<h/2+2;++j,++counter) {
			double xs = sq((i-ps.shape(0)/2 - x0) / sigma)+ sq((j-ps.shape(1)/2 - y0) / sigma);
            double e = std::exp(-0.5 * xs);
            double r = ps(i,j) - (scale * e + offset);
			SSE += r*r;
			sumY+= ps(i,j);
			sumY2+= sq(ps(i,j));
        }
    }
	SSD = sumY2 - sq(sumY) / size;
	double error = 1-SSE/SSD;
	sigma = std::abs(sigma);

	//std::cout<<"sigma: "<<sigma<<" scale: "<<scale<<" offset:"<<offset<< " x0: " << x0<<" y0: "<<y0<<" Error: "<<error<<std::endl;
	//std::cin.get();
	
	return error;
}
