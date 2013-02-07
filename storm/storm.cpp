/************************************************************************/
/*                                                                      */
/*                  ANALYSIS OF STORM DATA                              */
/*                                                                      */
/*         Copyright 2010 by Joachim Schleicher and Ullrich Koethe      */
/*                                                                      */
/*    Please direct questions, bug reports, and contributions to        */
/*    joachim.schleicher@iwr.uni-heidelberg.de                          */
/************************************************************************/

#define CHUNKSIZE 10

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>

#define R_INTERFACE_PTRS
#define R_NO_REMAP
#include <Rembedded.h>
#include <Rinterface.h>

#include "program_options_getopt.h"
#include "wienerStorm.hxx"
#include "configVersion.hxx"

#include <vigra/impex.hxx>
#include "myimportinfo.h"
#ifdef HDF5_FOUND
    #include <vigra/hdf5impex.hxx>
#endif

#include <vigra/timing.hxx>

void preventRConsoleWrite(const char* buf, int buflen)
{}

// MAIN
int main(int argc, char** argv) {
    // Read commandline Parameters
    std::map<char, double> params;
    std::map<char, std::string> files;
    if(parseProgramOptions(argc, argv, params, files)!=0) {
        return -1;
    }

    char *Rargv[] = {"REmbeddedStorm", "--silent"};
    R_SignalHandlers = FALSE;
    Rf_initEmbeddedR(sizeof(Rargv) / sizeof(Rargv[0]), Rargv);

    int factor = (int)params['g'];
    int roilen = (int)params['m'];
    float threshold = params['t'];
    std::string infile = files['i'];
    std::string outfile = files['o'];
    std::string coordsfile = files['c'];
    std::string filterfile = files['f'];
    std::string frames = files['F'];
    char verbose = (char)params['v'];

    if(verbose) {
        std::cout << "thr:" << threshold << " factor:" << factor << std::endl;
    }

    try
    {

        MultiArray<3,float> in;
        typedef MultiArrayShape<3>::type Shape;

        MyImportInfo info(infile, argv[0]);
        //~ in.reshape(info.shape());
        //~ readVolume(info, in);
        int stacksize = info.shape()[2];
        Size2D size2 (info.shapeOfDimension(0), info.shapeOfDimension(1)); // isnt' there a slicing operator?


        if(verbose) {
            std::cout << "Images with Shape: " << info.shape() << std::endl;
            std::cout << "Processing a stack of " << stacksize << " images..." << std::endl;
        }


        // found spots. One Vector over all images in stack
        // the inner set contains all spots in the image
        std::vector<std::set<Coord<float> > > res_coords(stacksize);
        BasicImage<float> filter(info.shapeOfDimension(0), info.shapeOfDimension(1)); // filter in fourier space
        DImage res((size2-Diff2D(1,1))*factor+Diff2D(1,1));
        // check if outfile is writable, otherwise throw error -> exit
        exportImage(srcImageRange(res), ImageExportInfo(outfile.c_str()));
        if(coordsfile!="") {
            std::ofstream cf (coordsfile.c_str());
            vigra_precondition(cf.is_open(), "Could not open coordinate-file for writing.");
            cf.close();
        }

        USETICTOC;
        TIC;  // measure the time

        // STORM Algorithmus
        //printIntensities(info);

        float a,offset, intercept;
        std::vector<float> parameterTrafo(6);
        findCorrectionCoefficients(info, parameterTrafo); 	//first 3 entries are parameters for the raw signal to poisson transformation
        													//the last 3 entries are parameters for the poission to gaussian with sigma = 1 transformation
        parameterTrafo[0] = 1;
		parameterTrafo[1] = 0;
		parameterTrafo[2] = 0;
        std::cout<<"a: "<<parameterTrafo[0]<<" b: "<<parameterTrafo[1]<<" intercept: "<<parameterTrafo[2]<<std::endl;
        //showPoisson(info, parameterTrafo);

        //parameterTrafo[3] = 1;
        //parameterTrafo[4] = 0;
        //parameterTrafo[5] = 3./8.;

        MultiArray<3, float> PoissonMeans(Shape3(info.shape()[0],info.shape()[1], 1));
        getPoissonDistributions(info, parameterTrafo[0],parameterTrafo[1], PoissonMeans);

        int w= info.shapeOfDimension(0), h =info.shapeOfDimension(1);
        //int vecw[] = {20,30,40,50,60,70,80};
        //int vech[] = {10,10,10,10,10,10,10};
        //printIntensities(info, vecw, vech, 7, parameterTrafo[0],parameterTrafo[1]);

        generateFilter(info, filter, filterfile, parameterTrafo);  // use the specified one or create wiener filter from the data
        //for(int i = 0; i < info.shapeOfDimension(0);i++){
        //	for(int j = 1; j < info.shapeOfDimension(1); j++){
        //		std::cout<<filter(i,j)<<" ";
        //	}
        //}
        //int sdf;
        //std::cin>>sdf;

        wienerStorm(info, filter, res_coords,parameterTrafo, PoissonMeans, threshold, factor, roilen, frames, verbose);

        // resulting image
        drawCoordsToImage<Coord<float> >(res_coords, res);

        int numSpots = 0;
        if(coordsfile != "") {
            numSpots = saveCoordsFile(coordsfile, res_coords, info.shape(), factor);
        }

        // end: done.
        TOC;
        std::cout << "detected " << numSpots << " spots." << std::endl;

        // some maxima are very strong so we scale the image as appropriate :
        double maxlim = 0., minlim = 0;
        findMinMaxPercentile(res, 0., minlim, 0.996, maxlim);
        std::cout << "cropping output values to range [" << minlim << ", " << maxlim << "]" << std::endl;
        if(maxlim > minlim) {
            transformImage(srcImageRange(res), destImage(res), ifThenElse(Arg1()>Param(maxlim), Param(maxlim), Arg1()));
        }
        exportImage(srcImageRange(res), ImageExportInfo(outfile.c_str()));



    }
    catch (vigra::StdException & e)
    {
        std::cout<<"There was an error:"<<std::endl;
        std::cout << e.what() << std::endl;
        return 1;
    }

    // prevent printing of R warnings
    void (*ptr_R_WriteConsole_old)(const char *, int) = ptr_R_WriteConsole;
    FILE *R_Consolefile_old = R_Consolefile;
    ptr_R_WriteConsole = preventRConsoleWrite;
    R_Consolefile = NULL;
    Rf_endEmbeddedR(0);
    ptr_R_WriteConsole = ptr_R_WriteConsole_old;
    R_Consolefile = R_Consolefile_old;

    return 0;
}
