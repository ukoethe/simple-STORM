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
#include <set>

#define R_INTERFACE_PTRS
#define R_NO_REMAP
#include <Rembedded.h>
#include <Rinterface.h>

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

    char *Rargv[] = {"REmbeddedStorm", "--silent", "--no-save"};
    R_SignalHandlers = FALSE;
    Rf_initEmbeddedR(sizeof(Rargv) / sizeof(Rargv[0]), Rargv);

    try
    {

        MultiArray<3,float> in;
        typedef MultiArrayShape<3>::type Shape;

        MyImportInfo info(argc, argv);

        if(info.verbose) {
		   std::cout << "thr:" << info.getThreshold() << " factor:" << info.getFactor() << std::endl;
	    }
        //~ in.reshape(info.shape());
        //~ readVolume(info, in);
        int stacksize = info.shape()[2];
        Size2D size2 (info.shape(0), info.shape(1)); // isnt' there a slicing operator?


        if(info.verbose) {
            std::cout << "Images with Shape: " << info.shape() << std::endl;
            std::cout << "Processing a stack of " << stacksize << " images..." << std::endl;
        }


        // found spots. One Vector over all images in stack
        // the inner set contains all spots in the image
        std::vector<std::set<Coord<float> > > res_coords(stacksize);

        DImage res((size2-Diff2D(1,1))*info.getFactor()+Diff2D(1,1));
        // check if outfile is writable, otherwise throw error -> exit
        exportImage(srcImageRange(res), ImageExportInfo(info.getOutfile().c_str()));
        if(!info.getCoordsfile().empty()) {
            std::ofstream cf (info.getCoordsfile().c_str());
            vigra_precondition(cf.is_open(), "Could not open coordinate-file for writing.");
            cf.close();
        }

        USETICTOC;
        TIC;  // measure the time

        // STORM Algorithmut
        //printIntensities(info);

        float a,offset, intercept;
        std::vector<float> parameterTrafo(6);
        findCorrectionCoefficients(info, parameterTrafo); 	//first 3 entries are parameters for the raw signal to poisson transformation
        													//the last 3 entries are parameters for the poission to gaussian with sigma = 1 transformation
        //parameterTrafo[0] = 1.5;
		//parameterTrafo[1] = 40;
//		parameterTrafo[2] = 0;
        std::cout<<"a: "<<parameterTrafo[0]<<" b: "<<parameterTrafo[1]<<" intercept: "<<parameterTrafo[2]<<std::endl;
        //showPoisson(info, parameterTrafo);

        //parameterTrafo[4] = 0;
        //parameterTrafo[5] = 3./8.;

        int w= info.shape(0), h =info.shape(1);
        int vecw[] = {16,30,40,50,60,16,26};
        int vech[] = {39,10,10,10,10,40,20};
        printIntensities(info, vecw, vech, 7, parameterTrafo[0],parameterTrafo[1]);

        wienerStorm(info, res_coords,parameterTrafo);

        // resulting image
        drawCoordsToImage<Coord<float> >(res_coords, res);

        int numSpots = 0;
        if(!info.getCoordsfile().empty()) {
            numSpots = saveCoordsFile(info.getCoordsfile(), res_coords, info.shape(), info.getFactor(), info.getPixelsize());
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
        exportImage(srcImageRange(res), ImageExportInfo(info.getOutfile().c_str()));



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
