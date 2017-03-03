// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpectrum.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TLine.h"

// C++ headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <string>

//NPTool headers
#include "TTiaraBarrelData.h"
#include "NPCalibrationSource.h"
#include "NPEnergyLoss.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

using namespace::std;

// const
const double striphalflength = 48.4;

//global variables

//functions

// MAIN
void CalibrateBarrelPosition(TString tripleAlphaFileName="./calibrationData/EXPT4/ER193_0.root", TString plotsFileName="./inspectBarrelHisto.root"){ 

/*
*
* In progress
*
*/

}
