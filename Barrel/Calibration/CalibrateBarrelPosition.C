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
#include "NPVDetector.h" // NPL::itoa fucntion
#include "NPCalibrationSource.h"
#include "NPEnergyLoss.h"
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#include "NPCalibrationManager.h"

#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

using namespace::std;

//global variables
int gStripNumber ;
int gDetectorNumber  ;
vector <double> gDataSum1,gDataSumerr1,gDataSum2,gDataSumerr2,gDataSum3,gDataSumerr3,gSlicePos,gZeroVector;
double gPos[2]; // [0] ->downstream strip position; [1] -> Upstream strip position
const double* gFinalCalParam;
const double* gFinalCalParamError;
NPL::EnergyLoss* gELossAlphaInSi ;
TH1F* gDownstreamOffset = new TH1F("DownstreamOffset","DownstreamOffset",200,-50,+50); // keV
TH1F* gUpstreamOffset = new TH1F("UpstreamOffset","UpstreamOffset",200,-50,+50); // keV

//functions
double FitPosition(TH1F* hist, TString pathToFittingBounds, int du, int detector, int strip); // function fitting the positions
double ErfFunction(Double_t *x, Double_t *par); // A special erf function used in the position fit
void ClearGlobalParameters(void); // clear global paramters, where calibration coefficients are stored
double ApplyCalibration(double U, double D); // uses calibration parameters to convert from upstream and downstream channel numbers to upstream and downstream energies
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString pathToMatchsticks, TString outputRootFileName); // uses the specified "alphaCalibrationFile" to make a file with histograms used to extract calibration coefficients
double fUpstream_E(double energy, unsigned short wedge, unsigned short ring); //Changes energy into linearized energy using matchstick data
double fDownstream_E(double energy, unsigned short wedge, unsigned short sector); //Changes energy into linearized energy using matchstick data



// MAIN
void CalibrateBarrelPosition(TString tripleAlphaFileName="/media/sh00319/DellPortableHardDrive/25Mg_dp/CalibrationFiles/NewAlphaData_2017/CombinedNewAlpha.root",
					           TString pathToMatchsticks="./Matchsticks_Calib_220118.txt",
									 	 TString pathToFittingBounds="./FittingBounds.txt"){

  //generate the outputFileName
  TString plotsFileName( tripleAlphaFileName( tripleAlphaFileName.Last('/')+1, tripleAlphaFileName.Length() ) );
  plotsFileName = "./inspectBarrelHistoPosition_250118_"+plotsFileName;

  //global variable
  gELossAlphaInSi = new NPL::EnergyLoss("He4_Si.SRIM","SRIM",100);

  //local variable
  TString CalibfName = "./Barrel_Position_Calib.txt";
  ofstream outputFile;
	outputFile.open(CalibfName.Data());
  TCanvas* canpos[8]; // initialises 8 canvases; 1 for each Barrel detector element, used for position fit

  TFile* fileToCalibrate;
  if(gSystem->AccessPathName(plotsFileName)){ //checks if the file exist already, condition is "true" if not
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
    cout << " No file to calibrate found - creating one now using triple alpha spectra..." << endl;
    fileToCalibrate = CreateFileToCalibrate(tripleAlphaFileName, pathToMatchsticks, plotsFileName);
	fileToCalibrate->Write();
  }
  else {
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
	cout << " The file " << plotsFileName << " is found in the present directory, it will be used for calibration " << endl;
	fileToCalibrate = new TFile(plotsFileName,"UPDATE");
  }
  fileToCalibrate->cd();


  for (int detector=1; detector<=8; detector++){
    gDetectorNumber = detector;
    TString namepos = Form("BarrelPos%d",detector);
    canpos[detector-1] = new TCanvas(namepos,namepos,650,650);
		canpos[detector-1]->Divide(4,2); // splits each canvas into 8; 2 for each strip in a detector

    for (int strip=1; strip<=4; strip++){
      gStripNumber = strip ;
		if ((detector==1 && strip == 3) || (detector==3) ||
			    (detector==5 && strip==1) || (detector==5 && strip==3) ||
			    (detector==6 && strip==2) || (detector==7 && strip==3)) {
				cout << "Detector " << detector << " and Strip " << strip << " is a broken channel. Skipping..." << endl;
				continue;
			}
      ClearGlobalParameters(); // clear parameters
      TString hname;

      // *************** Calibrate Positions *****************************
      for (int du=0; du<2; du++){
			  hname = Form("TIARABARREL_B%d_DOWNSTREAM%d_P",detector,strip);
				if (du == 1)
				  hname = Form("TIARABARREL_B%d_UPSTREAM%d_P",detector,strip);
				TH1F* currentHistPos = (TH1F*) fileToCalibrate->FindObjectAny(hname.Data());
				canpos[detector-1]->cd(2*(strip-1)+(du+1));
        if (currentHistPos && currentHistPos->Integral()>2000) {
					cout << "Position calibration: Working on " << hname << endl;
		      gPos[du] = FitPosition(currentHistPos, pathToFittingBounds, du, detector, strip);
          currentHistPos->Draw(); //Draw for inspection
          //cout << du << " " << gPos[du] << endl ;
				}
				else { // for when there is no histogram currentHistPos or if currentHistPos is empty
				  double kdummy;
          if (du==0) kdummy = -0.71; // default value
          else kdummy = +0.71;
          gPos[du]=kdummy;
        }
      }

      outputFile << detector << " " << strip << " " << "0" << " " << gPos[0] << endl;
      outputFile << detector << " " << strip << " " << "1" << " " << gPos[1] << endl;


    } //strip

  } //detector

  outputFile.close();


  cout << "...done!" << endl;

}
/*****************************************************************************************************************/
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString pathToMatchsticks, TString outputRootFileName){

  //Nptool data
  TTiaraBarrelData* barrelData = new TTiaraBarrelData;

//initiate matchstick calibrator
  CalibrationManager* Cal  = CalibrationManager::getInstance();
  Cal->AddFile(pathToMatchsticks.Data());
  for(int i = 0 ; i < 8 ; ++i){
    for( int j = 0 ; j < 4 ; ++j){
      Cal->AddParameter("TIARABARREL","B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_MATCHSTICK","TIARABARREL_B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_MATCHSTICK")   ;
      Cal->AddParameter("TIARABARREL","B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_MATCHSTICK","TIARABARREL_B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_MATCHSTICK")   ;
    }
  }
  Cal->LoadParameterFromFile();

  //initiate list of Histograms
  TH1F* barrelFrontStripP[8][4][2]; // 8 sides of the barrel, 4 strips in one side, up (1) and down (0)

  TString nameTitle; // same as NPTool calibration token
  for (int iSide =0; iSide<8 ; iSide++) {
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UPSTREAM%d_P",iSide+1,iStrip+1);
	      barrelFrontStripP[iSide][iStrip][1]= new TH1F (nameTitle,nameTitle,400,+0.3,+1);
		  nameTitle =Form("TIARABARREL_B%d_DOWNSTREAM%d_P",iSide+1,iStrip+1);
		  barrelFrontStripP[iSide][iStrip][0]= new TH1F (nameTitle,nameTitle,400,-1,-0.3);
		}
  }

  //Open TFile
  TFile* nptDataFile = new TFile(alphaCalibrationFile.Data(),"READ");
  if (nptDataFile == 0) {
    // if we cannot open the file, print an error message and return immediatly
    printf("Error: cannot open this file: %s \n",alphaCalibrationFile.Data());
    exit(-1);
    return NULL;
  }

  //Load Tree
  TTree* tree = (TTree*) nptDataFile->Get("T40Tree");
  //tree->Print(); // print a summary of the tree contents
  //Set Branch
  tree->SetBranchStatus("TiaraBarrel", true );
  tree->SetBranchStatus("fTiaraBarrel_*", true );
  tree->SetBranchAddress("TiaraBarrel", &barrelData );

	//Loop on tree and fill the histograms
	int entries = tree->GetEntries();
	cout << " INFO: Number of entries in tree: " << entries << endl;
	for(int i = 0 ; i < entries; i++) {
	  if (i%(entries/100)) printf("\r treated %2.f percent ",100.0*i/entries);
	  tree->GetEntry(i);
	  unsigned int sizeUE = barrelData->GetFrontUpstreamEMult();
	  unsigned int sizeDE = barrelData->GetFrontDownstreamEMult();

      //match and fill 2D histograms
	  for(unsigned int iU = 0 ; iU < sizeUE ; ++iU){
		  unsigned short sideU = barrelData->GetFrontUpstreamEDetectorNbr(iU);
		  unsigned short stripU  = barrelData->GetFrontUpstreamEStripNbr(iU);
	      for(unsigned int iD = 0 ; iD < sizeDE ; ++iD){
			  unsigned short sideD = barrelData->GetFrontDownstreamEDetectorNbr(iD);
			  unsigned short stripD  = barrelData->GetFrontDownstreamEStripNbr(iD);
			  if( sideU==sideD && stripU==stripD ){
				  double energyU = fUpstream_E(barrelData->GetFrontUpstreamEEnergy(iU),sideU,stripU);
				  double energyD = fDownstream_E(barrelData->GetFrontDownstreamEEnergy(iD),sideD,stripD);
				  double P = (energyU-energyD)/(energyD+energyU);
				  double E = energyD+energyU;
				  if(energyU>50 && energyD>50 && E>700 && P<-0.3)
				    barrelFrontStripP[sideU-1][stripU-1][0]->Fill(P);
				  if(energyU>50 && energyD>50 && E>700 && P>+0.3)
				    barrelFrontStripP[sideU-1][stripU-1][1]->Fill(P);
			 }
		  }
	  }
	}// end loop on tree
  nptDataFile->Close();

  TFile* output = new TFile(outputRootFileName, "RECREATE");
  output->cd();
  for (int iSide =0; iSide<8 ; iSide++) {
  	for(int iStrip=0 ; iStrip<4 ; iStrip++){
  		if (barrelFrontStripP[iSide][iStrip][0]->GetEntries()>0){
            barrelFrontStripP[iSide][iStrip][0]->Write();
            barrelFrontStripP[iSide][iStrip][1]->Write();
  		}
  	}// end of loop
  }
  output->Close();

  return output;
}



/*****************************************************************************************************************/
void ClearGlobalParameters(void){

	gDataSum1.clear();
	gDataSumerr1.clear();
	gDataSum2.clear();
	gDataSumerr2.clear();
	gDataSum3.clear();
	gDataSumerr3.clear();
	gSlicePos.clear();
	gZeroVector.clear();

	gFinalCalParam=NULL;
	gFinalCalParamError=NULL;
}



///////////////////////////////////////////////////////////////////////////////
double fUpstream_E(double energy, unsigned short side, unsigned short strip){
  static string name; name = "TIARABARREL/B" ;
  name+= NPL::itoa( side ) ;
  name+= "_UPSTREAM" ;
  name+= NPL::itoa( strip ) ;
  name+= "_MATCHSTICK";

  return CalibrationManager::getInstance()->ApplyCalibration(name,
      energy );
}
///////////////////////////////////////////////////////////////////////////////
double fDownstream_E(double energy, unsigned short side, unsigned short strip){
  static string name; name ="TIARABARREL/B" ;
  name+= NPL::itoa( side ) ;
  name+= "_DOWNSTREAM" ;
  name+= NPL::itoa( strip ) ;
  name+= "_MATCHSTICK";
  return CalibrationManager::getInstance()->ApplyCalibration(name,
      energy );
}



double ErfFunction(Double_t *x, Double_t *par) {
  // Par[0]: amplitude measured from the y = 0, amp is positive for "S" shape
  // Par[1]: inflexion point position, (what we are looking for)
  // Par[2]: steepness of the curve, (the smaller the steeper)

   double xx = x[0];
   double f = (par[0]/2)*TMath::Erf((xx-par[1])/par[2]) + TMath::Abs(par[0]/2);
   return f;
}

double FitPosition(TH1F* hist, TString pathToFittingBounds, int du, int detector, int strip){
  double position = -2 ; // default
	ifstream inputFile;
	inputFile.open(pathToFittingBounds);
	while (!inputFile.eof()) {
		int detectorNumber, stripNumber, directionIndicator;
		double lowerBound, upperBound;
		inputFile >> detectorNumber >> stripNumber >> directionIndicator >> lowerBound >> upperBound;
		if (detectorNumber==detector && stripNumber==strip && directionIndicator==du) {
			double rmin = lowerBound;
			double rmax = upperBound;
			//[!] du = {DownStream,UpStream} = {0,1}  => sign={+,-}
  		hist->Rebin();
  		//Set sign for Upstream by default, rmin, rmax are the same
  		double sign = -1 ;
  		// if downstream histo: change sign, swap rmin, rmax
  		if (du==0){
      	sign = +1 ;
   		}
   		// set starting values
   		double amp = sign*1000 ;
   		double pos = -sign*0.70;
   		double steepness = 0.03 ;

    	TF1* fitFunction = new TF1("fitFunction", ErfFunction ,rmin, rmax, 3); // 3 is number of parameters
			fitFunction->SetParameters(amp,pos,steepness);
			fitFunction->SetParNames("Amp","Pos","Steepness");
			fitFunction->SetParLimits(0,0,amp*5);
			//fitFunction->SetParLimits(1,pos-0.05,pos+0.05);
			fitFunction->SetParLimits(2,0.001,steepness*10);
			//fitFunction->FixParameter(2,steepness);

			hist->Fit(fitFunction,"RQM");
    	hist->SetLineColor(8);
    	TF1* fit = hist->GetFunction(fitFunction->GetName());
    	position = fit->GetParameter(1);
      break;
		}// closing if statement
	} // closing while loop
	if(inputFile.eof()) cout << " Detector: " << detector << " strip: " << strip << " (DS<0>/US<1>): " << du << " are not found." << endl;
	inputFile.close();

	return position;
}
