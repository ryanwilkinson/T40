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

// const
const double striphalflength = 48.4;

//global variables
vector <double> gDataSum1,gDataSumerr1,gDataSum2,gDataSumerr2,gDataSum3,gDataSumerr3,gPos,gZeroVector;
const double* gFinalCalParam;
const double* gFinalCalParamError;
NPL::EnergyLoss* gELossAlphaInSi ;

//functions
TH2D* FindHistogram(TString histname, TString filename); // finds a histogram "histname" in file "filename"
void SliceHistogram(TH2D* hist, int binstep,float ylow, float yup, float xmin, float xmax); // slices histogram "hist" between ylow and yup in steps of "binstep" bins, and select the x-range
vector<double> FitOneSlice(TH1D* hist); // fits a single slice of histogram "hist"
int Minimise(void); // numerical minimisation function which produces the final calibration parameters
double GetChiSquared(const double parameters[]); // returns a chi squared value to quantify the "goodness of fit" of a set of paramters "paramters[]" relative to the data
vector<double> CalculateEnergySum(const double parameters[], vector<double> p, double energy); // calculates total energy from EU, ED and BD
void ShowControl2DSpectra(TH2D* h2, TCanvas* c, int num); //  show spectra used to calculate calibration coefficients - useful for checking results
void ClearGlobalParameters(void); // clear global paramters, where calibration coefficients are stored
double PosToAngle(double pos); // converts a value of POS to an angle theta, used when calculating energy losses
double ApplyCalibration(double U, double D); // uses calibration parameters to convert from upstream and downstream channel numbers to upstream and downstream energies
TCanvas* ShowCalibration(int det, int strip); // plots E vs POS for a given detector and strip
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString outputRootFileName); // uses the specified "alphaCalibrationFile" to make a file with histograms used to extract calibration coefficients
//functions
double fUpstream_E(double energy, unsigned short wedge, unsigned short ring); //Changes energy into linearized energy using matchstick data
double fDownstream_E(double energy, unsigned short wedge, unsigned short sector); //Changes energy into linearized energy using matchstick data



// MAIN
void CalibrateBarrel(TString tripleAlphaFileName="../../../Tiara/TapeData/Root/POST/ER1_1.root", TString plotsFileName="./inspectBarrelHisto.root"){ // tripleAlphaFileName = run file with triple alpha spectra for the Barrel in it

  //global variable
  gELossAlphaInSi = new NPL::EnergyLoss("He4_Si.SRIM","SRIM",100);

  //local variable
  ofstream outputFile;
  outputFile.open("Barrel_Calib.txt");
  TCanvas* can[8]; // initialises 8 canvases; 1 for each Barrel detector element

  TFile* fileToCalibrate = new TFile(plotsFileName);
  if (fileToCalibrate->IsZombie()){
    cout << "No file to calibrate found - creating one now using triple alpha spectra..." << endl;
    fileToCalibrate = CreateFileToCalibrate(tripleAlphaFileName, plotsFileName);
  }
  //fileToCalibrate->Open(plotsFileName);
  //TString filename = "inspectBarrelHisto.root"; // the name of the root file made in the line above

  for (int detector=1; detector<=8; detector++){
    TString name = Form("Barrel%d",detector);
    can[detector-1] = new TCanvas(name,name,650,650);
    can[detector-1]->Divide(2,2); // splits each canvas into 4; one for each strip in a detector
    for (int strip=1; strip<=4; strip++){
      ClearGlobalParameters(); // clear parameters
	    TString hname = Form("TIARABARREL_B%d_PE%d_E",detector,strip); // histograms of (Upstream-Downstream)/(Upstream+Downstream) vs Upstream+Downstream
	    TH2D* h2 = FindHistogram(hname, plotsFileName);
      can[detector-1]->cd(strip);
      if (h2 && h2->Integral()) {
        cout << " Working on " << hname << endl;
			  SliceHistogram(h2,10,-0.67,0.67,600,1300); // slice between -0.67 and 0.67 with a step of 10 bins
			  int result = Minimise(); // performs the numerical minimisation and saves the final calibration values into global variable gFinalCalParam and gFinalCalParamError
        ShowControl2DSpectra(h2,can[detector-1],strip);
        double BDtune=1; // this variable allows one to change the degree of the ballistic deficit (BD) if necessary
			  outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_E " << gFinalCalParam[1] << " " << gFinalCalParam[0] << endl; // calibration coefficients and nptool tokens
			  outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_E " << gFinalCalParam[3] << " " << gFinalCalParam[2] << endl;
			  if(-gFinalCalParam[4]<0) { // BD should be positive
          cout << " WARNING: Ballistic deficit is negative, replacing by value zero "<<endl;
				  BDtune=0;
			  }
        outputFile << "TIARABARREL_BALLISTIC_B" << detector << "_STRIP" << strip << " " << 0 << " " << 0 << " " << -gFinalCalParam[4]*BDtune << endl;
			  outputFile << "TIARABARREL_B" << detector << "_STRIP" << strip << "_POS " << 0 << " " << 0.66 << endl;
        ShowCalibration(detector, strip)->Draw(); // on separate graph
      }
      else { // for when there is no histogram h2 or if h2 is empty - nptool tokens with default calibration parameters
		    outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_E 0 0 " << endl;
		    outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_E 0 0 " << endl;
	      outputFile << "TIARABARREL_BALLISTIC_B" << detector << "_STRIP" << strip << " 0 0 0 " << endl;
		    outputFile << "TIARABARREL_B" << detector << "_STRIP" << strip << "_POS " << 0 << " " << 0.66 << endl;
      }
    }
  }
  outputFile.close();
  cout << "...done!" << endl;

}
/*****************************************************************************************************************/
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString outputRootFileName){

  //Nptool data
  TTiaraBarrelData* barrelData = new TTiaraBarrelData;

//initiate matchstick calibrator
  CalibrationManager* Cal  = CalibrationManager::getInstance();
    Cal->AddFile("/home/mhd/Work/Tamu/T40/Matchsticks/Files/Matchsticks_Calib.txt");

  for(int i = 0 ; i < 8 ; ++i){
    for( int j = 0 ; j < 4 ; ++j){
      Cal->AddParameter("TIARABARREL","MATCHSTICK_B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_E","TIARABARREL_MATCHSTICK_B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_E")   ;
      Cal->AddParameter("TIARABARREL","MATCHSTICK_B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_E","TIARABARREL_MATCHSTICK_B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_E")   ;
    cout << "TIARABARREL_"<<"MATCHSTICK_B" << NPL::itoa(i+1) << "_UPSTREAM" <<NPL::itoa(j+1) <<"_E" << std::endl;
    }
  }
    Cal->LoadParameterFromFile();

  //initiate list of Histograms
  TH1F* barrelFrontStripP[8][4][2]; // 8 sides of the barrel, 4 strips in one side, up and down
  TH2F* barrelFrontStripDU[8][4];  // Downstream vs Upstream
  TH2F* barrelFrontStripPE[8][4]; // Upstream-Downstream/(sum) vs Upstream+Downstream
  TString nameTitle; // same as NPTool calibration token
  for (int iSide =0; iSide<8 ; iSide++) {
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UPSTREAM%d_PU",iSide+1,iStrip+1);
	    barrelFrontStripP[iSide][iStrip][0]= new TH1F (nameTitle,nameTitle,1500,50,1550);
		  nameTitle =Form("TIARABARREL_B%d_DOWNSTREAM%d_PD",iSide+1,iStrip+1);
		  barrelFrontStripP[iSide][iStrip][1]= new TH1F (nameTitle,nameTitle,1500,50,1550);
		}
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UD%d_E",iSide+1,iStrip+1);
		  barrelFrontStripDU[iSide][iStrip]= new TH2F (nameTitle,nameTitle,1500,50,1550,500,50,1550);
		  nameTitle =Form("TIARABARREL_B%d_PE%d_E",iSide+1,iStrip+1); // the ones we're interested in making
		  barrelFrontStripPE[iSide][iStrip]= new TH2F (nameTitle,nameTitle,1500,50,1550,500,-1,+1);
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
				  barrelFrontStripDU[sideU-1][stripU-1]->Fill(energyU,energyD);
				  barrelFrontStripPE[sideU-1][stripU-1]->Fill(energyU+energyD,(energyU-energyD)/(energyD+energyU));
				  if(energyD>0 && (energyD/(energyU+energyD)>0.90)) {barrelFrontStripP[sideU-1][stripU-1][1]->Fill((energyU-energyD)/(energyU+energyD));}
				}
			}
	  }
	}// end loop on tree
  nptDataFile->Close();

  TFile* output = new TFile(outputRootFileName, "RECREATE");
  output->cd();
  for (int iSide =0; iSide<8 ; iSide++) {
  	for(int iStrip=0 ; iStrip<4 ; iStrip++){
  		nameTitle = barrelFrontStripDU[iSide][iStrip]->GetTitle();
  		if (barrelFrontStripDU[iSide][iStrip]->GetEntries()>0){
  			barrelFrontStripPE[iSide][iStrip]->Write();
        barrelFrontStripDU[iSide][iStrip]->Write();
  		}
  	}// end of loop
  }
  output->Close();

  return output;
}

/*****************************************************************************************************************/
void SliceHistogram(TH2D* h2, int binstep, float ylow, float yup, float xmin, float xmax){

	TAxis* Yaxis = h2->GetYaxis();
	h2->GetXaxis()->SetRangeUser(xmin,xmax);
	int firstbin = Yaxis->FindBin(ylow); // determining upper and lower y bins
	int lastbin = Yaxis->FindBin(yup);
	int bin = firstbin;
	while(bin<lastbin){ // total of 31 bins over entire range considered currently
		float y1 = Yaxis->GetBinCenter(bin);
		float y2 = Yaxis->GetBinCenter(bin+binstep);
		TH1D* h1 = h2->ProjectionX(Form("range_%2.2f_%2.2f",y1,y2),bin,bin+binstep);
    h1->Rebin(2);
		vector <double> parameters = FitOneSlice(h1); // only works for slices with 3 peaks and more than 500 entries per slice
		//cout << " size " << parameters.size()<< endl;
		if (parameters.size()==6){
			gDataSum1.push_back(parameters[0]);
			gDataSumerr1.push_back(parameters[1]);
			gDataSum2.push_back(parameters[2]);
			gDataSumerr2.push_back(parameters[3]);
			gDataSum3.push_back(parameters[4]);
			gDataSumerr3.push_back(parameters[5]);
			gPos.push_back(0.5*(y1+y2));
		  gZeroVector.push_back(0);
		}
		bin+=binstep+1;
	}
}
/*****************************************************************************************************************/
vector <double> FitOneSlice(TH1D* Slice){

  //Slice->Rebin();
	vector <double> param;
	vector<double> PeakPositions;
	TSpectrum* Spec = new TSpectrum();
	int NumPeaksFound = Spec->Search(Slice, 3, "nodraw", 0.2);
	double* peaks = Spec->GetPositionX();

	for (int j=0 ; j<NumPeaksFound ; j++){
		//cout << "Peak at = " << peaks[j] << endl;
		PeakPositions.push_back(peaks[j]);
	}
	sort(PeakPositions.begin(), PeakPositions.end());
  if (NumPeaksFound==3 && Slice->Integral()>400 ){
		TF1* fittingfunc = new TF1("fittingfunc", "gaus(0)+gaus(3)+gaus(6)", PeakPositions.front()-20, PeakPositions.back()+20);
		fittingfunc->SetParameter(1,PeakPositions[0]);
		fittingfunc->SetParameter(4,PeakPositions[1]);
		fittingfunc->SetParameter(7,PeakPositions[2]);
		fittingfunc->SetParameter(2,7.5);
		fittingfunc->SetParameter(5,7.5);
		fittingfunc->SetParameter(8,7.5);
		Slice->Fit("fittingfunc", "nodrawQR");
		param.push_back(fittingfunc->GetParameter(1));
		param.push_back(fittingfunc->GetParameter(2));
		param.push_back(fittingfunc->GetParameter(4));
		param.push_back(fittingfunc->GetParameter(5));
		param.push_back(fittingfunc->GetParameter(7));
		param.push_back(fittingfunc->GetParameter(8));
		delete fittingfunc;
  }
 return param;
}
/*****************************************************************************************************************/
TH2D* FindHistogram(TString HistoName, TString filename){

  TFile* file = new TFile(filename.Data());
  TH2D* h2 = (TH2D*) file->FindObjectAny(HistoName.Data());

  return ((TH2D*) h2);
}
/*****************************************************************************************************************/
int Minimise(void){

  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Fumili");

	min->SetMaxFunctionCalls(1000000000);
	min->SetMaxIterations(1000000000);
	min->SetTolerance(0.001);

	ROOT::Math::Functor f(&GetChiSquared,5);
	double step[5] = {0.01,0.1,0.01,0.1,0.0001};
	double variable[5] = {6.5, -300, 6.5, -300, 0.01}; // change the other set
	// parameters: upstream gain, upstream offset, downstream gain, downstream offset, quadratic BD term, linear BD term, constant BD term
	min->SetFunction(f);
	// Set the free variables to be minimized
  for(unsigned int i = 0 ; i < 5 ; i++)
	min->SetVariable(i,Form("Par%i",i),variable[i], step[i]);

	//min->SetVariableLimits(0,5,7);
	//min->SetVariableLimits(1,-100,0);
	//min->SetVariableLimits(2,5,7);
	//min->SetVariableLimits(3,-100,0);
	//min->FixVariable(4);
	//min->FixVariable(0);
	//min->FixVariable(2);
	//min->FixVariable(1);
	//min->FixVariable(3);
	min->Minimize();

  //gFinalCalParam = variable ;
  gFinalCalParam = min->X();
  gFinalCalParamError = min->Errors();
  return 0;
}
/*****************************************************************************************************************/
double GetChiSquared(const double parameters[]){

	vector<double> calSumc1 = CalculateEnergySum(parameters, gPos, 5156.59);
	vector<double> calSumc2 = CalculateEnergySum(parameters, gPos, 5485.56);
	vector<double> calSumc3 = CalculateEnergySum(parameters, gPos, 5804.77);

	vector<double> bitsOfChi1, bitsOfChi2, bitsOfChi3;
  double ChiSquared=0;
  double diff = 0 ;
  for(unsigned int i=0; i<gPos.size(); i++){
    diff = (calSumc1[i]-gDataSum1[i]);
    bitsOfChi1.push_back(diff*diff/(gDataSumerr1[i]*gDataSumerr1[i]));
  }
  for(unsigned int i=0; i<gPos.size(); i++){
    diff = (calSumc2[i]-gDataSum2[i]);
    bitsOfChi2.push_back(diff*diff/(gDataSumerr2[i]*gDataSumerr2[i]));
  }
  for(unsigned int i=0; i<gPos.size(); i++){
    diff = (calSumc3[i]-gDataSum3[i]);
    bitsOfChi3.push_back(diff*diff/(gDataSumerr3[i]*gDataSumerr3[i]));
  }

  for(unsigned int i=0; i<bitsOfChi1.size(); i++){
    ChiSquared = ChiSquared+bitsOfChi1[i];
  }
  for(unsigned int i=0; i<bitsOfChi2.size(); i++){
    ChiSquared = ChiSquared+bitsOfChi2[i];
  }
  for(unsigned int i=0; i<bitsOfChi3.size(); i++){
    ChiSquared = ChiSquared+bitsOfChi3[i];
  }
//  cout << "Final value of Chi2 is " << ChiSquared << endl; // used for testing
  return ChiSquared;
}
/*****************************************************************************************************************/
vector<double> CalculateEnergySum(const double parameters[], vector<double> p, double energy){

	vector<double> calSum;
	double k = 0.66;

  for (unsigned int i = 0 ; i < p.size() ; i++){
	  // pass from x and y to SUM=x+y and POS=(y-x)
	  // where y: U = 0.5*SUM*(1+POS)   x: D = 0.5*SUM*(1-POS)
    double angle = PosToAngle(p[i]);
    double slow = gELossAlphaInSi->Slow(energy*keV,1*micrometer,angle)/keV;
	  double s = slow/(1-parameters[4]*(k*k-p[i]*p[i])) - parameters[1] - parameters[3]; // the ballistic deficit is subtracted here, this should lead to BD>0
	  s = 2*s/(p[i]*(parameters[0]-parameters[2]) + parameters[0] + parameters[2] );
	  calSum.push_back(s);
  }
  return calSum;
}
/*****************************************************************************************************************/
void ShowControl2DSpectra(TH2D* h2, TCanvas* can, int num){
	/*
  cout << "Minimised calibration parameters" << endl;
	cout << "Param \t Error " << endl;
	cout << gFinalCalParam[0] << " \t+/- " <<gFinalCalParamError[0] <<  endl;
	cout << gFinalCalParam[0] << " \t+/- " <<gFinalCalParamError[0] <<  endl;
	cout << gFinalCalParam[1] << " \t+/- " <<gFinalCalParamError[1] <<  endl;
	cout << gFinalCalParam[2] << " \t+/- " <<gFinalCalParamError[2] <<  endl;
	cout << gFinalCalParam[3] << " \t+/- " <<gFinalCalParamError[3] <<  endl;
	cout << gFinalCalParam[4] << " \t+/- " <<gFinalCalParamError[4] <<  endl;
    */

  can->cd(num);
  h2->Draw("colz");
  TGraphErrors* DataAlphaSet1 = new TGraphErrors(gDataSum1.size(), &gDataSum1[0], &gPos[0], &gDataSumerr1[0], &gZeroVector[0]);
  TGraphErrors* DataAlphaSet2 = new TGraphErrors(gDataSum2.size(), &gDataSum2[0], &gPos[0], &gDataSumerr2[0], &gZeroVector[0]);
  TGraphErrors* DataAlphaSet3 = new TGraphErrors(gDataSum3.size(), &gDataSum3[0], &gPos[0], &gDataSumerr3[0], &gZeroVector[0]);

	DataAlphaSet1->Draw("P same");
	DataAlphaSet1->SetMarkerStyle(20);
	DataAlphaSet1->SetMarkerColor(kRed);
	DataAlphaSet2->Draw("P same");
	DataAlphaSet2->SetMarkerStyle(20);
	DataAlphaSet2->SetMarkerColor(kBlue);
	DataAlphaSet3->Draw("P same");
	DataAlphaSet3->SetMarkerStyle(20);
	DataAlphaSet3->SetMarkerColor(kGreen);

	vector<double> calSum1=CalculateEnergySum(gFinalCalParam, gPos, 5156.59);
	vector<double> calSum2=CalculateEnergySum(gFinalCalParam, gPos, 5485.56);
	vector<double> calSum3=CalculateEnergySum(gFinalCalParam, gPos, 5804.77);

	TGraphErrors* calAlphaSet1 = new TGraphErrors(calSum1.size(), &calSum1[0], &gPos[0], &gZeroVector[0], &gZeroVector[0]);
	TGraphErrors* calAlphaSet2 = new TGraphErrors(calSum2.size(), &calSum2[0], &gPos[0], &gZeroVector[0], &gZeroVector[0]);
	TGraphErrors* calAlphaSet3 = new TGraphErrors(calSum3.size(), &calSum3[0], &gPos[0], &gZeroVector[0], &gZeroVector[0]);

	calAlphaSet1->Draw("P same");
	calAlphaSet1->SetMarkerStyle(20);
	calAlphaSet1->SetMarkerSize(0.7);
	calAlphaSet1->SetMarkerColor(kYellow);
	calAlphaSet2->Draw("P same");
	calAlphaSet2->SetMarkerStyle(20);
	calAlphaSet2->SetMarkerSize(0.7);
	calAlphaSet2->SetMarkerColor(kYellow);
	calAlphaSet3->Draw("P same");
	calAlphaSet3->SetMarkerStyle(20);
	calAlphaSet3->SetMarkerSize(0.7);
	calAlphaSet3->SetMarkerColor(kYellow);
}
/*****************************************************************************************************************/
void ClearGlobalParameters(void){

	gDataSum1.clear();
	gDataSumerr1.clear();
	gDataSum2.clear();
	gDataSumerr2.clear();
	gDataSum3.clear();
	gDataSumerr3.clear();
	gPos.clear();
	gZeroVector.clear();

	gFinalCalParam=NULL;
	gFinalCalParamError=NULL;
}
/*****************************************************************************************************************/
double PosToAngle(double pos){
	double x = (pos/0.66)*striphalflength; // in mm
	return TMath::ATan(x/33); // 33 mm is the distance from beam spot  (supposed at the center to the strip at 90 degree)
}
/*****************************************************************************************************************/
double ApplyCalibration(double Uch, double Dch){
	double k = 0.66;
	double pos = (Uch-Dch)/(Dch+Uch);
	double U = Uch*gFinalCalParam[0]+ gFinalCalParam[1];
	double D = Dch*gFinalCalParam[2]+ gFinalCalParam[3];
    //return (U+D);
    //return -gFinalCalParam[4](k*k-pos*pos);
	return (U+D)*(1-gFinalCalParam[4]*(k*k-pos*pos));  // the BD<0 thus the negative sugn adds it here!
}
/*****************************************************************************************************************/
TCanvas* ShowCalibration(int det, int strip){

  vector <double> calEnergy,pos;

  for (unsigned int i = 0 ; i < gDataSum1.size(); i++){
     //reconstruct the U (yaxis) vs D(xaxis) plane
     double uch = 0.5 * gDataSum1[i] * (1+gPos[i]); // these definitions are strict
     double dch = 0.5 * gDataSum1[i] * (1-gPos[i]);
     calEnergy.push_back(ApplyCalibration(uch, dch));
     pos.push_back(gPos[i]);
  }

  for (unsigned int i = 0 ; i < gDataSum2.size(); i++){
    double uch = 0.5 * gDataSum2[i] * (1+gPos[i]);
    double dch = 0.5 * gDataSum2[i] * (1-gPos[i]);
    calEnergy.push_back(ApplyCalibration(uch, dch));
    pos.push_back(gPos[i]);
  }

  for (unsigned int i = 0 ; i < gDataSum3.size(); i++){
    double uch = 0.5 * gDataSum3[i] * (1+gPos[i]);
    double dch = 0.5 * gDataSum3[i] * (1-gPos[i]);
    calEnergy.push_back(ApplyCalibration(uch, dch));
    pos.push_back(gPos[i]);
	}

  cout <<  "  size of points vector " << pos.size() << endl;

	TGraph* allAlpha = new TGraph(pos.size(), &calEnergy[0], &pos[0]);
	allAlpha->SetMarkerStyle(20);
	allAlpha->SetMarkerColor(kBlack);
  allAlpha->GetXaxis()->SetTitle("Energy (keV)");
  allAlpha->GetXaxis()->SetRangeUser(4500,6000);
  allAlpha->GetYaxis()->SetTitle("Pos (U-D)/(U+D)");
  allAlpha->GetYaxis()->SetRangeUser(-1,+1);
  allAlpha->SetTitle(Form("CalibratedSpectra_Detector%d_Strip%d",det,strip));

	TCanvas* can = new TCanvas(Form("Cal_det%d_strip%d",det,strip),Form("Cal_det%d_strip%d",det,strip),650,650);
  allAlpha->Draw("ap");
	TLine* lineE1 = new TLine (5156.6,-0.66,5156.6,0.66);
	lineE1->SetLineColor(kRed);
	lineE1->Draw();
	TLine* lineE2 = new TLine (5485.6,-0.66,5485.6,0.66);
	lineE2->Draw();
  lineE2->SetLineColor(kRed);
	TLine* lineE3 = new TLine (5804.8,-0.66,5804.8,0.66);
	lineE3->Draw();
  lineE3->SetLineColor(kRed);

  return can;
}


///////////////////////////////////////////////////////////////////////////////
double fUpstream_E(double energy, unsigned short side, unsigned short strip){
  static string name; name = "TIARABARREL/MATCHSTICK_B" ;
  name+= NPL::itoa( side ) ;
  name+= "_UPSTREAM" ;
  name+= NPL::itoa( strip ) ;
  name+= "_E";

  return CalibrationManager::getInstance()->ApplyCalibration(name,
      energy );
}
///////////////////////////////////////////////////////////////////////////////
double fDownstream_E(double energy, unsigned short side, unsigned short strip){
  static string name; name ="TIARABARREL/MATCHSTICK_B" ;
  name+= NPL::itoa( side ) ;
  name+= "_DOWNSTREAM" ;
  name+= NPL::itoa( strip ) ;
  name+= "_E"; 
  return CalibrationManager::getInstance()->ApplyCalibration(name,
      energy );
}
