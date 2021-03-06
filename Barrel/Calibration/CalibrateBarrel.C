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
TH2D* FindHistogram(TString histname, TString filename); // finds a histogram "histname" in file "filename"
double FitPosition(TH1F* hist, TString pathToFittingBounds, int du, int detector, int strip); // function fitting the positions
void SliceHistogram(TH2D* hist, int binstep,float ylow, float yup, float xmin, float xmax); // slices histogram "hist" between ylow and yup in steps of "binstep" bins, and select the x-range
vector<double> FitOneSlice(TH1D* hist); // fits a single slice of histogram "hist"
double ErfFunction(Double_t *x, Double_t *par); // A special erf function used in the position fit
int Minimise(void); // numerical minimisation function which produces the final calibration parameters
double GetChiSquared(const double parameters[]); // returns a chi squared value to quantify the "goodness of fit" of a set of paramters "paramters[]" relative to the data
vector<double> CalculateEnergySum(const double parameters[], vector<double> p, double energy); // calculates total energy from EU, ED and BD
void ShowControl2DSpectra(TH2D* h2, TCanvas* c, int num); //  show spectra used to calculate calibration coefficients - useful for checking results
void ClearGlobalParameters(void); // clear global paramters, where calibration coefficients are stored
double PosToAngle(double pos); // converts a value of POS to an angle theta, used when calculating energy losses
double ApplyCalibration(double U, double D); // uses calibration parameters to convert from upstream and downstream channel numbers to upstream and downstream energies
TCanvas* ShowCalibration(int det, int strip); // plots E vs POS for a given detector and strip
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString pathToMatchsticks, TString outputRootFileName); // uses the specified "alphaCalibrationFile" to make a file with histograms used to extract calibration coefficients
double fUpstream_E(double energy, unsigned short wedge, unsigned short ring); //Changes energy into linearized energy using matchstick data
double fDownstream_E(double energy, unsigned short wedge, unsigned short sector); //Changes energy into linearized energy using matchstick data



// MAIN
//void CalibrateBarrel(TString tripleAlphaFileName="/home/shuyaota/midas2nptool/root/EXPT5/R1_0.root",
void CalibrateBarrel(TString tripleAlphaFileName="/home/shuyaota/midas2nptool/root/EXPT5/ER395_0.root",
//void CalibrateBarrel(TString tripleAlphaFileName="/home/shuyaota/midas2nptool/root/EXPT5_New170725/R1_0.root",
//void CalibrateBarrel(TString tripleAlphaFileName="/media/gchristian/HD1/T4T/root/EXPT6/ER29_1.root",
					           TString pathToMatchsticks="/home/shuyaota/nptool/Projects/T40/Calibration/Matchsticks_Calib.txt",
									 	 TString pathToFittingBounds="./FittingBounds.txt"){
									 	 //TString pathToFittingBounds="/home/rw00227/nptool/Projects/T40/calibration/FittingBoundsER230.txt"){

  //generate the outputFileName
  TString plotsFileName="./ER29_1_inspectBarrelHisto.root";
  //TString plotsFileName="./calibration/ER230_0_inspectBarrelHisto.root";

  //global variable
  gELossAlphaInSi = new NPL::EnergyLoss("He4_Si.SRIM","SRIM",100);

  //local variable
  TString CalibfName( tripleAlphaFileName( tripleAlphaFileName.Last('/')+1, tripleAlphaFileName.Length() ) );
  CalibfName.ReplaceAll("root","txt");
  CalibfName = "./calibration/Barrel_Calib_"+CalibfName;
  ofstream outputFile;
	outputFile.open(CalibfName.Data());
  TCanvas* can[8]; // initialises 8 canvases; 1 for each Barrel detector element
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
    TString name = Form("Barrel%d",detector);
    can[detector-1] = new TCanvas(name,name,650,650);
    can[detector-1]->Divide(2,2); // splits each canvas into 4; one for each strip in a detector
    TString namepos = Form("BarrelPos%d",detector);
    canpos[detector-1] = new TCanvas(namepos,namepos,650,650);
		canpos[detector-1]->Divide(4,2); // splits each canvas into 8; 2 for each strip in a detector

    for (int strip=1; strip<=4; strip++){
      gStripNumber = strip ;
			if ((detector==1 && strip==3) || (detector==3) || (detector==5 && strip==1) || (detector==5 && strip==3) || (detector==6 && strip==2) || (detector==7 && strip==3)) {
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
      //store in file
      double k = (gPos[1] - gPos[0])/2;
      double d = (gPos[1] + gPos[0])/2;
      outputFile << "TIARABARREL_B" << detector << "_STRIP" << strip << "_POS " << d << " " << k << endl;

      // *************** Calibrate Energies *****************************
	    hname = Form("TIARABARREL_B%d_PE%d_E",detector,strip); // histograms of (Upstream-Downstream)/(Upstream+Downstream) vs Upstream+Downstream
	    TH2D* h2 = FindHistogram(hname, plotsFileName);
      can[detector-1]->cd(strip);
      int result = 0 ;
      if (h2 && h2->Integral()) {
        cout << "Energy calibration: Working on " << hname << endl;
				cout << "gPos[0] = " << gPos[0] << " & gPos[1] = " << gPos[1] << endl;
	      SliceHistogram(h2,10,gPos[0],gPos[1],800,1800); // slice between -0.67 and 0.67 with a step of 10 bins
		    result = Minimise(); // performs the numerical minimisation and saves the final calibration values into global variable gFinalCalParam and gFinalCalParamError
		    }
		  if(result==1){
        ShowControl2DSpectra(h2,can[detector-1],strip);
        double BDtune=1; // this variable allows one to change the degree of the ballistic deficit (BD) if necessary
	      outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_E " << gFinalCalParam[1] << " " << gFinalCalParam[0] << endl; // calibration coefficients and nptool tokens
		    outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_E " << gFinalCalParam[3] << " " << gFinalCalParam[2] << endl;
		    if(-gFinalCalParam[4]<0) { // BD = (-gFinalCalParam[4]) should be positive
          cout << " WARNING: Ballistic deficit is negative, replacing by value zero "<<endl;
		      BDtune=0;
		    }
        outputFile << "TIARABARREL_B" << detector << "_STRIP" << strip << "_BALLISTIC "<<0<< " " <<0<< " " << -gFinalCalParam[4]*BDtune << endl;
        //ShowCalibration(detector, strip)->Draw(); // on separate graph
      }
      else { // for when there is no histogram h2 or if h2 is empty - nptool tokens with default calibration parameters
	      outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_E 0 0 " << endl;
	      outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_E 0 0 " << endl;
        outputFile << "TIARABARREL_B" << detector << "_STRIP" << strip << "_BALLISTIC 0 0 0 " << endl;
      }

    } //strip

  canpos[detector-1]->Draw();
  } //detector

  outputFile.close();

  //offset check; useful for dead zone determination
  TCanvas* canoffset= new TCanvas("offset","offset",650,650);
  canoffset->Divide(1,2);
  canoffset->cd(1);
  gDownstreamOffset->Draw();
  canoffset->cd(2);
  gUpstreamOffset->Draw();
  canoffset->Draw();

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
  TH2F* barrelFrontStripUD[8][4];  // Downstream vs Upstream
  TH2F* barrelFrontStripPE[8][4]; // Upstream-Downstream/(sum) vs Upstream+Downstream
  TString nameTitle; // same as NPTool calibration token
  for (int iSide =0; iSide<8 ; iSide++) {
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UPSTREAM%d_P",iSide+1,iStrip+1);
	      barrelFrontStripP[iSide][iStrip][1]= new TH1F (nameTitle,nameTitle,400,+0.4,+1);
		  nameTitle =Form("TIARABARREL_B%d_DOWNSTREAM%d_P",iSide+1,iStrip+1);
		  barrelFrontStripP[iSide][iStrip][0]= new TH1F (nameTitle,nameTitle,400,-1,-0.4);
		}
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UD%d_E",iSide+1,iStrip+1);
		  barrelFrontStripUD[iSide][iStrip]= new TH2F (nameTitle,nameTitle,2100,-50,2050,2100,-50,2050);
		  nameTitle =Form("TIARABARREL_B%d_PE%d_E",iSide+1,iStrip+1); // the ones we're interested in making
		  barrelFrontStripPE[iSide][iStrip]= new TH2F (nameTitle,nameTitle,2100,-50,2050,500,-1,+1);
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
				  barrelFrontStripUD[sideU-1][stripU-1]->Fill(energyD,energyU);
				  barrelFrontStripPE[sideU-1][stripU-1]->Fill(E,P);
				  if(energyU>50 && energyD>50 && E>700 && P<-0.4)
				    barrelFrontStripP[sideU-1][stripU-1][0]->Fill(P);
				  if(energyU>50 && energyD>50 && E>700 && P>+0.4)
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
  		nameTitle = barrelFrontStripUD[iSide][iStrip]->GetTitle();
  		if (barrelFrontStripUD[iSide][iStrip]->GetEntries()>0){
  			barrelFrontStripPE[iSide][iStrip]->Write();
            barrelFrontStripUD[iSide][iStrip]->Write();
            barrelFrontStripP[iSide][iStrip][0]->Write();
            barrelFrontStripP[iSide][iStrip][1]->Write();
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
			gSlicePos.push_back(0.5*(y1+y2));
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
  if (NumPeaksFound==3 && Slice->Integral()>175 ){
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
	double variable[5] = {6.5, 0, 6.5, 0, 0.05}; // change the other set
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

	double minfail[5] = {0, -100, 0, -100, 0}; // default
	gFinalCalParam= minfail; // change the other set;
  gFinalCalParamError = minfail;

  int result = min->Minimize();
  //cout << " result  " << result << endl;
  //cin.get();
  if (result){
    gFinalCalParam = min->X();
    gFinalCalParamError = min->Errors();

    gDownstreamOffset->Fill(gFinalCalParam[3]);
    gUpstreamOffset->Fill(gFinalCalParam[1]);
  }

  return result;
}
/*****************************************************************************************************************/
double GetChiSquared(const double parameters[]){

	vector<double> calSumc1 = CalculateEnergySum(parameters, gSlicePos, 5156.59);
	vector<double> calSumc2 = CalculateEnergySum(parameters, gSlicePos, 5485.56);
	vector<double> calSumc3 = CalculateEnergySum(parameters, gSlicePos, 5804.77);

	vector<double> bitsOfChi1, bitsOfChi2, bitsOfChi3;
  double ChiSquared=0;
  double diff = 0 ;
  for(unsigned int i=0; i<gSlicePos.size(); i++){
    diff = (calSumc1[i]-gDataSum1[i]);
    bitsOfChi1.push_back(diff*diff/(gDataSumerr1[i]*gDataSumerr1[i]));
  }
  for(unsigned int i=0; i<gSlicePos.size(); i++){
    diff = (calSumc2[i]-gDataSum2[i]);
    bitsOfChi2.push_back(diff*diff/(gDataSumerr2[i]*gDataSumerr2[i]));
  }
  for(unsigned int i=0; i<gSlicePos.size(); i++){
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

  for (unsigned int i = 0 ; i < p.size() ; i++){
	  // pass from x and y to SUM=x+y and POS=(y-x)/SUM
	  // where y: U = 0.5*SUM*(1+POS)   x: D = 0.5*SUM*(1-POS)
    double k = (gPos[1] - gPos[0])/2;
    double d = (gPos[1] + gPos[0])/2;

    double angle = PosToAngle(p[i]);
    double slow = gELossAlphaInSi->Slow(energy*keV,(0.3)*micrometer,angle)/keV;
    //cout << energy << " " << slow << " " << angle/deg << endl ;
	  double s = slow/( 1 - parameters[4]*(k*k-pow(p[i]-d,2)) ) - parameters[1] - parameters[3]; // the ballistic deficit is subtracted here, this should lead to BD>0
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
  TGraphErrors* DataAlphaSet1 = new TGraphErrors(gDataSum1.size(), &gDataSum1[0], &gSlicePos[0], &gDataSumerr1[0], &gZeroVector[0]);
  TGraphErrors* DataAlphaSet2 = new TGraphErrors(gDataSum2.size(), &gDataSum2[0], &gSlicePos[0], &gDataSumerr2[0], &gZeroVector[0]);
  TGraphErrors* DataAlphaSet3 = new TGraphErrors(gDataSum3.size(), &gDataSum3[0], &gSlicePos[0], &gDataSumerr3[0], &gZeroVector[0]);

	DataAlphaSet1->Draw("P same");
	DataAlphaSet1->SetMarkerStyle(20);
	DataAlphaSet1->SetMarkerColor(kRed);
	DataAlphaSet1->SetLineColor(7);
	DataAlphaSet2->Draw("P same");
	DataAlphaSet2->SetMarkerStyle(20);
	DataAlphaSet2->SetMarkerColor(kBlue);
	DataAlphaSet2->SetLineColor(7);
	DataAlphaSet3->Draw("P same");
	DataAlphaSet3->SetMarkerStyle(20);
	DataAlphaSet3->SetMarkerColor(kGreen);
  DataAlphaSet3->SetLineColor(7);

	vector<double> calSum1=CalculateEnergySum(gFinalCalParam, gSlicePos, 5156.59);
	vector<double> calSum2=CalculateEnergySum(gFinalCalParam, gSlicePos, 5485.56);
	vector<double> calSum3=CalculateEnergySum(gFinalCalParam, gSlicePos, 5804.77);

	TGraphErrors* calAlphaSet1 = new TGraphErrors(calSum1.size(), &calSum1[0], &gSlicePos[0], &gZeroVector[0], &gZeroVector[0]);
	TGraphErrors* calAlphaSet2 = new TGraphErrors(calSum2.size(), &calSum2[0], &gSlicePos[0], &gZeroVector[0], &gZeroVector[0]);
	TGraphErrors* calAlphaSet3 = new TGraphErrors(calSum3.size(), &calSum3[0], &gSlicePos[0], &gZeroVector[0], &gZeroVector[0]);

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
	gSlicePos.clear();
	gZeroVector.clear();
  gPos[0]=0;
  gPos[1]=0;
	gFinalCalParam=NULL;
	gFinalCalParamError=NULL;
}
/*****************************************************************************************************************/
double PosToAngle(double pos){
	// All in mm
  double INNERBARREL_PCB_Width  = 27.76;
  double INNERBARREL_ActiveWafer_Length = 94.80;
  double INNERBARREL_ActiveWafer_Width = 24.0;
  double StripPitch = INNERBARREL_ActiveWafer_Width/4.0;

  // recaluculate the new length and the shift
  double k = (gPos[1] - gPos[0])/2;
  double d = (gPos[1] + gPos[0])/2;
  //Calculate the hit position as if it hits detector 3 (at 12 o'clock i.e. perpendiculart on the positive y-axis)
  double Z = (0.5*INNERBARREL_ActiveWafer_Length) * ((pos-d)/k);
  double Y = INNERBARREL_PCB_Width*(0.5+sin(45*deg));
  double X = (gStripNumber*StripPitch-0.5*INNERBARREL_ActiveWafer_Width)-(0.5*StripPitch);
  TVector3 HitPOS(X,Y,-Z);        // since RowPos = (U-D)/(U+D) => Downstream hit (i.e. Z>0) has RowPos<0, thus the sign
  TVector3 NormalOnDet(0,1,0);    // initiate with the normal on detector 3 at 12 o'clock (positive y-axis)
  //Rotate both vectors : Irrelevant, unless the source is not centered, so we will keep it
  HitPOS.RotateZ((3-gDetectorNumber)*45*deg);// looking downstream Detector 1 is at 3 o'clock (negative x-axis)
  NormalOnDet.RotateZ((3-gDetectorNumber)*45*deg);

  return( HitPOS.Angle(NormalOnDet) ) ;
}
/*****************************************************************************************************************/
double ApplyCalibration(double Uch, double Dch){
  // recaluculate the new length and the shift
  double k = (gPos[1] - gPos[0])/2;
  double d = (gPos[1] + gPos[0])/2;
	double pos = (Uch-Dch)/(Dch+Uch);
	double U = Uch*gFinalCalParam[0]+ gFinalCalParam[1];
	double D = Dch*gFinalCalParam[2]+ gFinalCalParam[3];
	return (U+D)*( 1 + (-gFinalCalParam[4])*(k*k-pow(pos-d,2)) );  // the BD<0 thus the negative sugn adds it here!
}
/*****************************************************************************************************************/
TCanvas* ShowCalibration(int det, int strip){

  vector <double> calEnergy,pos;

  for (unsigned int i = 0 ; i < gDataSum1.size(); i++){
     //reconstruct the U (yaxis) vs D(xaxis) plane
     double uch = 0.5 * gDataSum1[i] * (1+gSlicePos[i]); // these definitions are strict
     double dch = 0.5 * gDataSum1[i] * (1-gSlicePos[i]);
     calEnergy.push_back(ApplyCalibration(uch, dch));
     pos.push_back(gSlicePos[i]);
  }

  for (unsigned int i = 0 ; i < gDataSum2.size(); i++){
    double uch = 0.5 * gDataSum2[i] * (1+gSlicePos[i]);
    double dch = 0.5 * gDataSum2[i] * (1-gSlicePos[i]);
    calEnergy.push_back(ApplyCalibration(uch, dch));
    pos.push_back(gSlicePos[i]);
  }

  for (unsigned int i = 0 ; i < gDataSum3.size(); i++){
    double uch = 0.5 * gDataSum3[i] * (1+gSlicePos[i]);
    double dch = 0.5 * gDataSum3[i] * (1-gSlicePos[i]);
    calEnergy.push_back(ApplyCalibration(uch, dch));
    pos.push_back(gSlicePos[i]);
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
	TLine* lineE1 = new TLine (5156.6,gPos[0],5156.6,gPos[1]);
	lineE1->SetLineColor(kRed);
	lineE1->Draw();
	TLine* lineE2 = new TLine (5485.6,gPos[0],5485.6,gPos[1]);
	lineE2->Draw();
  lineE2->SetLineColor(kRed);
	TLine* lineE3 = new TLine (5804.8,gPos[0],5804.8,gPos[1]);
	lineE3->Draw();
  lineE3->SetLineColor(kRed);

  return can;
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
			cout << "*****/////***** k value for detector " << detector << " strip " << strip << " is " << position << endl;
      break;
		}// closing if statement
	} // closing while loop
	if(inputFile.eof()) cout << " Detector: " << detector << " strip: " << strip << " (DS<0>/US<1>): " << du << " are not found." << endl;
	inputFile.close();

	return position;
}
