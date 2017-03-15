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

//functions
TH1F* FindHistogram(TString histname, TString filename); // finds a histogram "histname" in file "filename"
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString pathToMatchsticks, TString outputRootFileName); // uses the specified "alphaCalibrationFile" to make a file with histograms used to extract calibration coefficients
double fUpstream_E(double energy, unsigned short wedge, unsigned short ring); //Changes energy into linearized energy using matchstick data
double fDownstream_E(double energy, unsigned short wedge, unsigned short sector); //Changes energy into linearized energy using matchstick data
double ErfFunction(Double_t *x, Double_t *par); // amn Erf function to fit the positions
TF1* FitPosition(TH1F* hist, int du, double rmin, double rmax); 
void GetFitWindow(TH1F* hist, int du, double &rmin, double &rmax); // in progress

// MAIN
void CalibrateBarrelPosition(TString tripleAlphaFileName/*to avoid conflict input the file name in the terminal*/, 
					 TString pathToMatchsticks="../../../T40/Matchsticks/Files/Matchsticks_Calib_dummy.txt",
					 TString plotsFileName="./inspectBarrelHisto2.root"){ //tripleAlphaFileName = run file with triple alpha spectra for the Barrel in it

  //local variable
  ofstream outputFile;
  outputFile.open("Barrel_Calib_Pos.txt");
  TCanvas* can[8]; // initialises 8 canvases; 1 for each Barrel detector element

  TFile* fileToCalibrate;
  if(gSystem->AccessPathName(plotsFileName)){ //checks if the file exist already, condition is "true" if not
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
    cout << " No file to calibrate found - creating one now using triple alpha spectra..." << endl;
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
    fileToCalibrate = CreateFileToCalibrate(tripleAlphaFileName, pathToMatchsticks, plotsFileName);
	fileToCalibrate->Write();
  }
  else {
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
	cout << " The file " << plotsFileName << " is found in the present directory, it will be used for position calibration " << endl;
	cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;			
	fileToCalibrate = new TFile(plotsFileName,"UPDATE");
  }
  fileToCalibrate->cd();

  // Calibration of position goes here 
	for (int detector=1; detector<=8; detector++){
		TString name = Form("Barrel%d",detector);
		can[detector-1] = new TCanvas(name,name,900,900);
		can[detector-1]->Divide(4,4); // splits each canvas into 8; 2 for each strip in a detector
		for (int strip=1; strip<=4; strip++){
			for (int du=0; du<2; du++){
				TString hname = Form("TIARABARREL_B%d_DOWNSTREAM%d_P",detector,strip); 
				if (du == 1) hname = Form("TIARABARREL_B%d_UPSTREAM%d_P",detector,strip);
				
         //ClearGlobalParameters(); // clear parameters
				TH1F* currentHist = (TH1F*) fileToCalibrate->FindObjectAny(hname.Data());
				can[detector-1]->cd(2*(strip-1)+(du+1));

				if (currentHist && currentHist->Integral()>2000) {
					cout << " Working on " << hname << endl;
          // Find the good window
					double rmin=0.5;
          double rmax=0.8;
		      //FitWindow(TH1F* hist, du, rmin, rmax); // this will give rmin and rmax
		      TF1* fitfunc = FitPosition(currentHist, du, rmin, rmax);
          //Draw for inspection
          currentHist->Draw();
		      //write in output file 
				  if (du == 0)
							outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_POS " << 0 << " " << fitfunc->GetParameter(1) << endl;
					else
							outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_POS " << 0 << " " << fitfunc->GetParameter(1) << endl;
				}
				else { // for when there is no histogram currentHist or if currentHist is empty - nptool tokens with default calibration parameters
				  if (du == 0)
							outputFile << "TIARABARREL_B" << detector << "_DOWNSTREAM" << strip << "_POS " << 0 << " " << 0 << endl;
					else
							outputFile << "TIARABARREL_B" << detector << "_UPSTREAM" << strip << "_POS " << 0 << " " << 0 << endl;
				}
			}
		}
	}

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
      Cal->AddParameter("TIARABARREL","MATCHSTICK_B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_E","TIARABARREL_MATCHSTICK_B"+NPL::itoa(i+1)+"_UPSTREAM"+NPL::itoa(j+1)+"_E")   ;
      Cal->AddParameter("TIARABARREL","MATCHSTICK_B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_E","TIARABARREL_MATCHSTICK_B"+NPL::itoa(i+1)+"_DOWNSTREAM"+NPL::itoa(j+1)+"_E")   ;
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
	      barrelFrontStripP[iSide][iStrip][1]= new TH1F (nameTitle,nameTitle,400,+0.2,+0.8);
		  nameTitle =Form("TIARABARREL_B%d_DOWNSTREAM%d_P",iSide+1,iStrip+1);
		  barrelFrontStripP[iSide][iStrip][0]= new TH1F (nameTitle,nameTitle,400,-0.8,-0.2);
		}
	  for(int iStrip=0 ; iStrip<4 ; iStrip++){
		  nameTitle =Form("TIARABARREL_B%d_UD%d_E",iSide+1,iStrip+1);
		  barrelFrontStripUD[iSide][iStrip]= new TH2F (nameTitle,nameTitle,1600,-50,1550,600,-50,1550);
		  nameTitle =Form("TIARABARREL_B%d_PE%d_E",iSide+1,iStrip+1); // the ones we're interested in making
		  barrelFrontStripPE[iSide][iStrip]= new TH2F (nameTitle,nameTitle,1600,-50,1550,500,-1,+1);
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
				  barrelFrontStripUD[sideU-1][stripU-1]->Fill(energyD,energyU);
				  barrelFrontStripPE[sideU-1][stripU-1]->Fill(energyU+energyD,(energyU-energyD)/(energyD+energyU));
				  if(energyU+energyD>650 && (energyD/(energyU+energyD)>0.60))
				    barrelFrontStripP[sideU-1][stripU-1][0]->Fill((energyU-energyD)/(energyU+energyD));
				  if(energyU+energyD>650 && (energyU/(energyU+energyD)>0.60))
				    barrelFrontStripP[sideU-1][stripU-1][1]->Fill((energyU-energyD)/(energyU+energyD)); 
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


double ErfFunction(Double_t *x, Double_t *par) {
  // Par[0]: amplitude measured from the y = 0, amp is positive for "S" shape 
  // Par[1]: inflexion point position, (what we are looking for) 
  // Par[2]: steepness of the curve, (the smaller the steeper) 

   double xx = x[0];
   double f = (par[0]/2)*TMath::Erf((xx-par[1])/par[2]) + TMath::Abs(par[0]/2);
   return f;
}

TH1F* FindHistogram(TString HistoName, TString filename){

  TFile* file = new TFile(filename.Data());
  TH1F* currentHist = (TH1F*) file->FindObjectAny(HistoName.Data());

  return ((TH1F*) currentHist);
}


TF1* FitPosition(TH1F* hist, int du, double rmin, double rmax){
//for du = {DownStream,UpStream} = {0,1}  => sign={+,-}
  
  hist->Rebin();

  //Set sign for Upstream by default, rmin, rmax are the same
  double sign = -1 ;
  // if downstream change sign, swap and rmin, rmax
  if (du==0){ 
      double rmincopy=rmin;
      sign = +1 ;
      rmin = -rmax;
      rmax = -rmincopy;
   }

   double amp = sign*50 ; 
   double pos = -sign*0.66 ; 
   double steepness = 0.03 ; 

    TF1* fitFunction = new TF1("fitFunction", ErfFunction ,rmin, rmax, 3); // 3 in nb of param
		fitFunction->SetParameters(amp,pos,steepness);

		fitFunction->SetParLimits(0,0,amp*2);
		fitFunction->SetParLimits(1,pos-0.1,pos+0.1);
		fitFunction->SetParLimits(2,0.005,steepness*3);
		//fitFunction->FixParameter(2,steepness);

		hist->Fit(fitFunction,"RM");
    hist->SetLineColor(8);
    TF1* fit = hist->GetFunction(fitFunction->GetName());

 return fit;
}

