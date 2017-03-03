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
TFile* CreateFileToCalibrate(TString alphaCalibrationFile, TString pathToMatchsticks, TString outputRootFileName); // uses the specified "alphaCalibrationFile" to make a file with histograms used to extract calibration coefficients
double fUpstream_E(double energy, unsigned short wedge, unsigned short ring); //Changes energy into linearized energy using matchstick data
double fDownstream_E(double energy, unsigned short wedge, unsigned short sector); //Changes energy into linearized energy using matchstick data


// MAIN
void CalibrateBarrelPosition(TString tripleAlphaFileName="../../../TapeData/Root/POST/ER1_1.root", 
					 TString pathToMatchsticks="../../../T40/Matchsticks/Files/Matchsticks_Calib_dummy.txt",
					 TString plotsFileName="./inspectBarrelHisto3.root"){ //tripleAlphaFileName = run file with triple alpha spectra for the Barrel in it

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
   /*



  */

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

