//c++
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

using namespace std;

//Root
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TSystem.h"
#include "TMath.h"

//NPTool
#include "TTiaraHyballData.h"
#include "NPVDetector.h" // NPL::itoa fucntion
#include "NPEnergyLoss.h"
#include "NPCalibrationSource.h"
#include "NPSiliconCalibrator.h"
#include "NPCalibrationManager.h"

//functions
double fRing_E(double energy, unsigned short wedge, unsigned short ring);
double fSector_E(double energy, unsigned short wedge, unsigned short sector);


void CalibrateHyball(TString tripleAlphaFileName="/home/rw00227/nptool/Projects/T40/calibration/calibrationData/postEXPT5/ER1.root",
					           TString pathToMatchsticks="/home/rw00227/nptool/Projects/T40/calibration/Matchsticks_Calib.txt"){

//generate the outputFileName
TString plotsFileName="/home/rw00227/nptool/Projects/T40/calibration/ER1_inspectHyballHisto.root";

//define dead layer in micrometer for SimpleCalibration method,
// in this case the return "value" is the distance from the pedestal
double deadLayer, HyballOuterRadius = 135.1 /* mm */, HyballInnerRadius = 32.6 /* mm */, SourceToHyballDistance = 147 /* mm */, NumberOfStrips = 16;
//define number of iteration for the alternative ZeroExtrapolation method,
// in this case the return "value" is the calculated effective dead layer thickness
unsigned int max_iteration = 1000; // default in nptool is 10000

//initiate output variables
vector < vector<double> > coeff;
vector < TString > nptToken; // calibration token name in NPTool
vector < TString > badchannels; //channels with bad data

//Nptool data
TTiaraHyballData* hyballData = new TTiaraHyballData ;
TFile* fileToCalibrate;

//spectra values
int lowerbound, upperbound, nobins;
    lowerbound = 500; //1200
    upperbound = 1500; //2500
    nobins = (upperbound - lowerbound)/2;

//initiate list of Histograms
TH1F* hyballRing[6][16];
TH1F* hyballSector[6][8];
TString nameTitle; // same as NPTool calibration token
for (int iWedge =0; iWedge<6 ; iWedge++) {
	for(int iRing=0 ; iRing<16 ; iRing++){
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_RING%d_E",iWedge+1,iRing+1);
		hyballRing[iWedge][iRing]= new TH1F (nameTitle,nameTitle,nobins,lowerbound,upperbound);
		}
	for(int iSector=0 ; iSector<8 ; iSector++){
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_SECTOR%d_E",iWedge+1,iSector+1);
		hyballSector[iWedge][iSector]= new TH1F (nameTitle,nameTitle,nobins,lowerbound,upperbound);
		}
	}

TH1F* hyballEnergyOffsetRing; // after Linearization, this should peak at zero
TH1F* hyballEnergyOffsetSector; // after Linearization, this should peak at zero
nameTitle =Form("TIARAHYBALL_Offset_Ring");
hyballEnergyOffsetRing= new TH1F (nameTitle,nameTitle,200,-1,1);
hyballEnergyOffsetRing->GetXaxis()->SetTitle("MeV");
nameTitle =Form("TIARAHYBALL_Offset_Sector");
hyballEnergyOffsetSector= new TH1F (nameTitle,nameTitle,200,-1,1);
hyballEnergyOffsetSector->GetXaxis()->SetTitle("MeV");

//initiate matchstick calibrator
  CalibrationManager* Cal  = CalibrationManager::getInstance();
  Cal->AddFile(pathToMatchsticks.Data());
  for(int i = 0 ; i < 6 ; ++i){
    for( int j = 0 ; j < 24 ; ++j){
      Cal->AddParameter("TIARAHYBALL", "D"+NPL::itoa(i+1)+"_STRIP_RING"+NPL::itoa(j+1)+"_MATCHSTICK","TIARAHYBALL_D"+NPL::itoa(i+1)+"_STRIP_RING"+NPL::itoa(j+1)+"_MATCHSTICK")   ;
    }
    for( int j = 0 ; j < 48 ; ++j){
      Cal->AddParameter("TIARAHYBALL", "D"+NPL::itoa(i+1)+"_STRIP_SECTOR"+NPL::itoa(j+1)+"_MATCHSTICK","TIARAHYBALL_D"+NPL::itoa(i+1)+"_STRIP_SECTOR"+NPL::itoa(j+1)+"_MATCHSTICK")   ;
    }
  }
  Cal->LoadParameterFromFile();

//initiate the source, calibrator etc...
NPL::CalibrationSource* alphaSource = new NPL::CalibrationSource();
alphaSource->Set_ThreeAlphaSource();
NPL::SiliconCalibrator* calibrator = new NPL::SiliconCalibrator();
NPL::EnergyLoss* ELossAlpha = new NPL::EnergyLoss("He4_Si.SRIM","SRIM",10); // need to be changed for 4He

  if(gSystem->AccessPathName(plotsFileName)){ //checks if the file exist already, condition is "true" if not
    cout << "No file to calibrate found - creating one now using triple alpha spectra..." << endl;
	//Open TFile
	TFile* nptDataFile = new TFile(tripleAlphaFileName.Data(),"READ");
	 if (nptDataFile == 0) {
		  // if we cannot open the file, print an error message and return immediatly
		  printf("Error: cannot open this file: %s \n",tripleAlphaFileName.Data());
		  return;
	   }
	nptDataFile->ls();

	//Load Tree
	TTree* tree = (TTree*) nptDataFile->Get("T40Tree");
	//tree->Print();
	//Set Branch
	tree->SetBranchStatus( "TiaraHyball" , true )             ;
	tree->SetBranchStatus( "fTiaraHyball_*" , true )          ;
	tree->SetBranchAddress( "TiaraHyball" , &hyballData )     ;

	//Loop on tree and fill the histograms
	int entries = tree->GetEntries();
	//entries = 1000000;
	cout << " INFO: Number of entries in tree: " << entries << endl;
	for(int i = 0 ; i < entries; i++) {
	  if (i%(entries/100)) printf("\r treated %2.f percent ",100.0*i/entries);
	  tree->GetEntry(i);
	  unsigned int sizeRingE = hyballData->GetRingEMult();
	  for(unsigned int i = 0 ; i < sizeRingE ; ++i){
		unsigned short wedge = hyballData->GetRingEDetectorNbr( i );
		unsigned short ring  = hyballData->GetRingEStripNbr( i );
		double energy = hyballData->GetRingEEnergy( i );
		if( wedge>0 && ring>0 && energy>lowerbound && energy<upperbound){
	    	double linenergy = fRing_E(energy,wedge,ring);
			hyballRing[wedge-1][ring-1]->Fill(linenergy);
			//cout << wedge << " " <<  ring << " " << energy << " " << linenergy <<  endl ;
		}
	  }
	  unsigned int sizeSectorE = hyballData->GetSectorEMult();
	  for(unsigned int i = 0 ; i < sizeSectorE ; ++i){
		unsigned short wedge  = hyballData->GetSectorEDetectorNbr( i );
		unsigned short sector = hyballData->GetSectorEStripNbr( i );
		double energy = hyballData->GetSectorEEnergy( i );
		if(wedge>0 && sector>0 && energy>lowerbound && energy<upperbound ){
	    	double linenergy = fSector_E(energy,wedge,sector);
			hyballSector[wedge-1][sector-1]->Fill(linenergy);
			//cout << wedge << " " <<  sector << " " << energy << endl ;
		}
	  }
	}// end loop on tree
	nptDataFile->Close();

	//Write all the data in the file
    fileToCalibrate = new TFile(plotsFileName,"RECREATE");
	fileToCalibrate->cd();
	for (int iWedge =0; iWedge<6 ; iWedge++) {
		for(int iRing=0 ; iRing<16 ; iRing++)
			hyballRing[iWedge][iRing]->Write();
		for(int iSector=0 ; iSector<8 ; iSector++)
			hyballSector[iWedge][iSector]->Write();
		}
	fileToCalibrate->Write();
  }
	else {
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		cout << " The file " << plotsFileName << " is found in the present directory, delete or change name to create a new calibration " << endl;
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		fileToCalibrate = new TFile(plotsFileName,"UPDATE");
		//return;
		}

	fileToCalibrate->cd();

//Pass Histograms to the calibrator and collect the callibration coeff
vector <double> coeffset; //simple set of coeffecients
TH1F* currentHist;

for (int iWedge =0; iWedge<6 ; iWedge++) {
//rings
	for(int iRing=0 ; iRing<16 ; iRing++){
		coeffset.clear();
        currentHist = NULL;
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_RING%d_E",iWedge+1,iRing+1);
		currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data());
		currentHist->SetName(nameTitle+"_fit");
		cout << "Number of entries (must be >100) is:  ["<<iWedge << "] [" << iRing<< "] "<< currentHist->GetEntries() << endl; // used for debugging
		if (currentHist->GetEntries()>100){
			TString pToken = Form("TIARAHYBALL/D%d_STRIP_RING%d_MATCHSTICK",iWedge+1,iRing+1); // Matchstick token
			double pedestal = Cal->GetPedestal(pToken.Data());

			// various cout statements below used for debugging
			cout << "\nDetector number: " << iWedge+1 << endl;
			cout << "Ring number: " << iRing+1 << endl;
			deadLayer = 0.59*micrometer;
			cout << "dead layer = " << deadLayer << endl;
			cout << "ring number = " << iRing+1 << endl;
			double heightOfStrip = HyballInnerRadius+((iRing+0.5)*((HyballOuterRadius-HyballInnerRadius)/NumberOfStrips));
			cout << "height of strip = " << heightOfStrip << endl;
			double theta = atan(heightOfStrip/SourceToHyballDistance);
			cout << "theta = " << theta << " radians, which in degrees is: " << theta*(180/M_PI) << endl;
			deadLayer = deadLayer/(cos(theta));
			cout << "new dead layer = " << deadLayer << endl;

			/*double value = calibrator->ZeroExtrapolation(
				currentHist,alphaSource, ELossAlpha,
				coeffset, pedestal, max_iteration,lowerbound,upperbound);*/
			double value = calibrator->SimpleCalibration(currentHist, alphaSource, ELossAlpha, coeffset, deadLayer, lowerbound, upperbound);
			double value2 = calibrator->SimpleCalibration(currentHist, alphaSource, ELossAlpha, coeffset, /* dead layer at 0 degrees */ 0.59*micrometer, lowerbound, upperbound); // used for debugging
			cout << "value (must be >=0 for non-zero calibration parameters) is " << value << endl; //used for debugging
			cout << "value using dead layer at zero degrees = " << value2 << endl;
			if (value>=0){
                currentHist->Write("",TObject::kOverwrite);
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle); // strip's token name in NPTool
			     //cout << "gain (must be >=0 for non-zero calibration parameters) is " << coeffset[1] << endl; //used for debugging
				if(coeffset[1]>0) hyballEnergyOffsetRing->Fill(coeffset[0]);
				}
			else {
                //error code for not enough peaks in spectra
                //push channel name to vector for outputting to screen
                badchannels.push_back(nameTitle);
				if (coeffset[1]==-1){
					coeffset[1]=0;
					}
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle);
				}
			}
		else {
                badchannels.push_back(nameTitle);
			coeffset.push_back(0); coeffset.push_back(0);
			coeff.push_back(coeffset);
			nptToken.push_back(nameTitle);
			}
		}

	//sectors
	for(int iSector=0 ; iSector<8 ; iSector++){
		coeffset.clear();
        currentHist = NULL;
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_SECTOR%d_E",iWedge+1,iSector+1);
		TH1F* currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data());
		currentHist->SetName(nameTitle+"_fit");
		//cout << "Number of entries (must be >100) is:  ["<<iWedge << "] [" << iSector << "] "<< currentHist->GetEntries() << endl; // used for debugging
		if(currentHist->GetEntries()>100){
			TString pToken = Form("TIARAHYBALL/D%d_STRIP_SECTOR%d_MATCHSTICK",iWedge+1,iSector+1);
			double pedestal = Cal->GetPedestal(pToken.Data());
			/*double value = calibrator->ZeroExtrapolation(
				currentHist, alphaSource, ELossAlpha,
				coeffset, pedestal, max_iteration,lowerbound,upperbound);*/
			double value = calibrator->SimpleCalibration(currentHist, alphaSource, ELossAlpha, coeffset, deadLayer, lowerbound, upperbound);
			cout << "value (must be >=0 for non-zero calibration parameters) is " << value << endl; // used for debugging
			if (value>=0){
                currentHist->Write("",TObject::kOverwrite);
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle); // strip's token name in NPTool
			    cout << "gain (must be >=0 for non-zero calibration parameters) is " << coeffset[1] << endl; //used for debugging
				if(coeffset[1]>0) hyballEnergyOffsetSector->Fill(coeffset[0]);
				}
			else {
                //error code for not enough peaks in spectra
                //push channel name to vector for outputting to screen
                badchannels.push_back(nameTitle);
				if (coeffset[1]==-1){
					coeffset[1]=0;
					}
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle);
				}
			}
		else {
			coeffset.push_back(0); coeffset.push_back(0);
			coeff.push_back(coeffset);
			nptToken.push_back(nameTitle);
			}
		}
	}

hyballEnergyOffsetRing->Write("",TObject::kOverwrite);
hyballEnergyOffsetSector->Write("",TObject::kOverwrite);

    //Print to screen any channels with bad spectra
    if(badchannels.size() > 0){
        int badch = badchannels.size();
        cout << "\nWARNING: THREE PEAKS NOT FOUND IN TRIPLE-ALPHA SPECTRUM FOR CHANNELS: " << endl;
        for(int i=0 ; i < badch ; i++){
            cout << badchannels[i] << endl;
        }
        cout << "SETTING GAIN AND OFFSET TO ZERO - BAD SPECTRUM.\nPLEASE CHECK THESE CHANNELS TO VERIFY THERE ARE NOT THREE PEAKS IN THIS SPECTRUM." << endl;
    }

fileToCalibrate->Close();

//write in a text file
  TString CalibfName( tripleAlphaFileName( tripleAlphaFileName.Last('/')+1, tripleAlphaFileName.Length() ) );
  CalibfName.ReplaceAll("root","txt");
  CalibfName = "./calibration/Hyball_Calib_"+CalibfName;
  ofstream myfile;
  myfile.open (CalibfName.Data());
	for(unsigned int i = 0 ; i < coeff.size() ; i++){
		myfile << nptToken.at(i)<<" " ;
		for(unsigned int j = 0 ; j < coeff.at(i).size() ; j++)
			myfile << coeff.at(i).at(j)<<" " ;
		myfile<<endl;
		}
	myfile.close();

}


  double fRing_E(double energy, unsigned short wedge, unsigned short ring){
    static string name; name = "TIARAHYBALL/D" ;
    name+= NPL::itoa( wedge ) ;
    name+= "_STRIP_RING" ;
    name+= NPL::itoa( ring ) ;
    name+= "_MATCHSTICK";
    return CalibrationManager::getInstance()->ApplyCalibration(name , energy );
  }

  double fSector_E(double energy, unsigned short wedge, unsigned short sector){
    static string name; name = "TIARAHYBALL/D" ;
    name+= NPL::itoa( wedge ) ;
    name+= "_STRIP_SECTOR" ;
    name+= NPL::itoa( sector ) ;
    name+= "_MATCHSTICK";
    return CalibrationManager::getInstance()->ApplyCalibration(name , energy );
  }
