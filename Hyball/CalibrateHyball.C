//c++
#include <iostream>
#include <fstream>

//Root
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

//NPTool
#include "TTiaraHyballData.h"
#include "TTiaraHyballPhysics.h"
#include "NPEnergyLoss.h"
#include "NPCalibrationSource.h"
#include "NPSiliconCalibrator.h"
#include "NPCalibrationManager.h"



void CalibrateHyball(TString pathToFile="./ER192_0-nptool.root"){

//initiate output variables
vector < vector<double> > coeff; 
vector < TString > nptToken; // calibration token name in NPTool
vector < TString > badchannels; //channels with bad data

//Nptool data 
TTiaraHyballData* hyballData = new TTiaraHyballData ;

//functions
double fRing_E(double energy, unsigned short wedge, unsigned short ring);
double fSector_E(double energy, unsigned short wedge, unsigned short sector);

//spectra values
int lowerbound, upperbound, nobins;

    lowerbound = 1200;
    upperbound = 2500;
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

//initiate matchstick calibrator
  CalibrationManager* Cal  = CalibrationManager::getInstance();
    Cal->CalibrationManager::AddFile("./Matchsticks_Calib.txt");

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
NPL::EnergyLoss* ELossAlphaInAl = new NPL::EnergyLoss("He4_Si.SRIM","SRIM",10); // need to be changed for 4He

//Open TFile
TFile* nptDataFile = new TFile(pathToFile.Data(),"READ");
 if (nptDataFile == 0) {
      // if we cannot open the file, print an error message and return immediatly
      printf("Error: cannot open this file: %s \n",pathToFile.Data());
      return;
   }
nptDataFile->ls();

//Load Tree
TTree* tree = (TTree*) nptDataFile->Get("T40Tree");
//tree->Print();
//Set Branch
  tree->SetBranchStatus( "TiaraHyball" , true )               ;
  tree->SetBranchStatus( "fTiaraHyball_*" , true )               ;
  tree->SetBranchAddress( "TiaraHyball" , &hyballData )      ;

	//Loop on tree and fill the histograms
	int entries = tree->GetEntries();
	//entries = 100000;
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

TFile output("inspectHisto.root","RECREATE");
output.cd();

//Pass Histograms to the calibrator and collect the callibration coeff
vector <double> coeffset; //simple set of coeffecients
for (int iWedge =0; iWedge<6 ; iWedge++) {
	
//rings
	for(int iRing=0 ; iRing<16 ; iRing++){
		coeffset.clear();
		nameTitle = hyballRing[iWedge][iRing]->GetTitle();
		//std::cout << "Number of entries (must be >300) is: " << hyballRing[iWedge][iRing]->GetEntries() << std::endl; // used for debugging
		if (hyballRing[iWedge][iRing]->GetEntries()>100){
			double value = calibrator->SimpleCalibration(hyballRing[iWedge][iRing], alphaSource, ELossAlphaInAl, coeffset,lowerbound,upperbound);
			//std::cout << "value (must be >=0 for non-zero calibration parameters) is " << value << std::endl; used for debugging
			if (value>=0){  
				hyballRing[iWedge][iRing]->Write();
				coeff.push_back(coeffset); 
				nptToken.push_back(nameTitle); // strip's token name in NPTool
				}
			else if (value==-3){
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
		nameTitle = hyballSector[iWedge][iSector]->GetTitle();
		//std::cout << "Number of entries (must be >100) is: " << hyballSector[iWedge][iSector]->GetEntries() << std::endl; // used for debugging
		if(hyballSector[iWedge][iSector]->GetEntries()>100){
			double value = calibrator->SimpleCalibration(hyballSector[iWedge][iSector], alphaSource, ELossAlphaInAl, coeffset, lowerbound,upperbound);
			//std::cout << "value (must be >=0 for non-zero calibration parameters) is " << value << std::endl; // used for debugging
			if (value>=0){
				hyballSector[iWedge][iSector]->Write();
				coeff.push_back(coeffset); 
				nptToken.push_back(nameTitle); // strip's token name in NPTool
				}
			else if (value==-3){
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

    //Print to screen any channels with bad spectra
    if(badchannels.size() > 0){
        int badch = badchannels.size();
        std::cout << "\nWARNING: THREE PEAKS NOT FOUND IN TRIPLE-ALPHA SPECTRUM FOR CHANNELS: " << std::endl;
        for(int i ; i < badch ; i++){
            std::cout << badchannels[i] << std::endl;
        }
        std::cout << "SETTING GAIN AND OFFSET TO ZERO - BAD SPECTRUM.\nPLEASE CHECK THESE CHANNELS TO VERIFY THERE ARE NOT THREE PEAKS IN THIS SPECTRUM." << std::endl;
    }

output.Close();

//write in a text file 
	ofstream myfile;
	myfile.open ("Hyball_Calib.txt");
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

