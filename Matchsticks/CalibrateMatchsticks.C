// Run using ".x CalibrateMatchsticks.C+"
// Run for every file with matchsticks data.


//c++
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
using namespace std;

//Root
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
//NPTool
#include "TTiaraHyballData.h"
#include "TTiaraBarrelData.h"


//functions
double MatchstickCalibration(TH1* histo, vector<double>& coeff);

//Main 
void CalibrateMatchsticks(TString pathToFile/*to avoid conflict input the file name in the terminal*/){ 

//initiate output variables
vector < vector<double> > coeff; 
vector < TString > nptToken; // calibration token name in NPTool
vector < TString > badchannels; //channels with bad data

//Nptool data 
TTiaraHyballData* hyballData = new TTiaraHyballData ;
TTiaraBarrelData* barrelData = new TTiaraBarrelData ;

//initiate list of Histograms
TH1F* barrelFrontStrip[8][4][2]; // 8 sides of the barrel, 4 strips in one side, up and down  
TH1F* hyballRing[6][16];
TH1F* hyballSector[6][8]; 
TString nameTitle; // same as NPTool calibration token

//create hists for barrel info
for (int iSide =0; iSide<8 ; iSide++) {
	for(int iStrip=0 ; iStrip<4 ; iStrip++){
		nameTitle =Form("TIARABARREL_MATCHSTICK_B%d_UPSTREAM%d_E",iSide+1,iStrip+1);
		barrelFrontStrip[iSide][iStrip][0]= new TH1F (nameTitle,nameTitle,4050,0,4050);
		nameTitle =Form("TIARABARREL_MATCHSTICK_B%d_DOWNSTREAM%d_E",iSide+1,iStrip+1);
		barrelFrontStrip[iSide][iStrip][1]= new TH1F (nameTitle,nameTitle,4050,0,4050);
		}
}

//create hist for hyball info
for (int iWedge =0; iWedge<6 ; iWedge++) {
	for(int iRing=0 ; iRing<16 ; iRing++){
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_RING%d_MATCHSTICK",iWedge+1,iRing+1);
		hyballRing[iWedge][iRing]= new TH1F (nameTitle,nameTitle,4050,0,4050);
		}
	for(int iSector=0 ; iSector<8 ; iSector++){
		nameTitle =Form("TIARAHYBALL_D%d_STRIP_SECTOR%d_MATCHSTICK",iWedge+1,iSector+1);
		hyballSector[iWedge][iSector]= new TH1F (nameTitle,nameTitle,4050,0,4050);
		}
	}


TString OutputfName( pathToFile( pathToFile.Last('/')+1, pathToFile.Length() ) );
OutputfName = "inspect_"+OutputfName;

  TFile* fileToCalibrate;
  if(gSystem->AccessPathName(OutputfName)){ //checks if the file exist already, condition is "true" if not
    cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
    cout << " No file to calibrate found - creating one now using " << pathToFile << endl;
    cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;

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

    //Set Barrel Branch
    tree->SetBranchStatus( "TiaraBarrel" , true )         ;
    tree->SetBranchStatus( "fTiaraBarrel_*" , true )      ;
    tree->SetBranchAddress( "TiaraBarrel" , &barrelData ) ;

    //Loop on tree and fill the histograms
    int entries = tree->GetEntries();
    //entries = 100000;
    cout << " INFO: Number of entries in Barrel tree: " << entries << endl;  
    for(int i = 0 ; i < entries; i++) { 
      if (i%(entries/100)) printf("\r treated %2.f percent of Barrel entries",100.0*i/entries);
      tree->GetEntry(i);
      unsigned int sizeUE = barrelData->GetFrontUpstreamEMult();
      unsigned int sizeDE = barrelData->GetFrontDownstreamEMult();

      // fill upstream 
      for(unsigned int i = 0 ; i < sizeUE ; ++i){
        unsigned short side = barrelData->GetFrontUpstreamEDetectorNbr( i );
        unsigned short strip  = barrelData->GetFrontUpstreamEStripNbr( i );
        double energy = barrelData->GetFrontUpstreamEEnergy( i );
        if( side>0 && strip>0 && energy>0 && energy<4050 ){
          barrelFrontStrip[side-1][strip-1][0]->Fill(energy);
          //cout << side << " " <<  strip << " " << energy << endl ; 
          }
        }

      // fill downstream
      for(unsigned int i = 0 ; i < sizeDE ; ++i){
        unsigned short side = barrelData->GetFrontDownstreamEDetectorNbr( i );
        unsigned short strip  = barrelData->GetFrontDownstreamEStripNbr( i );
        double energy = barrelData->GetFrontDownstreamEEnergy( i );
        if( side>0 && strip>0 && energy>0 && energy<4050 ){
          barrelFrontStrip[side-1][strip-1][1]->Fill(energy);
          //cout << side << " " <<  strip << " " << energy << endl ; 
         }
      }
     }

      //Set Hyball Branch
      tree->SetBranchStatus( "TiaraHyball" , true )               ;
      tree->SetBranchStatus( "fTiaraHyball_*" , true )               ;
      tree->SetBranchAddress( "TiaraHyball" , &hyballData )      ;

      for(int i = 0 ; i < entries; i++) {
        if (i%(entries/100)) printf("\r treated %2.f percent of Hyball entries",100.0*i/entries);
        tree->GetEntry(i);
        //fill ring information
        unsigned int sizeRingE = hyballData->GetRingEMult();
        for(unsigned int i = 0 ; i < sizeRingE ; ++i){
          unsigned short wedge = hyballData->GetRingEDetectorNbr( i );
          unsigned short ring  = hyballData->GetRingEStripNbr( i );
          double energy = hyballData->GetRingEEnergy( i );
          if( wedge>0 && ring>0 && energy>0 && energy<4050 ){
            hyballRing[wedge-1][ring-1]->Fill(energy);
            //cout << wedge << " " <<  ring << " " << energy << endl ; 
          }
        }
      //fill sector information
      unsigned int sizeSectorE = hyballData->GetSectorEMult();
      for(unsigned int i = 0 ; i < sizeSectorE ; ++i){
        unsigned short wedge  = hyballData->GetSectorEDetectorNbr( i );
        unsigned short sector = hyballData->GetSectorEStripNbr( i );
        double energy = hyballData->GetSectorEEnergy( i );
        if(wedge>0 && sector>0 && energy>0 && energy<4050 ){
          hyballSector[wedge-1][sector-1]->Fill(energy);
          //cout << wedge << " " <<  sector << " " << energy << endl ; 
        }
      }
    }// end loop on tree
    nptDataFile->Close();
    
    fileToCalibrate = new TFile(OutputfName,"RECREATE");
    //write all histograms
    for (int iSide =0; iSide<8 ; iSide++) {
	    for(int iStrip=0 ; iStrip<4 ; iStrip++){
		    barrelFrontStrip[iSide][iStrip][0]->Write();
		    barrelFrontStrip[iSide][iStrip][1]->Write();
		    }
    }
    for (int iWedge =0; iWedge<6 ; iWedge++) {
	    for(int iRing=0 ; iRing<16 ; iRing++){
		    hyballRing[iWedge][iRing]->Write();
		    }
	    for(int iSector=0 ; iSector<8 ; iSector++){
		    hyballSector[iWedge][iSector]->Write();
		    }
	    }
    fileToCalibrate->Write();
  }
  else {
	  cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
	  cout << " The file " << OutputfName << " is found in the present directory, it will be used for calibration " << endl;
	  cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
	  fileToCalibrate = new TFile(OutputfName,"UPDATE");
  }

fileToCalibrate->cd();

vector <double> coeffset; //simple set of coeffecients
TH1F* currentHist;

for (int iSide =0; iSide<8 ; iSide++) {
	for(int iStrip=0 ; iStrip<4 ; iStrip++){
    for(int iUD=0 ; iUD<2 ; iUD++){
      currentHist = NULL;
	    coeffset.clear();
	    TString nameTitle = Form("TIARABARREL_MATCHSTICK_B%d_UPSTREAM%d_E",iSide+1,iStrip+1);
	    if(iUD==1) nameTitle = Form("TIARABARREL_MATCHSTICK_B%d_DOWNSTREAM%d_E",iSide+1,iStrip+1);
	    currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data())->Clone();
  		currentHist->SetName(nameTitle+"_fit");
	    //cout << "Number of entries (must be >300) is: " << currentHist->GetEntries() << endl; // used for debugging
	    if (currentHist->GetEntries()>300){
         //pass histo to co-efficient finding function
		    double value = MatchstickCalibration(currentHist, coeffset);
		    //cout << "value (must be >=0 for non-zero calibration parameters) is " << value << endl; used for debugging
		    if (value>=0){  
			    currentHist->Write("",TObject::kOverwrite);//save hist to output file
			    coeff.push_back(coeffset); 
			    nptToken.push_back(nameTitle); // strip's token name in NPTool
			    }
		    else 
		      if (value==-3){
            //error code for not enough peaks in spectra
            //push channel name to vector for outputting to screen
            badchannels.push_back(nameTitle);
			      if (coeffset[1]==-1){
				      coeffset[1]=1;
				    }			
			      coeff.push_back(coeffset);
			      nptToken.push_back(nameTitle);
			    }
		  }
		    else {
			    coeffset.push_back(0); coeffset.push_back(1); coeffset.push_back(0);				
			    coeff.push_back(coeffset);
			    nptToken.push_back(nameTitle);
		  	}
      }
    }
}


for (int iWedge =0; iWedge<6 ; iWedge++) {
//rings
	for(int iRing=0 ; iRing<16 ; iRing++){
    currentHist = NULL;
    coeffset.clear();
    TString nameTitle = Form("TIARAHYBALL_D%d_STRIP_RING%d_MATCHSTICK",iWedge+1,iRing+1);
	  currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data())->Clone();
    currentHist->SetName(nameTitle+"_fit");		    
		//cout << "Number of entries (must be >300) is: " << currentHist->GetEntries() << endl; // used for debugging
		if (currentHist->GetEntries()>300){
      //pass histo to co-efficient finding function
			double value = MatchstickCalibration(currentHist, coeffset);
			//cout << "value (must be >=0 for non-zero calibration parameters) is " << value << endl; used for debugging
			if (value>=0){  
				currentHist->Write("",TObject::kOverwrite);//save hist to output file
				coeff.push_back(coeffset); 
				nptToken.push_back(nameTitle); // strip's token name in NPTool
				}
			else 
			  if (value==-3){
          //error code for not enough peaks in spectra
          //push channel name to vector for outputting to screen
          badchannels.push_back(nameTitle);
				  if (coeffset[1]==-1){
					  coeffset[1]=1;
					}			
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle);
				}
			}
		else {
			coeffset.push_back(0); coeffset.push_back(1); coeffset.push_back(0);				
			coeff.push_back(coeffset);
			nptToken.push_back(nameTitle);
			}
		}

	//sectors
	for(int iSector=0 ; iSector<8 ; iSector++){
	  currentHist = NULL;
    coeffset.clear();
    TString nameTitle = Form("TIARAHYBALL_D%d_STRIP_SECTOR%d_MATCHSTICK",iWedge+1,iSector+1);
	  currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data())->Clone();
    currentHist->SetName(nameTitle+"_fit");	
		//cout << "Number of entries (must be >100) is: " << currentHist->GetEntries() << endl; // used for debugging
		if(currentHist->GetEntries()>100){
      //pass histo to co-efficient finding function
			double value = MatchstickCalibration(currentHist, coeffset);
			//cout << "value (must be >=0 for non-zero calibration parameters) is " << value << endl; // used for debugging
			if (value>=0){
				currentHist->Write("",TObject::kOverwrite);//save hist to output file
				coeff.push_back(coeffset); 
				nptToken.push_back(nameTitle); // strip's token name in NPTool
				}
			else 
			if (value==-3){
        //error code for not enough peaks in spectra
        //push channel name to vector for outputting to screen
        badchannels.push_back(nameTitle);
				if (coeffset[1]==-1){
					coeffset[1]=1;
				}			
				coeff.push_back(coeffset);
				nptToken.push_back(nameTitle);
			}
		}
		else {
			coeffset.push_back(0); coeffset.push_back(1); coeffset.push_back(0);				
			coeff.push_back(coeffset);
			nptToken.push_back(nameTitle);
			}
		}
	}

    //Print to screen any channels with bad spectra
  if(badchannels.size() > 0){
    int badch = badchannels.size();
    cout << "\nWARNING: LESS THAN 7 PEAKS FOUND IN MATCHSTICKS SPECTRUM FOR CHANNELS: " << endl;
    for(int i=0 ; i < badch ; i++){
      cout << badchannels[i] << endl;
    }
    cout << "\nSETTING COEFFICIENTS TO ZERO - BAD SPECTRA.\nPLEASE CHECK THIS CHANNEL TO VERIFY THERE ARE SUFFICIENT PEAKS IN THIS SPECTRA." << endl;
  }

  fileToCalibrate->Close();

//write in a text file 
  TString CalibfName( pathToFile( pathToFile.Last('/')+1, pathToFile.Length() ) );
  CalibfName.ReplaceAll("root","txt");
  CalibfName = "Matchsticks_Calib_"+CalibfName;
  fstream fileopen;
  fileopen.open(CalibfName.Data());

  //Check whether there is currently a file with data of same type.
  bool open = true;
  if ( fileopen.peek() == ifstream::traits_type::eof() ){
      open = false;
  }
  fileopen.close();
  
  //if file has no data create a file and input all new data.
  if ( open == false ){
    ofstream myfile;
    myfile.open (CalibfName.Data());

      cout << "No file found; creating file "<< CalibfName.Data() << endl;
    for(unsigned int i = 0 ; i < coeff.size() ; i++){
	    myfile << nptToken.at(i)<<" " ;
	    for(unsigned int j = 0 ; j < coeff.at(i).size() ; j++)
		    myfile << coeff.at(i).at(j)<<" " ;
	        myfile<<endl;
          
	}
    myfile.close();}

  //if there is a file take data within file into consideration.
  else{
      //open file with data to be read
    ifstream readfile;
    readfile.open (CalibfName.Data());
    cout << "File " << CalibfName << " found; including data from file" << endl;
    string line; //element to be added to vectors
    vector<vector<string>> readvalues; //2D vector for data from text file
    
    //take information line by line from file and insert into 2D vector by column
    while(getline(readfile,line)){
      //some elements to create 2D vector
      vector<string> values;
      stringstream ls(line);
      string value;
      while(ls>>value){
          values.push_back(value);
      }
      if (values.size())
          readvalues.push_back(values);
    }
    readfile.close();

    ofstream myfile;
    myfile.open (CalibfName.Data());
    //below sections makes a judgement whether to write to new file data
    //from the previous file or newly calculated coefficients.
    double element1, element2;
    for(unsigned int j = 0 ; j < coeff.size() ; j++){
      if(readvalues[j][0] == nptToken[j]){
        //convert the data from file from string to double for comparison to new data.
        element1 = stod(readvalues[j][1]);
        element2 = stod(readvalues[j][2]);
        if((element1 == 0 && element2 == 1) || (coeff[j][1] != 0 && coeff[j][2] != 0)){
          myfile << nptToken.at(j) << " " ;
          for(unsigned int k = 0 ; k < coeff.at(j).size() ; k++)
            myfile << coeff.at(j).at(k)<<" " ;
          myfile<<endl;                
        }
        else{
          for(unsigned int k = 0 ; k < readvalues.at(j).size() ; k++)
            myfile << readvalues.at(j).at(k)<<" " ;
          myfile<<endl;
        }
      }
      else{
        cout << "Something wrong with current " << CalibfName << " file. File modified preventing match between tokens. Please either correct or delete current calib file." << endl;
      }
    }
    myfile.close(); 
  }
}


double MatchstickCalibration(TH1* histo, vector<double>& coeff){

   //adjust the histogram range to exclude first peak with high stat (noise)
   int nbin = histo->GetNbinsX();
   int bin200 = histo->FindBin(200); // 200 is an upper limit for noise peak
   int bin1000 = histo->FindBin(1000); // 200 is an upper limit for noise peak
   //cout << bin1000 << " " << histo->Integral(1,bin1000) << " " << histo->Integral(bin1000,nbin-1)  << endl ; 
   if(histo->Integral(bin1000,nbin)<500){ // nopeaks at higher channels, no need to proceed
    cout << "  nopeaks at higher channels, skipping " << histo->GetName() << endl; 
    coeff.clear();
    coeff.push_back(0);
    coeff.push_back(-1);
    coeff.push_back(0);
    return -2;
   }
  
  double xmax  = histo->GetXaxis()->GetXmax();
  double xmin  = histo->GetXaxis()->GetXmin();
  // Get the bin with maximum content below 200
  histo->GetXaxis()->SetRangeUser(xmin,200);
  int binlow  = histo->GetMaximumBin();
  double binclow  = histo->GetBinContent(binlow);
  // Get the bin with maximum content above 200
  histo->GetXaxis()->SetRangeUser(200,xmax);
  int binhigh  = histo->GetMaximumBin();
  double binchigh  = histo->GetBinContent(binhigh);

  // if the the noise peak is 3 times higher, readjust  
   if (binclow > 3*binchigh){ 
      //cout << " new x low " << histo->GetBinCenter(binlow)+50 << endl ; 
      histo->GetXaxis()->SetRangeUser(histo->GetBinCenter(binlow)+50, xmax);
   }
   //find peaks within the spectra using search from TSpectrum.h
   TSpectrum *matchstick_spectra = new TSpectrum(20);
   int PeakNbr = matchstick_spectra->Search(histo, 3,"", 0.001);
   //printf("Found %d candidate peaks to fit\n",PeakNbr);

   if(PeakNbr>7){

        //take channel number of peaks from tspectra
        Double_t *xpeaks = matchstick_spectra->GetPositionX();

        //rough peak position using the search function
        vector<double> PeakPosition;
        for(int i=0 ; i<PeakNbr ; i++){
            PeakPosition.push_back(xpeaks[i]);
        }

        //put peaks into ascending order of channel number
        sort(PeakPosition.begin(),PeakPosition.end());

        //vector for the channel number of peaks
        vector<double> PeakFitPosition;
        
        //Fit all the peaks in the spectra with gaussian
        for(int j=0 ; j<PeakNbr; j++){
            TF1* gausfit = new TF1("gausfit", "gaus(0)", PeakPosition[j]-30, PeakPosition[j]+30);
            histo->Fit(gausfit,"QR+");
            PeakFitPosition.push_back(gausfit->GetParameter(1));            
        }

        //Attempt to remove any noise peak at lowest channels
        //May need adjusting dependant on data
        //if(PeakFitPosition[10]-PeakFitPosition[9]<200 && PeakFitPosition[11]-PeakFitPosition[0]>500 && PeakFitPosition[9]-PeakFitPosition[8]<200){
        //PeakFitPosition.erase(PeakFitPosition.begin()+0);}


        //vector for the pulser voltages in mV, may need adjusting depending on the
        //physical input from the pulser
        vector<double> MatchstickVoltage;
        MatchstickVoltage.push_back(100);
        MatchstickVoltage.push_back(200);
        MatchstickVoltage.push_back(300);
        MatchstickVoltage.push_back(400);
        MatchstickVoltage.push_back(500);
        MatchstickVoltage.push_back(600);
        MatchstickVoltage.push_back(700);
        MatchstickVoltage.push_back(800);
        MatchstickVoltage.push_back(900);
        MatchstickVoltage.push_back(1000);
        MatchstickVoltage.push_back(1500);
        MatchstickVoltage.push_back(2000);
        MatchstickVoltage.push_back(2500);
        MatchstickVoltage.push_back(3000);
        MatchstickVoltage.push_back(3500);
        MatchstickVoltage.push_back(4000);
        MatchstickVoltage.push_back(4500);
        MatchstickVoltage.push_back(5000);
        MatchstickVoltage.push_back(5500);
        MatchstickVoltage.push_back(6000);
        MatchstickVoltage.push_back(6500);
        MatchstickVoltage.push_back(7000);
        MatchstickVoltage.push_back(7500);
        MatchstickVoltage.push_back(8000);
        MatchstickVoltage.push_back(9000);
        MatchstickVoltage.push_back(10000);

        //remove any pulser voltage not associated to a peak in the spectra
        //int MatchstickEndVoltage = PeakFitPosition.size();
        //MatchstickVoltage.erase(MatchstickVoltage.begin() + (MatchstickEndVoltage-1), MatchstickVoltage.end());
        MatchstickVoltage.resize(PeakFitPosition.size());
        
        //plot pulser voltage vs channel number, fit to determine function of the 
        // ADC's non-linearity and save graph to output root file.
        TGraph* graph = new TGraph(PeakFitPosition.size(), &PeakFitPosition[0], &MatchstickVoltage[0]);
        graph->GetXaxis()->SetTitle("Channel Number");
        graph->GetXaxis()->CenterTitle();
        graph->GetYaxis()->SetTitle("Pulser Voltage (mV)");
        graph->GetYaxis()->CenterTitle(); 
        graph->SetMarkerStyle(20); 

        //Ignore anything higher than 3980 from fit to avoid warping
        int peakfitsize = PeakFitPosition.size() - 1;
        int fitbegin = PeakFitPosition[0];
        int fitend = PeakFitPosition[peakfitsize];
        if(PeakFitPosition[peakfitsize]>3850){
            fitend = PeakFitPosition[(peakfitsize-1) ];}
        else{fitend = PeakFitPosition[peakfitsize];}

        TF1 *polfit = new TF1("polfit","pol2(0)",fitbegin,fitend);   

        graph->Fit(polfit,"QR");
        gStyle -> SetOptStat(0);
        gStyle -> SetOptFit(111);
        TString graphname = TString(histo->GetTitle())+"_gr";
        graph->Write(graphname.Data(),TObject::kOverwrite);

        //retrieve fit co-efficients and push into vector to be returned to main
        double polfit0, polfit1, polfit2;
        polfit0 = polfit->GetParameter(0);
        polfit1 = polfit->GetParameter(1);
        polfit2 = polfit->GetParameter(2);

        coeff.clear();
        coeff.push_back(polfit0);
        coeff.push_back(polfit1);
        coeff.push_back(polfit2);


        //rest is for debugging. produces residual voltages (fitted pulser voltage
        //vs acutal pulser voltage as a check of fit quality
        //saves graph to the output root file.
        /*vector<double> FitMatchstickVoltage;
        vector<double> ResidualVoltage;

        for(int k=0; k<MatchstickEndVoltage; k++){
        FitMatchstickVoltage.push_back((polfit2*(PeakFitPosition[k]*PeakFitPosition[k])) + (PeakFitPosition[k]*polfit1)+polfit0);
        ResidualVoltage.push_back(FitMatchstickVoltage[k]-MatchstickVoltage[k]);
        }


        TGraph* graph2 = new TGraph(PeakFitPosition.size(), &PeakFitPosition[0], &ResidualVoltage[0]);
        graph2->GetXaxis()->SetTitle("Channel Number");
        graph2->GetXaxis()->CenterTitle();
        graph2->GetYaxis()->SetTitle("Residual Voltage (mV)");
        graph2->GetYaxis()->CenterTitle();
        graph2->Fit("pol2","Q");
        graph2->Write(histo->GetTitle());*/

        //return value to prove coefficients have been found
        return(abs(polfit1/polfit2));}

    //if not enough peaks found return error code to main
    else{
        coeff.clear();
        coeff.push_back(0);
        coeff.push_back(-1);
        coeff.push_back(0);
        return -3;}

}



