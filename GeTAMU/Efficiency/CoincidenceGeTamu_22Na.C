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
#include "TGeTAMUPhysics.h"
#include "TTiaraBarrelPhysics.h"
#include "TTiaraHyballPhysics.h"
#include "NPVDetector.h" // NPL::itoa fucntion
#include "NPEnergyLoss.h"
#include "NPCalibrationSource.h"
#include "NPSiliconCalibrator.h"
#include "NPCalibrationManager.h"

#include "TFrame.h"

double gausfit(double *x, double *par)
{
float xx = x[0];
double f = par[0]*exp(-0.5*( (xx-par[1])*(xx-par[1]) ) / (par[2]*par[2]) ) + par[3];
f+= par[3] + par[4] * (xx - par[1]) + par[5] * (xx - par[1]) * (xx - par[1]);

return f;
}

//////////////////////////////////// READ ME FIRST ////////////////////////////////////////////////////////
/*
	This code is the first (1/2) step to get the absolute efficiency cuver for GeTAMU.
	Here we determine the absolute efficiency at one specific energy, which will be used to
	normalize other efficiency data (from 152Eu).
	For this purpose, you can use 60Co spectrum (1173 keV and 1332 keV) or 22Na spectrum 
	(two 511 keV). For 22Na spectrum, you need to set AngularCorrelation == false.

	Also NOTE! Root files are all .root after npanalysis (not midas2nptool).
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void CoincidenceGeTamu_22Na(TString Na22FileName="/home/shuyaota/nptool/Outputs/Analysis/EXPT6_ER74_0.root", bool AngularCorrelation=false)
{

	//generate the outputFileName
	TString plotsFileName( Na22FileName( Na22FileName.Last('/')+1, Na22FileName.Length() ) );
	plotsFileName= "inspect_"+plotsFileName;

	//initiate output variables
	vector < TString > badchannels; //channels with bad data

	//Nptool data
	TGeTAMUPhysics* getamuPhysics = new TGeTAMUPhysics() ;
	TFile* fileToCalibrate;


//////////////////////////////////////////////Angular Correlation Setting///////////////////////////////////////////////////////////////////

	double	CoincidenceHit[4][4];	// Coincidence hit of two gammas.
	double	InitialHit[2][4];	// Number of gammas hit the clover without coincidence.
	double	AngularCor[4][4];	// Angular Correlation
	double	GeometricalFactor[4][4];	// Almost same as Clover's geometrical coverage (i.e., 10% of 4pi, etc).
	double	SourceParticles = 1e7;	//Created gammas in the following simulation.

	if(AngularCorrelation)
	{
		//Angular correlation for two E2 gamma transitions (Thanks to Momo).
		//Assuming 1e7 gammas created at the target position. Clover1-4 are placed at the (measured) position of GeTAMU setup in T40 campaign.
		//For 60Co case, you can consider A: 1173 keV gamma, B: 1332 keV gamma. (since technically, 1173 keV is emitted first (2505keV --> 1332 --> GS), but frankly no difference).
		//By considering B = 1332 keV gamma, you measure A = 1173 keV gamma in coincidence with B, and determine the efficiency of A in a clover.
/*
	// This is a table.  
	 	A in Clover 1 	A in Clover 2 	A in Clover 3 	A in Clover 4 
	While B in Clover 1 	113102 	94608 	125177 	94608
	While B in Clover 2 	95001 	91928 	104257 	91928
	While B in Clover 3 	125294 	104356 	137457 	104356
	While B in Clover 4 	103644 	99496 	114182 	99496 

	// Singles. For example, assuming the first gamma is emitted isotropically, then A in Clover 1 / Source Particles ~ Geometrical efficiency of the Clover 1.
	    A in Clover 1 = 1.04039e+06

	    A in Clover 2 = 932421

	    A in Clover 3 = 1.14775e+06

	    A in Clover 4 = 1.01684e+06

	    B in Clover 1 = 1.04062e+06

	    B in Clover 2 = 932696

	    B in Clover 3 = 1.14724e+06

	    B in Clover 4 = 1.01729e+06
*/
		//CoincidenceHit[X][Y].	I.e., While 1173 keV gamma is detected in clover B, "efficiency of 1332 keV gamma (A)" in the clover A will be....
		//Basic equation is: Num (1332) at Clover X / Num (1173) at Clover B = eff (1332) at Clover A * CoincidenceHit[B][A] / 
		CoincidenceHit[0][0] = 113102;	//While B in Clover1, A in Clover1
		CoincidenceHit[0][1] = 94608;	//While B in Clover1, A in Clover2
		CoincidenceHit[0][2] = 125177;	//While B in Clover1, A in Clover3
		CoincidenceHit[0][3] = 94608;	//While B in Clover1, A in Clover4
		CoincidenceHit[1][0] = 95001;	//While B in Clover2, A in Clover1
		CoincidenceHit[1][1] = 91928;	//While B in Clover2, A in Clover2
		CoincidenceHit[1][2] = 104257;	//While B in Clover2, A in Clover3
		CoincidenceHit[1][3] = 91928;	//While B in Clover2, A in Clover4
		CoincidenceHit[2][0] = 125294;	//While B in Clover3, A in Clover1
		CoincidenceHit[2][1] = 104356;	//While B in Clover3, A in Clover2
		CoincidenceHit[2][2] = 137457;	//While B in Clover3, A in Clover3
		CoincidenceHit[2][3] = 104356;	//While B in Clover3, A in Clover4
		CoincidenceHit[3][0] = 103644;	//While B in Clover4, A in Clover1
		CoincidenceHit[3][1] = 99496;	//While B in Clover4, A in Clover2
		CoincidenceHit[3][2] = 114182;	//While B in Clover4, A in Clover3
		CoincidenceHit[3][3] = 99496;	//While B in Clover4, A in Clover4

		InitialHit[0][0] = 1.04039e+06;
		InitialHit[0][1] = 932421;
		InitialHit[0][2] = 1.14775e+06;
		InitialHit[0][3] = 1.01684e+06;
		InitialHit[1][0] = 1.04062e+06;
		InitialHit[1][1] = 932696;
		InitialHit[1][2] = 1.14724e+06;
		InitialHit[1][3] = 1.01729e+06;
	}

	//Case that you don't consider angular correlation of gammas. (like 22Na-->2 * 511 keV).
	else
	{
		for(int i=0; i<4; i++)
		{
			for(int j=0; j<4; j++)
			{
				CoincidenceHit[i][j] = 1;
			}
		}

		InitialHit[0][0] = 1.04039e+06;
		InitialHit[0][1] = 932421;
		InitialHit[0][2] = 1.14775e+06;
		InitialHit[0][3] = 1.01684e+06;
		InitialHit[1][0] = 1.04062e+06;
		InitialHit[1][1] = 932696;
		InitialHit[1][2] = 1.14724e+06;
		InitialHit[1][3] = 1.01729e+06;
	}

//////////////////////////////////////////////READ 60Co DATA File and Make Histograms///////////////////////////////////////////////////////////////////

	//60Co peaks and their branching ratio which will be used for efficiency measurements.
	int	NumOfPeak = 2;	//change this if you need

	//Note: The second gamma's efficiency is determined from the ratio to the first gamma. 
	//double	Epeaks[2] = {1173.228, 1332.492};	//in keV
	//double	Bratios[2] = {0.9986, 0.9998};
	//double	Epeaks[2] = {1332.492, 1173.228};	//in keV
	//double	Bratios[2] = {0.9998, 0.9986};
	double	Epeaks[2] = {511.0, 511.0};	//in keV
	double	Bratios[2] = {1.0, 1.0};

	//spectra values
	int lowerbound, upperbound, nobins;
    	lowerbound = 0; //in keV
    	upperbound = 3000; //in keV
    	nobins = 1000;

	//initiate list of Histograms
	TH1F* getamuAddBack_E[4];
	TH1F* getamuAddBack_E_Coincidence[4][4];
	TString nameTitle;

	for (int iClover =0; iClover<4 ; iClover++)
	{
		nameTitle =Form("GeTAMU_Clover%d_AddBack_E",iClover+1);
		getamuAddBack_E[iClover]= new TH1F (nameTitle,nameTitle,nobins,lowerbound,upperbound);

		for (int iClover2 =0; iClover2<4 ; iClover2++)
		{
			nameTitle =Form("GeTAMU_Clover%d_AddBack_E_CoincidenceWithClover%d",iClover+1, iClover2+1);
			getamuAddBack_E_Coincidence[iClover][iClover2]= new TH1F (nameTitle,nameTitle,nobins,lowerbound,upperbound);
		}

		
	}

	if(gSystem->AccessPathName(plotsFileName))
	{ //checks if the file exist already, condition is "true" if not
    		cout << "No file to calibrate found - creating one now..." << endl;

		//Open TFile
		TFile* nptDataFile = new TFile(Na22FileName.Data(),"OPEN");
	 	if (nptDataFile == 0)
		{
		  // if we cannot open the file, print an error message and return immediatly
		  printf("Error: cannot open this file: %s \n",Na22FileName.Data());
		  return;
	   	}

		nptDataFile->ls();

		//Load Tree
		TTree* tree = (TTree*) nptDataFile->Get("PhysicsTree");
		//Set Branch
		//tree->SetBranchStatus( "GeTAMU" , true )             ;
		tree->SetBranchAddress( "GeTAMU" , &getamuPhysics )     ;

		//Loop on tree and fill the histograms
		int entries = tree->GetEntries();
		cout << " INFO: Number of entries in tree: " << entries << endl;

		for(int i = 0 ; i < entries; i++)
		{
	  		if (i%(entries/100)) printf("\r treated %2.f percent ",100.0*i/entries);
	  		tree->GetEntry(i);

	  		unsigned int sizeAddBackE = getamuPhysics->AddBack_E.size();

	  		for(unsigned int j = 0 ; j < sizeAddBackE ; ++j)
			{
				unsigned short clover = getamuPhysics->AddBack_Clover[j];
				double energy = getamuPhysics->AddBack_E[j];

				if( clover>0 && energy>lowerbound && energy<upperbound)	getamuAddBack_E[clover-1]->Fill(energy);

				//Assuming coincidence of 1173 keV and 1332 keV gammas from 60Co in different Clovers.
				//if gamma energy is 1332 keV +- 10 keV
				if( clover>0 && energy>Epeaks[0]-10.0 && energy<Epeaks[0]+10.0)
				{
	  				for(unsigned int k = 0 ; k < sizeAddBackE ; ++k)
					{
						if(j==k)	continue;

						unsigned short clover_coincidence = getamuPhysics->AddBack_Clover[k];
						double energy_coincidence = getamuPhysics->AddBack_E[k];

						//if gamma energy is 1173 keV +- 10 keV
						if( clover_coincidence>0 && energy_coincidence>lowerbound && energy_coincidence<upperbound)	getamuAddBack_E_Coincidence[clover-1][clover_coincidence-1]->Fill(energy_coincidence);
					}
				}
	  		}
		}

		nptDataFile->Close();

		//Write all the data in the file
    		fileToCalibrate = new TFile(plotsFileName,"RECREATE");
		fileToCalibrate->cd();
		for (int iClover =0; iClover<4 ; iClover++)
		{
			getamuAddBack_E[iClover]->Write();

			for (int iClover2 =0; iClover2<4 ; iClover2++)	getamuAddBack_E_Coincidence[iClover][iClover2]->Write();
		}

		fileToCalibrate->Write();
  	}

	else 
	{
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		cout << " The file " << plotsFileName << " is found in the present directory, delete or change name to create a new calibration " << endl;
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		fileToCalibrate = new TFile(plotsFileName,"UPDATE");
		//return;
	}

	fileToCalibrate->cd();


//////////////////////////////////////////////FIT 60Co File /////////////////////////////////////////////////////////////////////////////////////////

	TH1F* currentHist;

   	char tgName[300];
        char sName[300];
        char ext[100];
        char imgName[100];
        char cName[100];

	//x: Gamma source energy, ex: error, y: efficiency, ey: error.
	double	x[4][4], ex[4][4], y[4][4], ey[4][4];
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			x[i][j] = j+1;
			ex[i][j] = 0.0;
			y[i][j] = -1.0;
			ey[i][j] = 0.0;
		}
	}

	double	SourceActivity[4] = {0.0, 0.0, 0.0, 0.0};	// This means the number of (e.g.) 1332 keV gamma without coincidence.
	double	BGActivityRatio[4] = {0.0, 0.0, 0.0, 0.0};	// This means the number of (e.g.) 1173 keV gamma to that of 1332 keV gamma. Needed to subtract "Chance coincidence (1332-1332, 1173-1173 coincidence which are impossible)"
	double	BGActivity;	//Chance coincidence (1332-1332keV, 1173-1173)

	for(unsigned int iClover =0; iClover<4 ; iClover++)
	{

		//NOTE!! In this loop, iClover2=-1 is the case without coincidence (singles), and iClover2=0,1,2,3 are coincidence with CL1,2,3,4.
		for(int iClover2 =-1; iClover2<4 ; iClover2++)
		{

	   		currentHist = NULL;

			if(iClover2<0)
			{
				nameTitle =Form("GeTAMU_Clover%d_AddBack_E",iClover+1);
				currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data());
				cout << "Number of entries (must be >100) is:  ["<<iClover << "] "<< currentHist->GetEntries() << endl; // used for debugging
			}

			else
			{
				nameTitle =Form("GeTAMU_Clover%d_AddBack_E_CoincidenceWithClover%d",iClover+1,iClover2+1);
				currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data());
				cout << "Number of entries (must be >100) is:  ["<<iClover << "] ["<<iClover2 << "] " << currentHist->GetEntries() << endl; // used for debugging
			}

			if (currentHist->GetEntries()>100)
			{

				for(unsigned int j=0; j<NumOfPeak; j++)
				{

					TString nameTitle2; 
					if(iClover2<0)	nameTitle2 =Form("GeTAMU_Clover%d_AddBack_E_%dkeV",iClover+1,(int)Epeaks[j]);
					else	nameTitle2 =Form("GeTAMU_Clover%d_AddBack_E_CoincidenceWithClover%d_%dkeV",iClover+1,iClover2+1,(int)Epeaks[j]);
					currentHist->SetName(nameTitle2+"_fit");
					currentHist->SetTitle(nameTitle2+"_fit");

					//Fitting function
					TF1	*f1;
					//For drawing the fit in the histogram
					TF1	*fga1;
					//For drawing the background in the histogram
					TF1	*fga2;

					f1 = new TF1("gausfit", gausfit, Epeaks[j]-15, Epeaks[j]+15, 6);
					fga1 = new TF1("fga1", "[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]) )+([3]+[4]*(x-[1])+[5]*(x-[1])*(x-[1]))", Epeaks[j]-15, Epeaks[j]+15);
					//background function (a+b*x+c*x*x)
					fga2 = new TF1("fga2", "[1]+[2]*(x-[0])+[3]*(x-[0])*(x-[0])", Epeaks[j]-15, Epeaks[j]+15);

	      				double Height = currentHist->GetBinContent(currentHist->FindBin(Epeaks[j]));
					double	BG_low = currentHist->GetBinContent(currentHist->FindBin(Epeaks[j]-30.0));
					double	BG_high = currentHist->GetBinContent(currentHist->FindBin(Epeaks[j]+30.0));

					//case the low energy gamma nearby the threshold (that case, no counts)
					if(BG_low == 0)	BG_low = BG_high;
					double	BG_ave = (BG_low + BG_high)/2.0;

					if(Height>BG_ave)	f1->SetParameters((double)(Height-BG_ave), (double)Epeaks[j], 3.0, BG_ave, 0.0, 0.0);
					else	f1->SetParameters(1, (double)Epeaks[j], 3.0, 0.0, 0.0, 0.0);

					//Fitting conditions. You can change these depending on your data.
					f1->FixParameter(0, (double)Height);
					//f1->FixParameter(1, Hgamma[s][j]+0.5);
					//f1->SetParLimits(1, Hge[s][t][j]-3.0,Hge[s][t][j]+3.0);
					f1->SetParLimits(1, Epeaks[j]-10.0,Epeaks[j]+10.0);
					f1->SetParLimits(2, 1.0, 10.0);
					//f1->FixParameter(2, 2.0);
					//f1->FixParameter(3, ((double)bg_low+(double)bg_high)/2.0);

					currentHist->Fit("gausfit", "Q", "", Epeaks[j]-15, Epeaks[j]+15);
					currentHist->Draw();
       		         		currentHist->Write(currentHist->GetName(),TObject::kOverwrite);

					//amplitude
					double ga1 = f1->GetParameter(0);
					//energy
					double ga2 = f1->GetParameter(1);
					//sigma
					double ga3 = f1->GetParameter(2);
					//bg param1
					double ga4 = f1->GetParameter(3);
					//bg param2
					double ga5 = f1->GetParameter(4);
					//bg param3
					double ga6 = f1->GetParameter(5);

					double ga7 = f1->GetParError(1);

	
					fga1->SetParameter(0, ga1);
					fga1->SetParameter(1, ga2);
					fga1->SetParameter(2, ga3);
					fga1->SetParameter(3, ga4);
					fga1->SetParameter(4, ga5);
					fga1->SetParameter(5, ga6);
					fga1->SetLineColor(4);
					double totalcounts = fga1-> Integral(Epeaks[j]-15, Epeaks[j]+15);

					fga2->SetParameter(0, ga2);
					fga2->SetParameter(1, ga4);
					fga2->SetParameter(2, ga5);
					fga2->SetParameter(3, ga6); double totalcounts_bg = fga2-> Integral(Epeaks[j]-15, Epeaks[j]+15);

					totalcounts -= totalcounts_bg;
					//error estimated by sqrt(x^2+y^2)
					double	totalcounts_err = pow(totalcounts+totalcounts_bg,0.5);

					totalcounts /= Bratios[j];
					totalcounts_err /= Bratios[j];

					//Singles gammas (i.e., No coincidence)
					if(iClover2==-1)
					{
						if(j==0)	SourceActivity[iClover] = totalcounts;
						//else if(j==1)	BGActivityRatio[iClover] = totalcounts / SourceActivity[iClover];	
						//cout << SourceActivity[iClover] << " " << BGActivityRatio[iClover] << endl;
					}

					//Coincidence gammas
					else if(iClover2>=0)
					{
						//if(j==0)	BGActivity = totalcounts * BGActivityRatio[iClover];
						if(j==0)	BGActivity = 0;
						else if(j==1)
						{
							x[iClover][iClover2]=iClover2+1;

							ex[iClover][iClover2]=0.0;	//error was set 0 because gamma source's energy is accurate enough. 

							y[iClover][iClover2]=(totalcounts - BGActivity)/SourceActivity[iClover];
							ey[iClover][iClover2]=totalcounts_err/SourceActivity[iClover];

							//Angular correlation must be corrected for. When we detect a gamma in CL1, we need to know how much gammas can be expected in CL2 (definitely <100%).
							//By correcting the angular correlation, we can get an intrinsic efficiency (i.e., assuming all gammas hit the Clover (but not detected))
							AngularCor[iClover][iClover2] = 1.0;
							
							//By correcting for the geometrical factor (angular coverage), you can convert the intrinsic efficiency to absolute efficiency (which we need for data correction).
							GeometricalFactor[iClover][iClover2] = InitialHit[1][iClover2] / SourceParticles;	//SourceParticles = 1e7 simulated gamma source

							//So this is the absolute efficiency and error.
							y[iClover][iClover2] = y[iClover][iClover2] / AngularCor[iClover][iClover2] * GeometricalFactor[iClover][iClover2];
							ey[iClover][iClover2] = ey[iClover][iClover2] / AngularCor[iClover][iClover2] * GeometricalFactor[iClover][iClover2];
	
							//cout << AngularCor[iClover][iClover2] << " " << GeometricalFactor[iClover][iClover2] << endl;
							//cout << y[iClover][iClover2] << endl;

						}
					}

				}	//end of 60Co sources loop
			}
		}	//iClover2 loop
	}	//iClover loop


//////////////////////////////////////////////MAKE A GRAPH OF 60Co EFFICIENCY AS A FUNCTION OF CLOVER # ///////////////////////////////////////////////////////////////////

	TCanvas *cD2;
  	TGraphErrors *grCal[4];
	for(unsigned int iClover =0; iClover<4 ; iClover++)
	{

		if(iClover==0)
		{
			sprintf(cName,"Ge_CoincidenceEfficiency_22Na");
			cD2 = new TCanvas(cName,cName, 200,40,1150, 1150);
			cD2->GetFrame()->SetFillColor(21);
			cD2->GetFrame()->SetBorderMode(1);
			cD2->GetFrame()->SetBorderSize(200); 
		}


		grCal[iClover] = new TGraphErrors(4,x[iClover], y[iClover], ex[iClover], ey[iClover]);
		sprintf(tgName,"Ge_CoincidenceEfficiency_22Na");
		//I didn't put legends, but you can do it
		//grCal[iClover]->GetXaxis()->SetTitle("Clover Number");
                grCal[iClover]->SetTitle(tgName);
		grCal[iClover]->GetXaxis()->SetTitle("Clover # (pink: by CL1, red: CL2, blue: CL3, green: CL4, black: average)");
		sprintf(tgName,"Absolute Efficiency of %dkeV",(int)Epeaks[1]);
		grCal[iClover]->GetYaxis()->SetTitle(tgName);
                grCal[iClover]->SetMarkerStyle(21);
               	grCal[iClover]->SetLineWidth(2);

		if(iClover==0)
		{
                	grCal[iClover]->SetMarkerColor(6);
                	grCal[iClover]->SetLineColor(6);
	  		grCal[iClover]->Draw("AP");

			grCal[iClover]->GetXaxis()->SetLimits(0.,5.);
			grCal[iClover]->GetYaxis()->SetRangeUser(0.,0.05);
		//	cD2->Update();
		}
		else if(iClover==1)	
		{
                	grCal[iClover]->SetLineColor(2);
                	grCal[iClover]->SetMarkerColor(2);
			grCal[iClover]->Draw("P & same");
		//	cD2->Update();
		}
		else if(iClover==2)	
		{
                	grCal[iClover]->SetLineColor(3);
                	grCal[iClover]->SetMarkerColor(3);
			grCal[iClover]->Draw("P & same");
		//	cD2->Update();
		}

		else if(iClover==3)
		{
                	grCal[iClover]->SetLineColor(4);
                	grCal[iClover]->SetMarkerColor(4);
			grCal[iClover]->Draw("P & same");

			sprintf(imgName,"Ge_CoincidenceEfficiency_22Na");
	        	sprintf(ext,"png");
   	        	//sprintf(sName,"/home/shuyaota/T40/GeTAMU/Efficiency/Graphs/%s.%s",imgName,ext);
   	        	sprintf(sName,"~/nptool/Outputs/%s.%s",imgName,ext);
			cout << sName << endl;
      			cD2->SaveAs(sName);
		}
	}

	fileToCalibrate->Close();


//////////////////////////////////////////////WRITE EFFICIENCY TO FILE ///////////////////////////////////////////////////////////////////

	//write in a text file. 
  	TString EfficiencyfName( Na22FileName( Na22FileName.Last('/')+1, Na22FileName.Length() ) );
  	EfficiencyfName.ReplaceAll("root","txt");
  	EfficiencyfName = "./GeTAMU_Efficiency_"+EfficiencyfName;
  	ofstream myfile;
  	myfile.open (EfficiencyfName.Data());

	//Use this output as input to EfficiencyGeTamu.C
	for(int iClover2 = 0 ; iClover2 < 4 ; iClover2++)
	{
		myfile << iClover2+1<<" " ;
		myfile << Epeaks[1] <<" " ;
		for(int iClover = 0 ; iClover < 4 ; iClover++)	
		{
			if(abs(iClover-iClover2) == 2)
			{
				myfile << y[iClover][iClover2]<<" " ;
				myfile << ey[iClover][iClover2]<<" " ;
			}
		}
		myfile<<endl;
	}
	myfile.close();

}
