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
	This code is the second (2/2) step to get the absolute efficiency cuver for GeTAMU.
	Here we normalize the efficiency from 152Eu to the absolute efficiency 
	obtained from 60Co and (if you have) 22Na data.

	Note you need to have Natural Background spectrum data too. 
	You can use a dammy spectrum, however, (ex., set the same name as 152Eu file),
	and then set BG_coeff = 0. By doing this, contribution from natural background is ignored. 

	Also NOTE! Root files are all .root after npanalysis (not midas2nptool).
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void EfficiencyGeTamu(TString Eu152FileName="/home/shuyaota/nptool/Outputs/Analysis/EXPT6_ER72_0.root", TString BGFileName="/home/shuyaota/nptool/Outputs/Analysis/EXPT6_ER77_0-ER77_1.root", TString pathToNormalizationFile="GeTAMU_Efficiency_EXPT6_ER73_0.txt")
//IF YOU HAVE 22Na Efficiency output too, use this line instead of above. That case, remember remove the comment out around line 90 about pathToNormalizationFile2.
//void EfficiencyGeTamu(TString Eu152FileName="/home/shuyaota/nptool/Outputs/Analysis/EXPT6_ER72_0.root", TString BGFileName="/home/shuyaota/nptool/Outputs/Analysis/EXPT6_ER77_0-ER77_1.root", TString pathToNormalizationFile="GeTAMU_Efficiency_EXPT6_ER73_0.txt", TString pathToNormalizationFile2="GeTAMU_Efficiency_EXPT6_ER74_0.txt")
{

	//SET how much you need to subtract natural BG spectrum from the current spectrum. 152Eu spectrum - BG_coeff[clover num] * BG
	//You can decide, e.g. by normalizing 1764 keV gamma (214Bi decay) or 2615 keV (208Tl) in BG to 152Eu. Note 1480 keV (40K) is overlapped by 152Eu 1476 keV gamma.
	//Try first with random numbers and then see the created spectrum. You'll see if you over/underestimated.
	double	BG_coeff[4] = {0.2,0.2,0.2,0.2};

	//generate the outputFileName
	TString plotsFileName( Eu152FileName( Eu152FileName.Last('/')+1, Eu152FileName.Length() ) );
	plotsFileName= "inspect_"+plotsFileName;

	TString plotsFileName2( BGFileName( BGFileName.Last('/')+1, BGFileName.Length() ) );
	plotsFileName2= "inspect_"+plotsFileName2;

	//initiate output variables
	vector < TString > badchannels; //channels with bad data

	//Nptool data
	TGeTAMUPhysics* getamuPhysics = new TGeTAMUPhysics() ;
	TFile* fileToCalibrate;
	TFile* fileToCalibrate_BG;

	//vector for Efficiency output
	vector < vector<double> > coeff;
	coeff.resize(4);

//////////////////////////////////////////////READ NORMALIZATION File///////////////////////////////////////////////////////////////////
	ifstream inputFile;
	inputFile.open(pathToNormalizationFile);
	int cloverNumber;
	double eng1, eng2, eff1, eff1_err, eff2, eff2_err;
	while (!inputFile.eof()) 
	{
		inputFile >> cloverNumber >> eng1 >> eff1 >> eff1_err >> eng2 >> eff2 >> eff2_err;
	}
	inputFile.close();

	double eng3, eff3, eff3_err;
	eng3 = 0.0;
	eff3 = 0.0;
	eff3_err = 0.0;

	//Delete the comment out if you have 22Na efficiency outputs too.
	/*
	inputFile.open(pathToNormalizationFile2);
	while (!inputFile.eof()) 
	{
		inputFile >> cloverNumber >> eng3 >> eff3 >> eff3_err;
	}
	inputFile.close();
	*/

//////////////////////////////////////////////READ Natural Background File///////////////////////////////////////////////////////////////////
	//spectra values
	int lowerbound, upperbound, nobins;
    	lowerbound = 0; //in keV
    	upperbound = 3000; //in keV
    	nobins = 1000;

	//initiate list of Histograms
	TH1F* getamuAddBack_E_BG[4];
	TString nameTitle_BG;

	for (int iClover =0; iClover<4 ; iClover++)
	{
		nameTitle_BG =Form("GeTAMU_Clover%d_AddBack_E_BG",iClover+1);
		getamuAddBack_E_BG[iClover]= new TH1F (nameTitle_BG,nameTitle_BG,nobins,lowerbound,upperbound);
	}

	if(gSystem->AccessPathName(plotsFileName2))
	{ //checks if the file exist already, condition is "true" if not
    		cout << "No file to calibrate found - creating one now..." << endl;

		//Open TFile
		TFile* nptDataFile_BG = new TFile(BGFileName.Data(),"OPEN");
	 	if (nptDataFile_BG == 0)
		{
		  // if we cannot open the file, print an error message and return immediatly
		  printf("Error: cannot open this file: %s \n",BGFileName.Data());
		  return;
	   	}

		nptDataFile_BG->ls();

		//Load Tree
		TTree* tree_BG = (TTree*) nptDataFile_BG->Get("PhysicsTree");
		//Set Branch
		//tree->SetBranchStatus( "GeTAMU" , true )             ;
		tree_BG->SetBranchAddress( "GeTAMU" , &getamuPhysics )     ;

		//Loop on tree and fill the histograms
		int entries_BG = tree_BG->GetEntries();
		cout << " INFO: Number of entries in BG tree: " << entries_BG << endl;

		for(int i = 0 ; i < entries_BG; i++)
		{
	  		if (i%(entries_BG/100)) printf("\r treated %2.f percent ",100.0*i/entries_BG);
	  		tree_BG->GetEntry(i);

	  		unsigned int sizeAddBackE_BG = getamuPhysics->AddBack_E.size();

	  		for(unsigned int j = 0 ; j < sizeAddBackE_BG ; ++j)
			{
				unsigned short clover = getamuPhysics->AddBack_Clover[j];
				double energy = getamuPhysics->AddBack_E[j];

			if( clover>0 && energy>lowerbound && energy<upperbound)	getamuAddBack_E_BG[clover-1]->Fill(energy);
	  		}
		}

		nptDataFile_BG->Close();

		//Write all the data in the file
    		fileToCalibrate_BG = new TFile(plotsFileName2,"RECREATE");
		fileToCalibrate_BG->cd();
		for (int iClover =0; iClover<4 ; iClover++)	getamuAddBack_E_BG[iClover]->Write();

		fileToCalibrate_BG->Write();
  	}

	else 
	{
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		cout << " The file " << plotsFileName2 << " is found in the present directory, delete or change name to create a new calibration " << endl;
		cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "<<endl;
		fileToCalibrate_BG = new TFile(plotsFileName2,"UPDATE");
		//return;
	}

	fileToCalibrate_BG->cd();
	//fileToCalibrate_BG->Close();


//////////////////////////////////////////////READ 152EU File///////////////////////////////////////////////////////////////////
	//spectra values
    	lowerbound = 0; //in keV
    	upperbound = 3000; //in keV
    	nobins = 1000;

	//initiate list of Histograms
	TH1F* getamuAddBack_E[4];
	TString nameTitle;

	for (int iClover =0; iClover<4 ; iClover++)
	{
		nameTitle =Form("GeTAMU_Clover%d_AddBack_E",iClover+1);
		getamuAddBack_E[iClover]= new TH1F (nameTitle,nameTitle,nobins,lowerbound,upperbound);
	}

	if(gSystem->AccessPathName(plotsFileName))
	{ //checks if the file exist already, condition is "true" if not
    		cout << "No file to calibrate found - creating one now using triple alpha spectra..." << endl;

		//Open TFile
		TFile* nptDataFile = new TFile(Eu152FileName.Data(),"OPEN");
	 	if (nptDataFile == 0)
		{
		  // if we cannot open the file, print an error message and return immediatly
		  printf("Error: cannot open this file: %s \n",Eu152FileName.Data());
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
	  		}
		}

		nptDataFile->Close();

		//Write all the data in the file
    		fileToCalibrate = new TFile(plotsFileName,"RECREATE");
		fileToCalibrate->cd();
		for (int iClover =0; iClover<4 ; iClover++)	getamuAddBack_E[iClover]->Write();

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

	//fileToCalibrate->cd();


//////////////////////////////////////////////FIT 152EU File /////////////////////////////////////////////////////////////////////////////////////////

	TH1F* currentHist;
	TH1F* currentHist_BG;

	for(unsigned int iClover =0; iClover<4 ; iClover++)
	{
		//152Eu peaks and their branching ratio which will be used for efficiency measurements.
		int	NumOfPeak = 11;	//change this if you need
		double	Epeaks[11] = {121.78, 244.70, 344.28, 411.11, 443.97, 778.90, 867.37, 964.08, 1085.9, 1112.1, 1408.0};	//in keV
		double	Bratios[11] = {0.2837, 0.0753, 0.2657, 0.02238,	0.03125, 0.1297,0.04214,0.1463,	0.1013,	0.1354,	0.21};

   		currentHist = NULL;

		nameTitle =Form("GeTAMU_Clover%d_AddBack_E",iClover+1);
		currentHist = (TH1F*) fileToCalibrate->Get(nameTitle.Data());
		cout << "Number of entries (must be >100) is:  ["<<iClover << "] "<< currentHist->GetEntries() << endl; // used for debugging

		nameTitle_BG =Form("GeTAMU_Clover%d_AddBack_E_BG",iClover+1);
		currentHist_BG = (TH1F*) fileToCalibrate_BG->Get(nameTitle_BG.Data());
		cout << "Number of entries (must be >100) is:  ["<<iClover << "] "<< currentHist_BG->GetEntries() << endl; // used for debugging

		currentHist->Add(currentHist_BG,-BG_coeff[iClover]);

		if (currentHist->GetEntries()>100)
		{

			double	x[20],y[20],ex[20],ey[20];
        		char cName[100];

			for(unsigned int j=0; j<NumOfPeak; j++)
			{
				TString nameTitle2; 
				nameTitle2 =Form("GeTAMU_Clover%d_AddBack_E_%dkeV",iClover+1,(int)Epeaks[j]);
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

				//f1->SetParLimits(0, (double)gesumenergy[s][Hge[s][t][j]]-bg_ave-100.0, (double)gesumenergy[s][Hge[s][t][j]]-bg_ave+100);


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
				fga2->SetParameter(3, ga6);
				double totalcounts_bg = fga2-> Integral(Epeaks[j]-15, Epeaks[j]+15);

				totalcounts -= totalcounts_bg;
				double	totalcounts_err = pow(totalcounts,0.5);
				//cout << ga1 << " " << ga2 << " " << ga3 << " " << ga5 << " " << ga6 << " " << ga7 << " " << totalcounts << " " << totalcounts_bg << endl;

				totalcounts /= Bratios[j];
				totalcounts_err /= Bratios[j];

				//tentatively. relative efficiency. SourceActivity doesn't mean anything AT ALL but normalizing the efficiency and put into the range of the graph.
				double	SourceActivity = 370000.0 * 5;	//37 kBq * 5sec.
				x[j]=Epeaks[j];
				ex[j]=0.0;
				y[j]=totalcounts/SourceActivity;
				ey[j]=totalcounts_err/SourceActivity;

				cout << x[j] << " " << y[j] << " " << ey[j] << " " << totalcounts << endl;
			}
			//end of 152Eu sources loop


//////////////////////////////////////////////GET 152EU EFFICIENCY ///////////////////////////////////////////////////////////////////

   		     	char tgName[100];
        		char sName[300];
        		char ext[100];
        		char imgName[100];
    			char iPath[300] = "/home/shuyaota/T40/GeTAMU/Efficiency/Graphs/";


			TF1	*lFit, *lFit2;
			TCanvas *cD2;
  			TGraphErrors *grCal;

			double	fVals[3];

			sprintf(cName,"Ge_RelativeEfficiency_CL%d", iClover+1);
			cD2 = new TCanvas(cName,cName, 200,40,1150, 1150);
			cD2->GetFrame()->SetFillColor(21);
			cD2->GetFrame()->SetBorderMode(1);
			cD2->GetFrame()->SetBorderSize(200); 

			//cD2->SetLogx();
			//cD2->SetLogy();

			grCal = new TGraphErrors(NumOfPeak,x, y, ex, ey);
			sprintf(tgName,"Ge_RelativeEfficiency_cl%d",iClover+1);
			grCal->GetXaxis()->SetTitle("Gamma Energy (keV)");
			grCal->GetYaxis()->SetTitle("Relative Efficiency");

			sprintf(cName,"lFit_cl%d",iClover);

                	lFit = new TF1(cName,"exp([0] + [1] * log(x/1000.0) + [2] * pow(log(x/1000),2.0))",100,10000.);

                	lFit->SetParameter(0,1.0);
                	lFit->SetParameter(1,1.0);
                	lFit->SetParameter(2,1.0);
                	//grCal->Fit(cName,"RNQ");
                	grCal->Fit(cName,"Q","", 100.0,3000.0);
	  		grCal->Draw();
                	lFit->GetParameters(&fVals[0]);
                	//cout << fVals[0] << " " << fVals[1] << " " << fVals[2] << endl;

                	grCal->SetTitle(tgName);
                	grCal->SetLineColor(1);
                	grCal->SetLineWidth(2);
                	grCal->SetMarkerColor(1);
                	grCal->SetMarkerStyle(21);

	  		grCal->Draw("AP");
	  		lFit->Draw("same");

			grCal->GetXaxis()->SetLimits(0.,3000.);
			grCal->GetYaxis()->SetRangeUser(0.,1.);

			cD2->Update();

			sprintf(imgName,"Ge_RelativeEfficiency_cl%d",iClover+1);
	        	sprintf(ext,"png");
   	        	//sprintf(sName,"/home/shuyaota/T40/GeTAMU/Efficiency/Graphs/%s.%s",imgName,ext);
   	        	sprintf(sName,"~/nptool/Outputs/%s.%s",imgName,ext);
			cout << sName << endl;
      			cD2->SaveAs(sName);

//////////////////////////////////////////////NORMALIZE TO 60Co EFFICIENCY///////////////////////////////////////////////////////////////////

			double	NormEnergy = eng2;	//60Co
			double	NormEfficiency = eff2;	//Efficiency@NormEnergy
			double	NormFactor = NormEfficiency / exp(fVals[0] + fVals[1] * log(NormEnergy/1000.0) + fVals[2] * pow(log(NormEnergy/1000),2.0));

			for(unsigned int j=0;j<NumOfPeak;j++)
			{
				y[j] *= NormFactor;
				ey[j] *= NormFactor;
			}

			x[NumOfPeak] = eng1;
			ex[NumOfPeak] = 0.0;
			y[NumOfPeak] = eff1;
			ey[NumOfPeak] = eff1_err;
			x[NumOfPeak+1] = eng2;
			ex[NumOfPeak+1] = 0.0;
			y[NumOfPeak+1] = eff2;
			ey[NumOfPeak+1] = eff2_err;
			NumOfPeak +=2;

			if(eng3>0)
			{
				x[NumOfPeak] = eng3;
				ex[NumOfPeak] = 0.0;
				y[NumOfPeak] = eff3;
				ey[NumOfPeak] = eff3_err;
				NumOfPeak ++;
			}

			TCanvas *cD3;
  			TGraphErrors *grCal2;

			sprintf(cName,"Ge_FinalEfficiency_cl%d_normalized", iClover+1);
			cD3 = new TCanvas(cName,cName, 200,40,1150, 1150);
			cD3->GetFrame()->SetFillColor(21);
			cD3->GetFrame()->SetBorderMode(1);
			cD3->GetFrame()->SetBorderSize(200); 

			grCal2 = new TGraphErrors(NumOfPeak,x, y, ex, ey);
			sprintf(tgName,"Ge_FinalEfficiency_cl%d",iClover+1);
			grCal2->GetXaxis()->SetTitle("Gamma Energy (keV)");
			grCal2->GetYaxis()->SetTitle("Absolute Efficiency");

			sprintf(cName,"lFit2_cl%d",iClover);


                	lFit2 = new TF1(cName,"exp([0] + [1] * log(x/1000.0) + [2] * pow(log(x/1000),2.0))",100,10000.);

	                lFit2->SetParameter(0,fVals[0]);
       	       		lFit2->SetParameter(1,fVals[1]);
       		        lFit2->SetParameter(2,fVals[2]);
                	grCal2->Fit(cName,"Q","", 100.0,3000.0);
		  	grCal2->Draw();
       		        lFit2->GetParameters(&fVals[0]);
                	//cout << fVals[0] << " " << fVals[1] << " " << fVals[2] << endl;
			for(unsigned int j=0; j<3; j++)	coeff[iClover].push_back(fVals[j]);

             		grCal2->SetTitle(tgName);
                	grCal2->SetLineColor(1);
                	grCal2->SetLineWidth(2);
                	grCal2->SetMarkerColor(1);
                	grCal2->SetMarkerStyle(21);

		  	grCal2->Draw("AP");
		  	lFit2->Draw("same");

			grCal2->GetXaxis()->SetLimits(0.,3000.);
			grCal2->GetYaxis()->SetRangeUser(0.,0.1);

			cD3->Update();

			sprintf(imgName,"Ge_FinalEfficiency_cl%d_Normalized",iClover+1);
	   		sprintf(ext,"png");
   	        	//sprintf(sName,"/home/shuyaota/T40/GeTAMU/Efficiency/Graphs/%s.%s",imgName,ext);
   	        	sprintf(sName,"~/nptool/Outputs/%s.%s",imgName,ext);
			//cout << sName << endl;
      			cD3->SaveAs(sName);

			//fprintf(fout, "%i %lf %lf\n", s, fVals[1],fVals[0]);

		}
	}

	fileToCalibrate->Close();


	//write in a text file
  	TString EfficiencyfName( Eu152FileName( Eu152FileName.Last('/')+1, Eu152FileName.Length() ) );
  	EfficiencyfName.ReplaceAll("root","txt");
  	EfficiencyfName = "./GeTAMU_Efficiency_"+EfficiencyfName;
  	ofstream myfile;
  	myfile.open (EfficiencyfName.Data());
	for(unsigned int iClover = 0 ; iClover < 4 ; iClover++)
	{
		myfile << iClover+1<<" " ;
		for(unsigned int j = 0 ; j < coeff.at(iClover).size() ; j++)	myfile << coeff.at(iClover).at(j)<<" " ;
		myfile<<endl;
	}
	myfile.close();

}
