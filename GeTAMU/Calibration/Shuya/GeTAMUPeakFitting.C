#include "TROOT.h"
#include "TMath.h"
#include "TTree.h"
#include "TEventList.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
//#include "/home/ota2/TMVA-v4.2.0/inc/TMVA/Timer.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include "TLeaf.h"
#include "TSystem.h"
//by Shuya 160922
#include "TApplication.h"
#include "TInterpreter.h"
#include "TSystem.h"
// NPTOOL
//#include "/home/t40/packages/nptool/NPLib/include/TTiaraBarrelData.h"
//#include "/home/t40/packages/nptool/NPLib/include/TTiaraHyballData.h"
//#include "/home/t40/packages/nptool/NPLib/include/TGeTAMUData.h"
#include "/home/shuyaota/nptool/NPLib/include/TTiaraBarrelData.h"
#include "/home/shuyaota/nptool/NPLib/include/TTiaraHyballData.h"
#include "/home/shuyaota/nptool/NPLib/include/TGeTAMUData.h"
//by Shuya 160923
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include "TStyle.h"
//by Shuya 161102
#include "TApplication.h"
//#include "/home/t40/packages/nptool/NPLib/include/TMicromegaData.h"


///////////////////////////////NOTE!! YOU HAVE TO CHANGE PARAMETERS DEPENDING ON GAMMA SOURCES///////////////////////////////////
//Int_t EventNum
//Ege[]
//Input Data FileName
//lFit (range of fitting)
//lFit (range of drawing (if any))


////Search "Change" will take you to these variables.

/////////////////////////////////////////////////////////////////////////////////////////////////////////

using std::vector;

TGeTAMUData*      GeTAMU = new TGeTAMUData();



//Skewed Gaussian based on Radware fitting (NOTE! This is just a polinomial background and no Skewing function is included yet).
double skewedgausfit(double *x, double *par)
{
float xx = x[0];
double f = par[0]*exp(-0.5*( (xx-par[1])*(xx-par[1]) ) / (par[2]*par[2]) );
f += par[3] + par[4]*(xx-par[0]);

return f;
}

double gausfit(double *x, double *par)
{
float xx = x[0];
double f = par[0]*exp(-0.5*( (xx-par[1])*(xx-par[1]) ) / (par[2]*par[2]) );
f += par[3];
return f;
}


int midas2nptool_Ge_PeakFitting()
{
FILE *fin;
FILE *fout;
FILE *fout2;

int debug = 0;
int i,j,k,q;
int range = 0;
int maxrange = 0;
int rangeproj[1000] = {0};
int deuteronsum[20][3000] = {{0}};
int deuterontotal[3000] = {0};
int protonsum[3000] = {0};
float theta = 0;
short int array[3000] = {0};

//I just added to fixe an error (TBuffa...) by SHuya 160921.
if (!(gInterpreter->IsLoaded("vector")))
  gInterpreter->ProcessLine("#include <vector>");

//gSystem->Load("libdict.so");
gSystem->Load("libTree.so");


//Ge
vector<vector<vector<int> > >gesegmentspec;
gesegmentspec.resize(4);
for(i=0;i<4;i++)
{
	gesegmentspec[i].resize(4);
	for(j=0;j<4;j++)	gesegmentspec[i][j].resize(5000);
}

vector<vector<vector<int> > >gecrystalspec;
gecrystalspec.resize(4);
for(i=0;i<4;i++)
{
	gecrystalspec[i].resize(4);
	for(j=0;j<4;j++)	gecrystalspec[i][j].resize(5000);
}

printf("Vector creation complete\n");

int fileloop = 0;
int filenum = 0;
int ncols = 0;
int strlength = 0;
FILE *batch;
char temp[400] = "";
char subdir[400] = "";
char finname[400] = "";
char finname0[400] = "";
char finname1[400] = "";
int	bg_low=0;
int	bg_high=0;
double	bg_ave=0;

Double_t	x[20], y[20], z[20], ex[20], ey[20], ez[20];

//////////////////////////////////////////////////////////////////////
//Change depending on files PART 1 OUT OF 3!!!!
//These are the number of peaks you are expecting for calibration.

Int_t EventNum;

//137Cs
//EventNum = 1;

//60Co
EventNum = 2;

//152Eu
//EventNum = 4;
//////////////////////////////////////////////////////////////////////


vector<double> ege;
ege.resize(EventNum);

vector<double> Ege;
Ege.resize(EventNum);

//////////////////////////////////////////////////////////////////////
//Change depending on files PART 2 OUT OF 3!!!!
//Gamma energies.

//137Cs
//Ege[0] = 661.657;

//60Co
Ege[0] = 1173.228;
Ege[1] = 1332.492;

//152Eu
//Ege[0] = 121.8;
//Ege[1] = 244.7;
//Ege[2] = 344.3;	
//Ege[3] = 1408.0;
//////////////////////////////////////////////////////////////////////


vector<vector<vector <int> > > Hge;
Hge.resize(4);
for(i=0;i<4;i++)
{
	Hge[i].resize(7);
	for(j=0;j<7;j++)
	{
		Hge[i][j].resize(EventNum);
	}
}

for(i=0;i<4;i++)
{
	for(j=0;j<7;j++)
	{
		for(int k=0;k<EventNum;k++)	Hge[i][j][k]=0;
	}
}



//////////////////////////////////////////////////////////////////////
//Change depending on files PART 3 OUT OF 3!!!!
//Data file (root file converted from MIDAS by midas2nptool) and its path.

char fName[500];
char rPath[500];

//EXPT4's 137Cs run.
//sprintf(fName,"ER189_0");
//EXPT4's 60Co run.
sprintf(fName,"ER190_0");
//EXPT4's 152Eu run.
//sprintf(fName,"ER191_0");
sprintf(rPath,"/home/shuyaota/midas2nptool/root/EXPT4/%s.root",fName);

//These are Output directories for Output files and figures. Create these directories for your convenience.
//char iPath[300] = "/home/shuyaota/midas2nptool/Graphs/EXPT4";
//char iPath2[300] = "/home/shuyaota/midas2nptool/Outputs/EXPT4";
char iPath[300] = "/home/shuyaota/midas2nptool/Graphs/TEST";
char iPath2[300] = "/home/shuyaota/midas2nptool/Outputs/TEST";
//////////////////////////////////////////////////////////////////////



char sName[500];
char ext[10];
char imgName[800];
char iName2[300];

//This is the final output you get (i.e., ADC value vs Energy)
sprintf(iName2,"%s/%s_calparam.dat",iPath2,fName);
fout = fopen(iName2,"w");
printf("Opened the output data file\n");



//////////////////////////////////////////////////////////////////////
//Reading the Input root file

TFile *finR = new TFile(rPath,"OPEN");
if(finR->IsOpen())	printf("root file opened for reading.\n");
else
{
   	printf("root file didn't open.\n");
	exit (1);
}


	TTree *tin = (TTree*)finR->Get("T40Tree");
	tin->SetBranchAddress("GeTAMU",&GeTAMU);


	for(int k = 0;k< tin->GetEntries();++k)
	{
		//cout << k << endl;
		tin->GetEntry(k);
	
		//Ge Segments
		for(int p=0;p<(int)(GeTAMU->GetMultiplicitySegmentE());p++)
		{
			short	geclover = GeTAMU->GetSegmentCloverNbrE(p);
			//cout << geclover << endl;
			geclover--;
			short	gesegment = GeTAMU->GetSegmentSegmentNbrE(p);
			gesegment--;
			//cout << geclover << " " << gesegment << endl;
			//cout << GeTAMU->GetGeSegmentEnergy(p) << endl;
			int	energy = (int)(GeTAMU->GetSegmentEnergy(p));

			//cout << geclover << " " << gesegment << " " << energy << endl;
    			if(geclover >= 0 && geclover < 4 && gesegment >= 0 && gesegment < 3 && energy >= 0 && energy < 5000)	gesegmentspec[geclover][gesegment][energy]++; // (DetNbr, StripNbr, Energy)

		}

		//Ge Crystals
		for(int p=0;p<(int)(GeTAMU->GetMultiplicityCoreE());p++)
		{
			short	geclover = GeTAMU->GetCoreCloverNbrE(p);
			geclover--;
			short	gecrystal = GeTAMU->GetCoreCrystalNbrE(p);
			gecrystal--;
			int	energy = (int)(GeTAMU->GetCoreEnergy(p));

			//cout << geclover << " " << gecrystal << " " << energy << endl;
    			if(geclover >= 0 && geclover < 4 && gecrystal >= 0 && gecrystal < 4 && energy >= 0 && energy < 5000)	gecrystalspec[geclover][gecrystal][energy]++; // (DetNbr, StripNbr, Energy)

		}
	}

finR->Close();
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
//Finding the peaks and fitting

TH1I	*h9[4][10];
TH1I	*h10[4][10];
//sigma means threshold. By Shuya 160923.
 double sigma = 0.2;
  Int_t nfound[4][10];
  const int nPeaks = EventNum;
char	brN[400];
  char cName[500];
  TCanvas *cD0;

//4 clovers, 4 crystals, 20 (potential gamma sources)
TH1I	*h11[4][10][20];
  TCanvas *cD[300];


//s==0; clover1, 1;clover2, 2;clover3, 3;clover4
for(int s=0;s<4;s++)
{

//if you want to skip some clovers, remove the comment out.
//if(s==0)	continue;


//These are options so you can analyze the clovers that you just need.
if(s<4)
//if(s==3)
{

	//3 segments (t=0-2), 4 cores (t=3-6)
	//if you want to analyze only 1 core (or segment), remove the comment out 
	//for(int t=5;t<6;t++)
	for(int t=0;t<7;t++)
	{

		//broken segment (Clover3's Right segment)
		if(s==2 && t==2)	continue;
		cout << "starting # " << s << " Ge analysis" << endl;


		char rootFileName[300];
		if(t<3)	sprintf(rootFileName,"%s/midas2nptool_Ge_PeakFitting_%s_cl%d_segment%d.root",iPath2,fName,s+1,t+1);
		else if(t<7)	sprintf(rootFileName,"%s/midas2nptool_Ge_PeakFitting_%s_cl%d_crystal%d.root",iPath2,fName,s+1,t-3+1);

		TFile *foutR0 = new TFile(rootFileName,"RECREATE");
		if ( foutR0->IsOpen() ) 
		{
			printf("root file opened for writing.\n");
		}	 
		else 
		{
    			printf("root file didn't open.\n");
    			exit (1);
		}

		if(t<3)	sprintf(brN,"GeSpec_cl%d_segment%d",s+1,t+1);
		else if(t<7)	sprintf(brN,"GeSpec_cl%d_crystal%d",s+1,t-3+1);
		h9[s][t] = new TH1I(brN, brN, 5000, 0, 5000);
		h9[s][t]->GetXaxis()->SetTitle("Channel (ch)");
		h9[s][t]->GetYaxis()->SetTitle("Counts");

		for(int k=0;k<5000;k++)
		{
			//segment
			if(t<3)
			{
				for(int l=0;l<gesegmentspec[s][t][k];l++)	h9[s][t]->Fill(k);
			}
			//crystal
			else if(t<7)
			{
				for(int l=0;l<gecrystalspec[s][t-3][k];l++)	h9[s][t]->Fill(k);
			}
		}

	  	cD0 = new TCanvas(brN,brN, 200,40,1150, 1150);
	  	cD0->GetFrame()->SetFillColor(21);
	  	cD0->GetFrame()->SetBorderMode(1);
	  	cD0->GetFrame()->SetBorderSize(200); 

//  for (int i = 0; i < 32; i++ ) {
 //   hIn[i]->SetAxisRange(1000,5000,"X");

    		TSpectrum *sp = new TSpectrum(nPeaks);
    		nfound[s][t] = sp->Search(h9[s][t],2,"",sigma);
    		//    cout << nfound[i] << endl;
    		if (nfound[s][t] < nPeaks ) {
      			nfound[s][t] = sp->Search(h9[s][t],2,"",sigma/2);
    		}

		for(int i=0;i<nfound[s][t];i++)	Hge[s][t][i]=sp->GetPositionX()[i];



/////////////////////////// IMPORTANT!!!!! Use this for 152Eu calibration because it is very difficulult to find peaks correctly using TSpectrum //////////////////////////////////////
/*
	//for 152Eu + 60Co

	//For 152Eu, you may need to set the rough estimated ADC value here for each peak.
	//set nfound[s][t] too, because probably TSpectrum's peak search is not working correctly.

		//4 is the number of peaks I set for 152Eu, but you can find and fit peaks as many as you like.
		//nfound[s][t] = 4;
		nfound[s][t] = EventNum;

		//Comment by Shuya 170320. Give the exact peak position! (the most frequent position).
		//Clover 1
		if(s==0)
		{
			if(t==0)
			{
				//clover#==s, channel==t, peak number=0-3
				Hge[s][t][0] = 148;	
				Hge[s][t][1] = 234;	
				Hge[s][t][2] = 304;	
				Hge[s][t][3] = 1049;	
			}
			else if(t==1)
			{
				Hge[s][t][0] = 164;	
				Hge[s][t][1] = 248;	
				Hge[s][t][2] = 316;	
				Hge[s][t][3] = 1039;	
			}
			else if(t==2)
			{
				Hge[s][t][0] = 140;	
				Hge[s][t][1] = 226;	
				Hge[s][t][2] = 297;	
				Hge[s][t][3] = 1048;	
			}
			else if(t==3)
			{
				Hge[s][t][0] = 169;	
				Hge[s][t][1] = 253;	
				Hge[s][t][2] = 320;	
				Hge[s][t][3] = 1041;	
			}
			else if(t==4)
			{
				Hge[s][t][0] = 156;	
				Hge[s][t][1] = 241;	
				Hge[s][t][2] = 310;	
				Hge[s][t][3] = 1045;	
			}
			else if(t==5)
			{
				Hge[s][t][0] = 154;	
				Hge[s][t][1] = 239;	
				Hge[s][t][2] = 309;	
				Hge[s][t][3] = 1049;	
			}
			else if(t==6)
			{
				Hge[s][t][0] = 162;	
				Hge[s][t][1] = 246;	
				Hge[s][t][2] = 314;	
				Hge[s][t][3] = 1044;	
			}
		}

		//Clover 2
		else if(s==1)
		{
			if(t==0)
			{
				Hge[s][t][0] = 140;	
				Hge[s][t][1] = 227;	
				Hge[s][t][2] = 296;	
				Hge[s][t][3] = 1045;	
			}
			else if(t==1)
			{
				Hge[s][t][0] = 140;	
				Hge[s][t][1] = 227;	
				Hge[s][t][2] = 298;	
				Hge[s][t][3] = 1048;	
			}
			else if(t==2)
			{
				Hge[s][t][0] = 141;	
				Hge[s][t][1] = 227;	
				Hge[s][t][2] = 298;	
				Hge[s][t][3] = 1046;	
			}
			else if(t==3)
			{
				Hge[s][t][0] = 179;	
				Hge[s][t][1] = 261;	
				Hge[s][t][2] = 327;	
				Hge[s][t][3] = 1035;	
			}
			else if(t==4)
			{
				Hge[s][t][0] = 172;	
				Hge[s][t][1] = 255;	
				Hge[s][t][2] = 322;	
				Hge[s][t][3] = 1040;	
			}
			else if(t==5)
			{
				Hge[s][t][0] = 146;	
				Hge[s][t][1] = 233;	
				Hge[s][t][2] = 302;	
				Hge[s][t][3] = 1048;	
			}
			else if(t==6)
			{
				Hge[s][t][0] = 162;	
				Hge[s][t][1] = 246;	
				Hge[s][t][2] = 314;	
				Hge[s][t][3] = 1041;	
			}
		}

		//Clover 3
		else if(s==2)
		{
			if(t==0)
			{
				//In some channel, the lowest peak of 152Eu is below threshold and no data.
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 260;	
				Hge[s][t][2] = 327;	
				Hge[s][t][3] = 1043;	
			}
			else if(t==1)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 243;	
				Hge[s][t][2] = 311;	
				Hge[s][t][3] = 1042;	
			}
			else if(t==3)
			{
				Hge[s][t][0] = 141;	
				Hge[s][t][1] = 228;	
				Hge[s][t][2] = 298;	
				Hge[s][t][3] = 1049;	
			}
			else if(t==4)
			{
				Hge[s][t][0] = 167;	
				Hge[s][t][1] = 251;	
				Hge[s][t][2] = 319;	
				Hge[s][t][3] = 1044;	
			}
			else if(t==5)
			{
				Hge[s][t][0] = 145;	
				Hge[s][t][1] = 231;	
				Hge[s][t][2] = 302;	
				Hge[s][t][3] = 1048;	
			}
			else if(t==6)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 247;	
				Hge[s][t][2] = 315;	
				Hge[s][t][3] = 1042;	
			}
		}

		//Clover 4
		else if(s==3)
		{
			if(t==0)
			{
				Hge[s][t][0] = 143;	
				Hge[s][t][1] = 231;	
				Hge[s][t][2] = 303;	
				Hge[s][t][3] = 1068;	
			}
			else if(t==1)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 236;	
				Hge[s][t][2] = 306;	
				Hge[s][t][3] = 1046;	
			}
			else if(t==2)
			{
				Hge[s][t][0] = 177;	
				Hge[s][t][1] = 260;	
				Hge[s][t][2] = 329;	
				Hge[s][t][3] = 1050;	
			}
			else if(t==3)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 235;	
				Hge[s][t][2] = 304;	
				Hge[s][t][3] = 1046;	
			}
			else if(t==4)
			{
				Hge[s][t][0] = 163;	
				Hge[s][t][1] = 248;	
				Hge[s][t][2] = 316;	
				Hge[s][t][3] = 1046;	
			}
			else if(t==5)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 245;	
				Hge[s][t][2] = 314;	
				Hge[s][t][3] = 1044;	
			}
			else if(t==6)
			{
				Hge[s][t][0] = 0;	
				Hge[s][t][1] = 245;	
				Hge[s][t][2] = 313;	
				Hge[s][t][3] = 1044;	
			}
		}
*/
///Cut until here if you don't use this code for 152Eu.
		

	//Draw the peak search results.
		h9[s][t]->Draw();
  		char sName[500];
		if(t<3)	sprintf(imgName,"PeakSearch_%s_cl%d_segment%d", fName, s+1,t+1);
		else if(t<7)	sprintf(imgName,"PeakSearch_%s_cl%d_crystal%d", fName, s+1,t-3+1);
		sprintf(ext,"png");
		sprintf(sName,"%s/%s.%s",iPath,imgName,ext);
    		cD0->SaveAs(sName);


	//Sorting the peak search results (because sometimes the found peaks's order is not in order)

//sorting routine (by bubble sort method)
//for(i=0;i<18;i++)     cout << x[i] << " " << y[i] << " " << ex[i] << " " << ey[i] << endl;

	      int p, q, flag = 1;    // set flag to 1 to start first pass
	      Double_t temp, temp2, temp3, temp4;             // holding variable
	      int numLength = EventNum;
	      for(int i = 1; (i <= numLength) && flag; i++)
     	      {
	          flag = 0;
        	  for (int j=0; j < (numLength -1); j++)
	         {
        	       //if (Hge[s][t][j+1] > Hge[s][t][j])      // ascending order simply changes to <
        	       if (Hge[s][t][j+1] < Hge[s][t][j])      // ascending order simply changes to <
             		 {
              		      temp = Hge[s][t][j];             // swap elements
                	    Hge[s][t][j] = Hge[s][t][j+1];
                	    Hge[s][t][j+1] = temp;
                  	  flag = 1;               // indicates that a swap occurred.
          	         }
         	 }
        	}


//////////////////////////////// MAKING HISTOGRAM FOR EACH FOUND PEAK //////////////////////////////////////////
//Making histograms around the peak position (+- 100 ch)

		//for(j=0;j<EventNum;j++)
		for(int j=0;j<nfound[s][t];j++)
		{


			if(t<3)	sprintf(brN,"E%dV_GeSpec_Cl%d_Segment%d",(int)Ege[j], s+1, t+1);
			else if(t<7)	sprintf(brN,"E%dV_GeSpec_Cl%d_Crystal%d",(int)Ege[j], s+1, t-3+1);
			h11[s][t][j] = new TH1I(brN, brN, 200, Hge[s][t][j]-100, Hge[s][t][j]+100);
			h11[s][t][j]->GetXaxis()->SetTitle("Channel (ch)");
			h11[s][t][j]->GetYaxis()->SetTitle("Counts/ch");

			for(k=Hge[s][t][j]-100;k<Hge[s][t][j]+100;k++)
			{
				//segment
				if(t<3)
				{
					for(int l=0;l<gesegmentspec[s][t][k];l++)	h11[s][t][j]->Fill(k);
				}
				//crystal
				else if(t<7)
				{
					for(int l=0;l<gecrystalspec[s][t-3][k];l++)	h11[s][t][j]->Fill(k);
				}
			}


			TF1	*f1;
			TF1	*fga1;
			//by Shuya 160320. Replaced with skewed gaussian fit
			f1 = new TF1("gausfit", gausfit, Hge[s][t][j]-15, Hge[s][t][j]+15, 4);
			fga1 = new TF1("fga1", "[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]) )+[3]", Hge[s][t][j]-15, Hge[s][t][j]+15);
			//For Skewed Gaussian, try followings.
			//f1 = new TF1("skewedgausfit", skewedgausfit, Hge[s][t][j]-15, Hge[s][t][j]+15, 5);
			//fga1 = new TF1("fga1", "[0]*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]) )+[3]+[4]*x", Hge[s][t][j]-15, Hge[s][t][j]+15);

			gStyle->SetOptFit(111);
			gStyle->SetOptStat(1001111);

			if(j==0)
			{
			    sprintf(cName,"cQ%i",2*s);
			    //cD[(j-90)/6] = new TCanvas(cName,cName, 200,40,840,840);
			    cD[2*s] = new TCanvas(cName,cName, 200,40,1150, 1150);
			    cD[2*s]->GetFrame()->SetFillColor(21);
			    cD[2*s]->GetFrame()->SetBorderMode(1);
			    cD[2*s]->GetFrame()->SetBorderSize(200); 
			    //cD[2*s]->Divide(2,3);
			    cD[2*s]->Divide(2,4);
			}

			//else if(j==6)
			else if(j==8)
			{
			    sprintf(cName,"cQ%i",2*s+1);
			    //cD[(j-90)/6] = new TCanvas(cName,cName, 200,40,840,840);
			    cD[2*s+1] = new TCanvas(cName,cName, 200,40,1150, 1150);
			    cD[2*s+1]->GetFrame()->SetFillColor(21);
			    cD[2*s+1]->GetFrame()->SetBorderMode(1);
			    cD[2*s+1]->GetFrame()->SetBorderSize(200); 
			    //cD[2*s+1]->Divide(2,3);
			    cD[2*s+1]->Divide(2,4);
			}

			//if(j/6==0)	cD[2*s]->cd(j%6+1);
			//else if(j/6==1)	cD[2*s+1]->cd(j%6+1);
			//if(j/6==0)	cD[2*s]->cd(j%6+1);
			//else if(j/6==1)	cD[2*s+1]->cd(j%6+1);
			if(j/8==0)	cD[2*s]->cd(j%8+1);
			else if(j/8==1)	cD[2*s+1]->cd(j%8+1);

	

//////////////////////////////// PEAK FITTING PART //////////////////////////////////////////
//Fix the fitting setting according to your need. 

			/*
			bg_low = 0.0;
			bg_high = 0.0;
			bg_ave = 0.0;
			*/


			//segment
			if(t<3)
			{
				//Evaluating the BG from +-30 ch off the peaks. 
				bg_low = gesegmentspec[s][t][Hge[s][t][j]-30];
				bg_high = gesegmentspec[s][t][Hge[s][t][j]+30];
				//case the low energy gamma nearby the threshold (like 121 keV 152Eu)
				if(bg_low == 0)	bg_low = bg_high;
				bg_ave = (bg_low + bg_high)/2.0;

				//by Shuya 170320
				if(gesegmentspec[s][t][Hge[s][t][j]]>bg_ave)	f1->SetParameters((double)gesegmentspec[s][t][Hge[s][t][j]]-((double)bg_low+(double)bg_high)/2.0, (double)Hge[s][t][j]+0.5, 2.0, ((double)bg_low+(double)bg_high)/2.0);
				else	f1->SetParameters(1, (double)Hge[s][t][j]+0.5, 2.0, 0.0);
				//if(gesegmentspec[s][t][Hge[s][t][j]]>bg_ave)	f1->SetParameters((double)gesegmentspec[s][t][Hge[s][t][j]]-((double)bg_low+(double)bg_high)/2.0, (double)Hge[s][t][j]+0.5, 2.0, bg_ave,((double)bg_high-(double)bg_low)/60.0);
				//else	f1->SetParameters(1, (double)Hge[s][t][j]+0.5, 2.0, 0.0,1.0);

				f1->SetParLimits(0, (double)gesegmentspec[s][t][Hge[s][t][j]]-bg_ave-200.0, (double)gesegmentspec[s][t][Hge[s][t][j]]-bg_ave+200);
			}

			//crystal
			else if(t<7)
			{
				//Evaluating the BG from +-30 ch off the peaks. 
				bg_low = gecrystalspec[s][t-3][Hge[s][t][j]-15];
				bg_high = gecrystalspec[s][t-3][Hge[s][t][j]+15];
				//case the low energy gamma nearby the threshold
				if(bg_low == 0)	bg_low = bg_high;
				bg_ave = (bg_low + bg_high)/2.0;

				//by Shuya 170320
				if(gecrystalspec[s][t-3][Hge[s][t][j]]>bg_ave)	f1->SetParameters((double)gecrystalspec[s][t-3][Hge[s][t][j]]-((double)bg_low+(double)bg_high)/2.0, (double)Hge[s][t][j]+0.5, 2.0, ((double)bg_low+(double)bg_high)/2.0);
				else	f1->SetParameters(1, (double)Hge[s][t][j]+0.5, 2.0, 0.0);
				//if(gecrystalspec[s][t-3][Hge[s][t][j]]>bg_ave)	f1->SetParameters((double)gecrystalspec[s][t-3][Hge[s][t][j]]-((double)bg_low+(double)bg_high)/2.0, (double)Hge[s][t][j]+0.5, 2.0, bg_ave,((double)bg_high-(double)bg_low)/30.0);
				//else	f1->SetParameters(1, (double)Hge[s][t][j]+0.5, 2.0, 0.0,1.0);

				f1->SetParLimits(0, (double)gecrystalspec[s][t-3][Hge[s][t][j]]-bg_ave-100, (double)gecrystalspec[s][t-3][Hge[s][t][j]]-bg_ave+100);
			}


			f1->SetParLimits(1, Hge[s][t][j]-10.0,Hge[s][t][j]+10.0);
			f1->SetParLimits(2, 1.0, 5.0);
			//f1->FixParameter(0, (double)gammaspec[j][Egamma[s]]);
			//f1->FixParameter(1, Hgamma[s][j]+0.5);
			//f1->SetParLimits(1, Hge[s][t][j]-3.0,Hge[s][t][j]+3.0);
			//f1->SetParLimits(2, 1.0, 10.0);
			//f1->FixParameter(2, 2.0);
			//f1->FixParameter(3, ((double)bg_low+(double)bg_high)/2.0);

			f1->SetParNames("Constant", "Centroid", "Sigma", "Background");

			//by Shuya 170320
			if(t<3)	h11[s][t][j]->Fit("gausfit", "Q", "", Hge[s][t][j]-15, Hge[s][t][j]+15);
			else if(t<7)	h11[s][t][j]->Fit("gausfit", "Q", "", Hge[s][t][j]-10, Hge[s][t][j]+10);
			//if(t<3)	h11[s][t][j]->Fit("skewedgausfit", "Q", "", Hge[s][t][j]-15, Hge[s][t][j]+15);
			//else if(t<7)	h11[s][t][j]->Fit("skewedgausfit", "Q", "", Hge[s][t][j]-10, Hge[s][t][j]+10);

			h11[s][t][j]->Draw();

			//Amplitude
			double ga1 = f1->GetParameter(0);
			//Peak position
			double ga2 = f1->GetParameter(1);
			//Standard deviation
			double ga3 = f1->GetParameter(2);
			double chi2 = f1->GetChisquare();
			//Peak position error
			double ga4 = f1->GetParError(1);
			//Standard deviation error
			double ga5 = f1->GetParError(2);
	
			fga1->SetParameter(0, ga1);
			fga1->SetParameter(1, ga2);
			fga1->SetParameter(2, ga3);
			fga1->SetParameter(3, ga4);
			fga1->SetLineColor(4);

			//Data points to make ADC channel (x) vs Energy (y) plot.
			x[j]=ga2;
			ex[j]=ga4;
			y[j]=Ege[j];
			ey[j]=ege[j];
			z[j]=ga3;
			ez[j]=ga5;

	//Comment by Shuya 161108. End of nfound[] loop
		}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// CALIBRATION PART //////////////////////////////////////////
//Note!! This is a temporarary calibration using a single run (or single source). You should use another code to merge all the peak data from several runs.


		TF1	*lFit;
		char	fitName[400];
		double	fVals[2];

		//These are output files for each channel.
		if(t<3)	sprintf(cName,"%s/%s_cl%d_segment%d_calsource.dat",iPath2,fName, s+1,t+1);
		else if(t<7)	sprintf(cName,"%s/%s_cl%d_crystal%d_calsource.dat",iPath2,fName, s+1,t-3+1);
		fout2 = fopen(cName,"w");
		printf("Opened the output data file2\n");

	    //  for(i=0;i<2;i++)
	      for(i=0;i<1;i++)
		{
			//sprintf(imgName,"Calibration_%s-%s_ch%d_%d",finname0,finname1, s, i);
			if(t<3)	sprintf(imgName,"Calibration_%s_cl%d_segment%d_%d", fName, s+1, t+1,i);
			else if(t<7)	sprintf(imgName,"Calibration_%s_cl%d_crystal%d_%d", fName, s+1, t-3+1,i);
		        sprintf(ext,"png");
		        sprintf(sName,"%s/%s.%s",iPath,imgName,ext);
    			cD[2*s+i]->SaveAs(sName);
		}


		TCanvas *cD2;
  		TGraphErrors *grCal;
       		char tgName[100];

		if(t<3)	sprintf(cName,"%s_cl%d_segment%d_calibrationsource",fName,s+1,t+1);
		else if(t<7)	sprintf(cName,"%s_cl%d_crystal%d_calibrationsource",fName,s+1,t-3+1);
		cD2 = new TCanvas(cName,cName, 200,40,1150, 1150);
		cD2->GetFrame()->SetFillColor(21);
		cD2->GetFrame()->SetBorderMode(1);
		cD2->GetFrame()->SetBorderSize(200); 
		//if(((j-90)/6)!=35)	cD[(j-90)/6]->Divide(2,3);
		
		//grCal = new TGraphErrors(EventNum,x, y, ex, ey);
		grCal = new TGraphErrors(nfound[s][t],x, y, ex, ey);
    		if(t<3)	sprintf(tgName,"%s_cl%d_segment%d_calibrationsource",fName,s+1,t+1);
    		else if(t<7)	sprintf(tgName,"%s_cl%d_crystal%d_calibrationsource",fName,s+1,t-3+1);
		//grCal->GetXaxis()->SetTitle("Pulser Voltage (mV)");
		grCal->GetXaxis()->SetTitle("ADC channel (ch)");
		grCal->GetYaxis()->SetTitle("Gamma Energy (keV)");


	//by Shuya 160923
                sprintf(fitName,"lFit_%i_%i",s+1,t+1);
                //lFit = new TF1(fitName,"pol1",0,4000.0);
                ///lFit = new TF1(fitName,"[0]+[1]*x",0,4092.0);

	//You need to adjust depending on files
                //lFit = new TF1(fitName,"[0]+[1]*x",1500,10000.0);
                //lFit = new TF1(fitName,"[0]+[1]*x",500,10000.0);
                lFit = new TF1(fitName,"[0]+[1]*x",150.,5000.0);
                lFit->SetParameter(0,0);
                lFit->SetParameter(1,0.3);
                //lFit->SetParameter(3,0.0);
                grCal->Fit(fitName,"RNQ");
                //grCal->Fit(fitName,"Q","", x[0], x[17]);
                lFit->GetParameters(&fVals[0]);
                //fVals[0]=lFit->GetParameter(0);
                //fVals[1]=lFit->GetParameter(1);
                //cout << fVals[0] << " " << fVals[1] << endl;
                cout << fVals[0] << " " << fVals[1] << endl;
                //cout << fVals[0] << " " << fVals[1] << " " << fVals[2] << " " << fVals[3] << endl;
                //grCal->GetFunction("pol1")->SetLineWidth(2);


                grCal->SetTitle(tgName);
                grCal->SetLineColor(1);
                grCal->SetLineWidth(2);
                grCal->SetMarkerColor(1);
                grCal->SetMarkerStyle(21);


	  	//grCal[k]->Draw("AP");
	  	grCal->Draw("AP");
                //lFit->Draw("same");

	//Adjust depending on files
		//TF1 *lFit2 = new TF1(fitName,"[0]+[1]*x",0.,10000.);
		TF1 *lFit2 = new TF1(fitName,"[0]+[1]*x",0.,5000.);
                lFit2->SetLineColor(kRed);
		lFit2->SetParameter(0,fVals[0]);
		lFit2->SetParameter(1,fVals[1]);
		lFit2->Draw("same");
	//Adjust depending on files
		//grCal->GetXaxis()->SetLimits(0.,10000.);
		grCal->GetXaxis()->SetLimits(0.,5000.);
		//grCal->GetYaxis()->SetRangeUser(0.,1000.);
		grCal->GetYaxis()->SetRangeUser(0.,10000.);
		cD2->Update();

      		//sprintf(imgName,"%s-%s_ch%d",finname0,finname1,s);
      		if(t<3)	sprintf(imgName,"CalibrationResult_%s_cl%d_segment%d",fName,s+1,t+1);
      		else if(t<7)	sprintf(imgName,"CalibrationResult_%s_cl%d_crystal%d",fName,s+1,t-3+1);
	        sprintf(ext,"png");
   	        sprintf(sName,"%s/%s.%s",iPath,imgName,ext);
      		cD2->SaveAs(sName);


		//Peak and corresponding energy are saved in this output files for each channel (crystal and segment). 
		//for(i=0;i<EventNum;i++)	fprintf(fout2, "%lf %lf %lf %lf\n", x[i], y[i], ex[i], ey[i]);
		for(i=0;i<nfound[s][t];i++)	fprintf(fout2, "%lf %lf %lf %lf %lf %lf\n", x[i], y[i], z[i], ex[i], ey[i], ez[i]);
		fclose(fout2);

		//This is an output file for the calibration parameters (a,b) temporarily obtained using this run.
		fprintf(fout, "%i %i %lf %lf\n", s, t, fVals[0],fVals[1]);

	   foutR0->Write();
	   foutR0->Close();

	//this is for t's loop
	}

//this is for if(s==XX)
  }

}

fclose(fout);
return 0;
}
