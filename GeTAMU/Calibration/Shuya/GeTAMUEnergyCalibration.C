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


int midas2nptool_Ge_EnergyCalibration()
{

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


printf("Vector creation complete\n");

int fileloop = 0;
int filenum = 0;
int ncols = 0;
int strlength = 0;
FILE *batch;
char tempname[400] = "";
char subdir[400] = "";
char foutname[400] = "";
//by Shuya 140610
char finname0[400] = "";
char finname1[400] = "";
int	bg_low=0;
int	bg_high=0;
double	bg_ave=0;

int ch[32] = {0};
int ch2[32] = {0};


double	x[100], ex[100],y[100],ey[100],z[100],ez[100];


///////////////////////////////NOTE!! YOU HAVE TO CHANGE SOME VARIABLES DEPENDING ON YOUR NEED ///////////////////////////////////////////////
////Search "Change" will take you to these variables.



////////////////////////////////////CHANGE DEPENDING ON YOUR NEED PART 1 OUT OF 2/////////////////////////////////////////////////////////////
//Change the outputfile name as you like. This Output gives you a table of calibration coefficients for each Ge channel. 

//Output file
FILE *fout;
//sprintf(foutname, "/home/shuyaota/midas2nptool/Outputs/EXPT4/ER189_0-190_0-191_0_FinalCalibrationParam.dat");
sprintf(foutname, "/home/shuyaota/midas2nptool/Outputs/TEST/ER189_0-190_0-191_0_FinalCalibrationParam.dat");
fout = fopen(foutname,"w");
printf("Data file to process %s\n", foutname);
	        
//Output figures
//char iPath[300] = "/home/shuyaota/midas2nptool/Graphs/EXPT4";
char iPath[300] = "/home/shuyaota/midas2nptool/Graphs/TEST";
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//4 clovers
for(int s=0;s<4;s++)
{


	//segment 0-2, core 3-6
	for(int t=0;t<7;t++)
	{
		//Dead channel (CL3, Segment Right).
		if(s==2 && t==2)
		{
			fprintf(fout, "%i %i %lf %lf\n", s, t, 0.0,0.0);
			continue;
		}


	        char fileName[100];
		double adc[32] = {0.};
		double adc_err[32] = {0.};
		double offset[32] = {0.};
		double slope[32] = {0.};
		double mean[32] = {0.};
		double energy[32] = {0.};
		double energy_err[32] = {0.};
		double width[32] = {0.};
		double width_err[32] = {0.};

		double	zero_offset[32] = {0.};

		int w=0;

////////////////////////////////////CHANGE DEPENDING ON YOUR NEED PART 2 OUT OF 3/////////////////////////////////////////////////////////////
		//Change number of files depending on your need. 
		filenum=2;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(int u=0;u<filenum;u++)
		{
	

////////////////////////////////////CHANGE DEPENDING ON YOUR NEED PART 3 OUT OF 3/////////////////////////////////////////////////////////////
			//Change input filenames (output from midas2nptool_Ge_PeakFitting.C)
			if(u==0)	sprintf(fileName,"ER189_0");
			else if(u==1)	sprintf(fileName,"ER190_0");
			else if(u==2)	sprintf(fileName,"ER191_0");

       	 		FILE *fincal;
		        if(t<3)	sprintf(finname0, "/home/shuyaota/midas2nptool/Outputs/TEST/%s_cl%d_segment%d_calsource.dat", fileName,s+1,t+1);
		        else if(t<7)	sprintf(finname0, "/home/shuyaota/midas2nptool/Outputs/TEST/%s_cl%d_crystal%d_calsource.dat", fileName,s+1,t-3+1);
		        fincal = fopen(finname0,"r");
		        printf("Data file to process %s\n", finname0);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			while(!feof(fincal))
			{
       			 	fscanf(fincal,"%lf %lf %lf %lf %lf %lf", &adc[w], &energy[w], &width[w], &adc_err[w], &energy_err[w], &width_err[w]);
				w++;
			}
			w--;

       	 		fclose(fincal);
		}



	//sorting routine (by bubble sort method) (Because you are merging different data using different gamma sources.)

	      int flag = 1;    // set flag to 1 to start first pass
	      Double_t temp, temp2, temp3, temp4;             // holding variable
	      //int numLength = 8;
	      int numLength = w;
	      for(int i = 1; (i <= numLength) && flag; i++)
     	      {
	          flag = 0;
        	  for (int j=0; j < (numLength -1); j++)
	         {
        	       //if (Hge[s][t][j+1] > Hge[s][t][j])      // ascending order simply changes to <
        	       if (adc[j+1] < adc[j])      // ascending order simply changes to <
             		 {
              		      temp = adc[j];             // swap elements
                	    adc[j] = adc[j+1];
                	    adc[j+1] = temp;

              		      temp = energy[j];             // swap elements
                	    energy[j] = energy[j+1];
                	    energy[j+1] = temp;

              		      temp = width[j];             // swap elements
                	    width[j] = width[j+1];
                	    width[j+1] = temp;

              		      temp = adc_err[j];             // swap elements
                	    adc_err[j] = adc_err[j+1];
                	    adc_err[j+1] = temp;

              		      temp = energy_err[j];             // swap elements
                	    energy_err[j] = energy_err[j+1];
                	    energy_err[j+1] = temp;

              		      temp = width_err[j];             // swap elements
                	    width_err[j] = width_err[j+1];
                	    width_err[j+1] = temp;

                  	  flag = 1;               // indicates that a swap occurred.
          	         }
         	 }
        	}

		for(int i=0;i<numLength;i++)
		{
			x[i] = adc[i];
			ex[i] = adc_err[i];
			y[i] = energy[i];
			ey[i] = energy_err[i];
			z[i] = width[i];
			ez[i] = width_err[i];
		     cout << x[i] << " " << y[i] << " " << z[i] << " " << ex[i] << " " << ey[i] << " " << ez[i] << endl;
		}
		



		TF1	*lFit;
		TCanvas *cD2;
  		TGraphErrors *grCal;
      		char tgName[100];
      		char cName[100];
	      	char sName[100];
      		char ext[100];
	        char imgName[100];
		double	fVals[2];

		if(t<3)	sprintf(cName,"Ge_FinalCalibration_cl%d_segment%d", s+1, t+1);
		else if(t<7)	sprintf(cName,"Ge_FinalCalibration_cl%d_crystal%d", s+1, t-3+1);
		//sprintf(cName,"Ge_FinalCalibration_cl%d_crystal_",j);
		cD2 = new TCanvas(cName,cName, 200,40,1150, 1150);
		cD2->GetFrame()->SetFillColor(21);
		cD2->GetFrame()->SetBorderMode(1);
		cD2->GetFrame()->SetBorderSize(200); 
		
		//grCal = new TGraphErrors(numLength,x, y, ex, ey);
		//grCal = new TGraphErrors(8,x, y, ex, ey);
		grCal = new TGraphErrors(numLength,x, y, ex, ey);
		//sprintf(tgName,"Ge_FinalCalibration_cl%d_crystal%d",j);
		if(t<3)	sprintf(tgName,"Ge_FinalCalibration_cl%d_segment%d",s+1,t+1);
		else if(t<7)	sprintf(tgName,"Ge_FinalCalibration_cl%d_crystal%d",s+1,t-3+1);
		grCal->GetXaxis()->SetTitle("ADC channel (ch)");
		grCal->GetYaxis()->SetTitle("Energy (MeV)");

                if(t<3)	sprintf(cName,"lFit_cl%d_segment%d",s+1,t+1);
                else if(t<7)	sprintf(cName,"lFit_cl%d_crystal%d",s+1,t-3+1);
                lFit = new TF1(cName,"[0]+[1]*x",0.,5000.);
                //lFit->SetParameter(0,zero_offset[i]);
                //lFit->SetParameter(1,slope[i]);

                lFit->SetParameter(0,0);
                lFit->SetParameter(1,0.3);
                grCal->Fit(cName,"RNQ");
                //grCal->Fit(fitName,"Q","", x[0], x[17]);
                lFit->GetParameters(&fVals[0]);
                //fVals[0]=lFit->GetParameter(0);
                //fVals[1]=lFit->GetParameter(1);
                //cout << fVals[0] << " " << fVals[1] << endl;
                cout << fVals[0] << " " << fVals[1] << endl;

                grCal->SetTitle(tgName);
                grCal->SetLineColor(1);
                grCal->SetLineWidth(2);
                grCal->SetMarkerColor(1);
                grCal->SetMarkerStyle(21);

	  	grCal->Draw("AP");
	  	lFit->Draw("same");

		grCal->GetXaxis()->SetLimits(0.,5000.);
		//grCal->GetYaxis()->SetLimits(0.,15.);
		grCal->GetYaxis()->SetRangeUser(0.,10000.);
		cD2->Update();

		//sprintf(imgName,"Ge_FinalCalibration_cl%d_crystal%d",j);
		if(t<3)	sprintf(imgName,"Ge_FinalCalibration_cl%d_segment%d",s+1,t+1);
		else if(t<7)	sprintf(imgName,"Ge_FinalCalibration_cl%d_crystal%d",s+1,t-3+1);
	        sprintf(ext,"png");
   	        sprintf(sName,"%s/%s.%s",iPath,imgName,ext);
      		cD2->SaveAs(sName);

		fprintf(fout, "%i %i %lf %lf\n", s, t, fVals[1],fVals[0]);
		}
	}

	fclose(fout);
	return 0;
}
