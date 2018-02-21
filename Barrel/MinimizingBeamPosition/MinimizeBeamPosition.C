/*
This code is used to minimise the position of the beam spot
from observing the inelastic channel in the barrel at a specific energy

input: 
1. Angle of emission 
(typically from an inelastic kinematic line calculation 
at the right beam energy)
2. text file of the observed
strip angle fwhm_angle

*/

// ROOT headers
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TVector3.h"
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
#include <vector>
#include <string>

//NPTool headers
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"

#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

using namespace::std;


class clsLine3D{
	//members
  public:
	TVector3 fPoint; // point on the line
	TVector3 fDir; // directing vector
	
	public:
	//functions
	clsLine3D(TVector3 p,TVector3 l){
		fDir = l; 
		fPoint = p;
		}
	~clsLine3D();


	TVector3 GetPoint(double t){ 
		return (t*fDir+fPoint);
		}
		
		
	double GetPointLineDistance(TVector3 p){
	  TVector3 v1 = fPoint-p; // vector to arb. point on line
	  TVector3 v2 = (v1.Dot(fDir))*fDir; 
	  v1 = v1-v2;
		return v1.Mag();
		}	
		
};


class clsPlane3D{
  public:
	//members
  TVector3 fPoint; // point on the line
	TVector3 fNorm; // directing vector
	public:
	//functions
	clsPlane3D(TVector3 p, TVector3 n){
	  fNorm = n;  //normal vector
		fPoint = p; // point on plane
		}

	clsPlane3D(TVector3 p1, TVector3 p2, TVector3 p3){
	  TVector3 p12 = p2-p1;
	  TVector3 p13 = p3-p1;
	  TVector3 n = p12.Cross(p13);
	  //
	  fNorm = n.Unit();  //normal vector
		fPoint = p1; // point on plane
		}
				
	~clsPlane3D();
};

//Global variables
  //fixed
  double gAngle;
  vector <int> gBarrelNumber;
  vector <int> gStripNumber;
  vector <double> gHyperStrip;
  vector <double> gPhiOfStrip;// Phi of a strip [0,360]
  vector <double> gDataAngle;
  vector <double> gDataAngleErr;
  vector <double> gDataHitOnStripZ; //Calculated from gDataAngle and assumed beam position
  vector <double> gDataHitOnStripZErrPlus; // this correspond to Angle - AngleSigma
  vector <double> gDataHitOnStripZErrMinus;
  vector <double> gDataHitPhi;      //Calculated from gDataAngle and assumed beam position
  vector <double> gNull;
  vector <double> gAssumedAngle;
  vector <double> gAssumedPhiAngle;
  vector <TVector3> gPointOnStrip;
  vector <TVector3> gHitOnStrip;
  vector <double> gHitOnStripZ;
  vector <double> gHitOnStripY;
  vector <double> gHitOnStripX;
  vector <double> gPointOnStripZ;
  vector <double> gPointOnStripY;
  vector <double> gPointOnStripX;
  
  //TH3F* gHist3D = new TH3F("3d","3d",100,-50,50,100,-50,50,200,-100,100);
    
  //fixed
  TVector3 gUserBeamSpot; 
  
  //final result
	TVector3 gFinalMinimPosition;
  TVector3 gFinalMinimPositionError;
  
//Functions
TVector3 GetLinePlaneIntersect(clsLine3D* l,clsPlane3D* p);
int Minimise(void); // numerical minimisation function which produces the final calibration parameters
double GetChiSquare(const double parameters[]); // returns a chi squared value 
TVector3 GetHitOnStrip(TVector3 Beamspot, int StripNumber, double angle); // calculate impact position on strip middle line
double GetHitThetaAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate angle from position assuming beam spot
double GetHitPhiAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate phi angle [0,360] from position assuming beam spot
TVector3 GetNormalOnDetector(double hyperstrip);
TVector3 GetPointOnStrip(double hyperstrip, double pos=0); // pos = [-1;+1]

// MAIN
void MinimizeBeamPosition(double Angle=-1, TString data="sample_BarrelAngles_Mg25.txt", double ax=0, double ay=0, double az=0){

  if(Angle==-1) {
    cout << " ERROR ---- Provide Angle ! ----- "<< endl; 
    cout << " To execute:\n\n .x MinimizeBeamPosition.C++( <angle-degree>, <input-ascii-file> ) \n"<< endl; 
    exit(-1); 
  }
  else {
    cout << endl; 
    cout << "Angle provided: " << Angle << endl;
    cout << "File provided: " << data << endl;
    cout << "Assumed Source position input (mm) : " << ax << " " << ay << " " << az << endl;
    cout << "NB0: The Assumed Source position input \"MUST\" be the position of the beam spot used\nduring the data analysis \"AND\" calculated relative to the geometrical barrel centroid.\n" <<endl; 
    }
  
  gAngle = Angle*TMath::DegToRad();
  gUserBeamSpot.SetXYZ(ax,ay,az);
  
// Read the data file and store values in c-vectors
  int barrel,strip, hstrip;
  double angle, sangle;
  string line;
  ifstream myfile (data.Data());
  if (myfile.is_open()) {
    while ( myfile>>barrel){
    if(barrel==0) {
      getline (myfile,line);
      continue ; // skip "comments" line starting with zero
      }
    myfile>>strip>>angle>>sangle;
    if(!(angle>0)) continue ; 
      gBarrelNumber.push_back(barrel);
      gStripNumber.push_back(strip);
      gHyperStrip.push_back((barrel-1)*4+strip);
      gDataAngle.push_back(angle);
      gDataAngleErr.push_back(sangle);
      gNull.push_back(0); // for plotting purposes
      //cout << " Reading " << barrel << " " << strip << " " << (barrel-1)*4+strip << " " <<  angle << " " << sangle << endl; 
    }
    myfile.close();
  }
  else { cout << "Unable to open file"; exit(-1);} 
 
 
  TNtuple *gHitScatterPlot = new TNtuple("gHitScatterPlot", "gHitScatterPlot", "x:y:z:hyperstrip");
  gHitScatterPlot->Fill(gUserBeamSpot.X(), gUserBeamSpot.Y(), gUserBeamSpot.Z(),-1);
   
// Construct the absolute positions of the central line along every strip
for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
  TVector3 temp = GetPointOnStrip(gHyperStrip[i],0);
  gPointOnStrip.push_back(temp);
  gPointOnStripZ.push_back(temp.Z());
  gPointOnStripY.push_back(temp.Y());
  gPointOnStripX.push_back(temp.X()); 
  double phi_360 = temp.Phi()*TMath::RadToDeg();
  if(temp.Y()<0)  phi_360= 360+phi_360;  
  gPhiOfStrip.push_back(phi_360);
  
  //Fill this point in the Ntuple 
  gHitScatterPlot->Fill(temp.X(), temp.Y(), temp.Z(),gHyperStrip[i]);
  
  //Calculate z and phi from data with respect to the assumed beam spot
    clsLine3D* CentralStripLine = new clsLine3D(temp, TVector3(0,0,1)/*Strips are always along this vector*/ ); //Calculate central line of a strip
    double distance = CentralStripLine->GetPointLineDistance(gUserBeamSpot);//Calculate the shortest (perp) distance from the beam spot to the line
    double z0 = distance/TMath::Tan(gDataAngle[i]*deg);//Calculate z wrt to beamspot z, and fill the absolute value (in case the beam spot is not at 0,0,0
    z0+=gUserBeamSpot.Z();
    gDataHitOnStripZ.push_back(z0);
   //repeat same procedure for Errors
    double z = distance/TMath::Tan((gDataAngle[i]-gDataAngleErr[i])*deg);
    z+=gUserBeamSpot.Z();
    gDataHitOnStripZErrPlus.push_back(fabs(z-z0));
    z = distance/TMath::Tan((gDataAngle[i]+gDataAngleErr[i])*deg);
    z+=gUserBeamSpot.Z();
    gDataHitOnStripZErrMinus.push_back(fabs(z-z0));
    
    
  //Fill this point in the Ntuple 
  temp = GetPointOnStrip(gHyperStrip[i],-1);
  gHitScatterPlot->Fill(temp.X(), temp.Y(), temp.Z(),gHyperStrip[i]);
  temp = GetPointOnStrip(gHyperStrip[i],+1);
  gHitScatterPlot->Fill(temp.X(), temp.Y(), temp.Z(),gHyperStrip[i]);
  
  
  //inspect
  //cout << gHyperStrip[i] << "  " <<  temp.Phi()*TMath::RadToDeg()<< endl ; 
}


// perform the numerical minimisation
  int result = 0 ;
  result = Minimise(); 
  cout << "Minimizing Result " << result << endl << endl;
  
  //Fill in the 3D after minimisation
  gHitScatterPlot->Fill(gFinalMinimPosition.X(), gFinalMinimPosition.Y(), gFinalMinimPosition.Z(),0);
  for (unsigned i=0; i < gHitOnStrip.size() ; i++)
    gHitScatterPlot->Fill(gHitOnStrip[i].X(), gHitOnStrip[i].Y(), gHitOnStrip[i].Z(), gHyperStrip[i]);
    
  //Show results
  TCanvas* canInspection= new TCanvas("canInspection","canInspection",650,650);
  canInspection->Divide(2,2);
  
  TGraphErrors* grStripvsData = new TGraphErrors(gDataAngle.size(),&gHyperStrip[0],&gDataAngle[0],&gNull[0],&gDataAngleErr[0]);
  grStripvsData->SetTitle("Theta (minimised) vs HyperStrip");
  grStripvsData->SetMarkerColor(kBlack);
  grStripvsData->SetMarkerStyle(33);
  grStripvsData->SetMarkerSize(2.3);
  
  TGraphErrors* grStripvsAssumedData = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gAssumedAngle[0],&gNull[0],&gNull[0]);
  grStripvsAssumedData->SetMarkerColor(kRed);
  grStripvsAssumedData->SetMarkerStyle(20);
  
  TGraph* grThetavsPhi = new TGraphErrors(gAssumedAngle.size(),&gAssumedPhiAngle[0],&gAssumedAngle[0]);
  grThetavsPhi->SetTitle("Theta (minimised) vs Phi (minimised)");
  grThetavsPhi->SetMarkerColor(kRed);
  grThetavsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsPhi = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPhiOfStrip[0]);
  grStripvsPhi->SetTitle(" Phi (of Strip) vs HyperStrip");
  grStripvsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsZ = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripZ[0]);
  grStripvsZ->SetTitle(" HitPosition(Z) vs HyperStrip");
  grStripvsZ->SetMarkerStyle(20);
  grStripvsZ->SetMarkerColor(kRed);
  
  TGraphAsymmErrors* grStripvsDataZ = new TGraphAsymmErrors(gAssumedAngle.size(),&gHyperStrip[0],&gDataHitOnStripZ[0],&gNull[0],&gNull[0],&gDataHitOnStripZErrMinus[0],&gDataHitOnStripZErrPlus[0]);
  grStripvsDataZ->SetTitle(" HitPosition(Z) (Data and minimized) vs HyperStrip");
  grStripvsDataZ->SetMarkerStyle(33);
  grStripvsDataZ->SetMarkerColor(kBlack);
  grStripvsDataZ->SetMarkerSize(2.3);


  TGraph* grStripvsMinimizedPhi = new TGraphErrors(gAssumedAngle.size(),&gPhiOfStrip[0],&gAssumedPhiAngle[0]);
  grStripvsMinimizedPhi->SetTitle(" Phi (minimised) vs Phi (of Strip)");
  grStripvsMinimizedPhi->SetMarkerStyle(20);
  grStripvsMinimizedPhi->SetMarkerColor(kRed);

   TGraph* grHyperStripvsMinimizedPhi = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gAssumedPhiAngle[0]);
  grHyperStripvsMinimizedPhi->SetTitle(" Phi (minimised) vs Phi (of Strip)");
  grHyperStripvsMinimizedPhi->SetMarkerStyle(20);
  grHyperStripvsMinimizedPhi->SetMarkerColor(kRed);
   
  TGraph* grStripvsY = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripY[0]);
  grStripvsY->SetTitle(" HitPosition(Y) vs HyperStrip");
  grStripvsY->SetMarkerStyle(20);
  grStripvsY->SetMarkerColor(kRed);
  
  TGraph* grStripvsX = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripX[0]);
  grStripvsX->SetTitle(" HitPosition(X) vs HyperStrip");
  grStripvsX->SetMarkerStyle(20);
  grStripvsX->SetMarkerColor(kRed);
  
  TGraph* grStripvsPointZ = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPointOnStripZ[0]);
  grStripvsPointZ->SetTitle(" Point on Strip Position(Z) vs HyperStrip");
  grStripvsPointZ->SetMarkerStyle(28);
  grStripvsPointZ->SetMarkerColor(3);

  TGraph* grStripvsPointY = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPointOnStripY[0]);
  grStripvsPointY->SetTitle(" Point on Strip Position(Y) vs HyperStrip");
  grStripvsPointY->SetMarkerStyle(20);
  grStripvsPointY->SetMarkerColor(3);
  
  TGraph* grStripvsPointX = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPointOnStripX[0]);
  grStripvsPointX->SetTitle(" Point on Strip Position(X) vs HyperStrip");
  grStripvsPointX->SetMarkerStyle(20);
  
  TGraphAsymmErrors* grDataZvsPhi = new TGraphAsymmErrors(gAssumedAngle.size(),&gPhiOfStrip[0],&gDataHitOnStripZ[0],&gNull[0],&gNull[0],&gDataHitOnStripZErrMinus[0],&gDataHitOnStripZErrPlus[0]);
  grDataZvsPhi->SetTitle("HitPosition(Z) (Data and minimised) vs Phi (of Strip)");
  grDataZvsPhi->SetMarkerStyle(33);
  grDataZvsPhi->SetMarkerColor(kBlack);
  grDataZvsPhi->SetMarkerSize(2.3);
  
 
  TGraph* grMinimZvsPhi = new TGraphErrors(gAssumedAngle.size(),&gPhiOfStrip[0],&gHitOnStripZ[0]);
  grMinimZvsPhi->SetTitle(" HitPosition(Z) vs HyperStrip");
  grMinimZvsPhi->SetMarkerStyle(20);
  grMinimZvsPhi->SetMarkerColor(kRed);
  
  canInspection->cd(1);
  //grStripvsMinimizedPhi->Draw("ap");
  grHyperStripvsMinimizedPhi->Draw("ap");
  canInspection->cd(2);
  grStripvsPhi->Draw("ap");
  canInspection->cd(3);
  grStripvsY->Draw("alp");
  canInspection->cd(4);
  grStripvsX->Draw("alp");
  /*
  canInspection->cd(7);
  grStripvsPointZ->Draw("alp");
  canInspection->cd(8);
  grStripvsPointY->Draw("alp");
  canInspection->cd(9);
  grStripvsPointX->Draw("alp");
  grStripvsPointY->Draw("same p");
  */
  canInspection->Draw();
    
  TCanvas* canResult= new TCanvas("canResult","canResult",650,650);
  canResult->Divide(2,2);
  canResult->cd(1);
  grStripvsData->Draw("alp"); 
  grStripvsAssumedData->Draw("p same"); 
  canResult->cd(2);
  grThetavsPhi->Draw("ap");
  canResult->cd(3);
  grStripvsDataZ->Draw("alp");
  grStripvsZ->Draw("p same");
  canResult->cd(4);
  grDataZvsPhi->Draw("ap");
  grMinimZvsPhi->Draw("p same");
  canResult->Draw();


  TCanvas* canInspection3D= new TCanvas("scatter3D","scatter3D",650,650);
  canInspection3D->Divide(2,2);
  gHitScatterPlot->SetMarkerStyle(20);
  canInspection3D->cd(1);
  gHitScatterPlot->Draw("y:x","hyperstrip>0","");
  gHitScatterPlot->GetHistogram()->SetTitle("GREEN(1) -> RED(2)");
  gHitScatterPlot->SetMarkerColor(kGreen);
  gHitScatterPlot->Draw("y:x","z==0 && hyperstrip>=1 && hyperstrip<=4","same"); 
  gHitScatterPlot->SetMarkerColor(kRed);
  gHitScatterPlot->Draw("y:x","z==0 && hyperstrip>=5 && hyperstrip<=8","same"); 
  gHitScatterPlot->SetMarkerStyle(47);
  gHitScatterPlot->SetMarkerColor(kMagenta); gHitScatterPlot->Draw("y:x","hyperstrip==-1","same");
  gHitScatterPlot->SetMarkerStyle(34);
  gHitScatterPlot->SetMarkerColor(kCyan); gHitScatterPlot->Draw("y:x","hyperstrip==0","same");
  

  canInspection3D->cd(2);
  gHitScatterPlot->SetMarkerStyle(1);
  gHitScatterPlot->Draw("y:x:z>>hisyxz","",""); /*>>hxyz(100,-40,40,100,-40,40,100,-40,40)*/ 
  gHitScatterPlot->SetMarkerStyle(7);
  gHitScatterPlot->SetMarkerColor(kBlue);
  gHitScatterPlot->Draw("y:x:z","z==0","same"); 
  gHitScatterPlot->SetMarkerColor(kGreen);
  gHitScatterPlot->Draw("y:x:z","z==0 && hyperstrip>=1 && hyperstrip<=4","same"); 
  gHitScatterPlot->SetMarkerColor(kRed);
  gHitScatterPlot->Draw("y:x:z","z==0 && hyperstrip>=5 && hyperstrip<=8","same"); 
  gHitScatterPlot->SetMarkerStyle(20);
  gHitScatterPlot->SetMarkerColor(kMagenta); 
  gHitScatterPlot->Draw("y:x:z","hyperstrip==-1","same");
  gHitScatterPlot->SetMarkerColor(kCyan); 
  gHitScatterPlot->Draw("y:x:z","hyperstrip==0","same");
  gHitScatterPlot->SetMarkerColor(kBlack);
  gHitScatterPlot->Draw("y:x:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");

  
  canInspection3D->cd(3);
  gHitScatterPlot->SetMarkerStyle(1);
  gHitScatterPlot->Draw("y:z>>hisyz","",""); 
  gHitScatterPlot->SetMarkerStyle(20);
  gHitScatterPlot->SetMarkerColor(kMagenta); gHitScatterPlot->Draw("y:z","hyperstrip==-1","same");
  gHitScatterPlot->SetMarkerColor(kCyan); gHitScatterPlot->Draw("y:z","hyperstrip==0","same");
  gHitScatterPlot->SetMarkerColor(kBlack);
  gHitScatterPlot->Draw("y:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");

  
  canInspection3D->cd(4);
  gHitScatterPlot->SetMarkerStyle(1);
  gHitScatterPlot->Draw("x:z>>hisyz","",""); 
  gHitScatterPlot->SetMarkerStyle(20);
  gHitScatterPlot->SetMarkerColor(kMagenta); gHitScatterPlot->Draw("x:z","hyperstrip==-1","same");
  gHitScatterPlot->SetMarkerColor(kCyan); gHitScatterPlot->Draw("x:z","hyperstrip==0","same");
  gHitScatterPlot->SetMarkerColor(kBlack);
  gHitScatterPlot->Draw("x:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");

  
  
}

/*****************************************************************************************************************/
int Minimise(void){

   ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  // ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");

	min->SetMaxFunctionCalls(1000000000);
	min->SetMaxIterations(1000000000);
	min->SetTolerance(0.001);

	ROOT::Math::Functor f(&GetChiSquare,3);
	double step[3] = {0.1,0.1,0.01};
	double variable[3] = {0,0,0}; // parameters: X,Y,Z of beam spot
	min->SetFunction(f);
	// Set the free variables to be minimized
  for(unsigned int i = 0 ; i < 3 ; i++)
	min->SetVariable(i,Form("Par%i",i),variable[i], step[i]);
	//min->SetVariableLimits(0,-5,+5);
	//min->FixVariable(0);

	gFinalMinimPosition.SetXYZ(-10, -10, -10);
  gFinalMinimPositionError.SetXYZ(-10, -10, -10);

  int result = min->Minimize();

  if (result){
  const double * minPos = min->X();
  const double * minPosErr = min->Errors();
    gFinalMinimPosition.SetXYZ(minPos[0],minPos[1],minPos[2]);
    gFinalMinimPositionError.SetXYZ(minPosErr[0],minPosErr[1],minPosErr[2]);
  }
  cout  << "RESULT \n";
   cout << "BeamSpot (mm)    : " << gFinalMinimPosition.X() << " " << gFinalMinimPosition.Y() << " " << gFinalMinimPosition.Z() << "\n" ;
   cout << "BeamSpotErr (mm) : " << gFinalMinimPositionError.X() << " " << gFinalMinimPositionError.Y() << " " << gFinalMinimPositionError.Z() << "\n" ;
   cout << "NB1: Above are the \"ABSOLUTE\" positions, calculated relative to the geometrical centroid of the barrel." <<endl;
   cout << "NB2: The geometrical centroid of the barrel is \"NOT\" necessary at the center of the experiment.\n" <<endl; 
  return result;
}

/*****************************************************************************************************************/ 
double GetChiSquare(const double parameters[]){

  double ChiSquare=0;
  double diff = 0 ;
  gAssumedAngle.clear();
  gAssumedPhiAngle.clear();
  gHitOnStrip.clear();
  gHitOnStripZ.clear();
  gHitOnStripY.clear();
  gHitOnStripX.clear(); 
  //create a beam spot and another point along the beam
  TVector3 BeamSpot(parameters[0],parameters[1],parameters[2]);
  // this is another point along the beam to define plane of action
  TVector3 AlongBeam(parameters[0],parameters[1],parameters[2]+1); // 1 is arbitrary
 
 // Construct the absolute positions of the central line along every strip
  for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
    //construct plane of action and get it's normal vector
    clsPlane3D*  ActionPlane = new clsPlane3D(gPointOnStrip[i], BeamSpot, AlongBeam); // order count
    TVector3 norm = ActionPlane->fNorm;
    //construct line of action
    TVector3 direction(0,0,1); // initiate direction of hit
    direction.Rotate(gAngle, norm); // rotate in the plane of action
    clsLine3D* ActionLine = new clsLine3D(BeamSpot, direction);
    //Construct the plane of the strip
    TVector3 NormalOnStrip = GetNormalOnDetector(gHyperStrip[i]); // strips in the same detector has same normal
    clsPlane3D* StripPlane = new clsPlane3D(gPointOnStrip[i], NormalOnStrip);
    //Get Intersection
    TVector3 HitPosition = GetLinePlaneIntersect(ActionLine,StripPlane);
    gHitOnStrip.push_back(HitPosition);
    gHitOnStripZ.push_back(HitPosition.Z());
    gHitOnStripY.push_back(HitPosition.Y());
    gHitOnStripX.push_back(HitPosition.X());
    //Calculate assumed angle (or the angle generated by the analysis)
    double angle = GetHitThetaAngle(HitPosition,gUserBeamSpot);
    gAssumedAngle.push_back(angle*TMath::RadToDeg());
    gAssumedPhiAngle.push_back(GetHitPhiAngle(HitPosition,gUserBeamSpot));
  } 
  
  for(unsigned int i=0; i<gHyperStrip.size(); i++){
    //cout << " Assumed " << gAssumedAngle[i] << endl; 
    diff = (gAssumedAngle[i]-gDataAngle[i]);
    //ChiSquare += diff*diff; // all have the same weight
    ChiSquare += diff*diff/(gDataAngleErr[i]*gDataAngleErr[i]);
  }
  ChiSquare = ChiSquare/(gHyperStrip.size()-3); // 3 is the parameters size 
  
  //cout << "   Chi2 is " << ChiSquare << endl; // used for testing
  
  return ChiSquare;
}

/*****************************************************************************************************************/
TVector3 GetNormalOnDetector(double hyperstrip){
  // get the detector number
  int det = (((int)hyperstrip-1)/4) + 1;
  //Calculate the hit position as if it hits detector 3 (at 12 o'clock i.e. perpendicular on the positive y-axis)
  TVector3 NormalOnDet(0,1,0);    // initiate with the normal on detector 3 at 12 o'clock (positive y-axis)
  NormalOnDet.RotateZ((3-det)*45*deg);
 
  return( NormalOnDet ) ;
}

/*****************************************************************************************************************/
TVector3 GetPointOnStrip(double hyperstrip, double pos){ // pos = [-1;+1]
  // get the detector and strip number
  int det   = (((int)hyperstrip-1)/4) + 1;
  int strip = (((int)hyperstrip-1)%4) + 1;
  //cout << " Extracting: " << det << " " << strip << " " << hyperstrip << endl ;
   
	// All in mm
  double INNERBARREL_PCB_Width  = 27.76;
  double INNERBARREL_ActiveWafer_Length = 94.80;
  double INNERBARREL_ActiveWafer_Width = 24.0;
  double StripPitch = INNERBARREL_ActiveWafer_Width/4.0;
  //Calculate the hit position as if it hits detector 3 (at 12 o'clock i.e. perpendicular on the positive y-axis)
  double Z = (0.5*INNERBARREL_ActiveWafer_Length) * (pos);
  double Y = INNERBARREL_PCB_Width*(0.5+sin(45*deg));
  double X = ((strip-2.5)*StripPitch);
  TVector3 aPos(X,Y,Z);        
  aPos.RotateZ((3-det)*45*deg);// looking downstream, Detector 1 is at 3 o'clock (negative x-axis)
 
  return( aPos ) ;
}


/*****************************************************************************************************************/
double GetHitThetaAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate angle from position assuming beam spot
	
	TVector3 zaxis(0,0,1);
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(zaxis) ) ;
}

/*****************************************************************************************************************/
double GetHitPhiAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate angle from position assuming beam spot
	
	TVector3 xaxis(1,0,0);
	TVector3 particle = HitPosition - BeamSpot;
	double angle = particle.Angle(xaxis)*TMath::RadToDeg();
	if (particle.Y()<0) angle = 360-angle;
	
  return( angle ) ;
}


/*****************************************************************************************************************/
TVector3 GetLinePlaneIntersect(clsLine3D* line, clsPlane3D* plane){
	TVector3 l = line->fDir;
	TVector3 l0 = line->fPoint;
	TVector3 n = plane->fNorm;
	TVector3 p0 = plane->fPoint;
	//solve the Line3D paramteric variable
	double t = ((p0-l0).Dot(n))/(l.Dot(n));
return (line->GetPoint(t));
}






