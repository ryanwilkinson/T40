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
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
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
	clsLine3D(TVector3 lp,TVector3 l){
		fDir = l; 
		fPoint = lp;
		}
	~clsLine3D();

	TVector3 GetPoint(double t){ 
		return (t*fDir+fPoint);
		}
};


class clsPlane3D{
  public:
	//members
  TVector3 fPoint; // point on the line
	TVector3 fNorm; // directing vector
	public:
	//functions
	clsPlane3D(TVector3 np, TVector3 n){
	  fNorm = n;  //normal vector
		fPoint = np; // point on plane
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
  vector <double> gPhiStrip; // Phi of a strip
  vector <double> gDataAngle;
  vector <double> gDataAngleSigma;
  vector <double> gNull;
  vector <double> gAssumedAngle;
  vector <double> gAssumedPhiAngle;
  vector <TVector3> gPointOnStrip;
  vector <TVector3> gHitOnStrip;
  vector <double> gHitOnStripZ;
  vector <TVector3> gNormalOnDetector;
  //recycled with every iteration
  //...
  //finale result
	TVector3 gFinalMinimPosition;
  TVector3 gFinalMinimPositionError;
  
//Functions
TVector3 GetLinePlaneIntersect(clsLine3D* l,clsPlane3D* p);
int Minimise(void); // numerical minimisation function which produces the final calibration parameters
double GetChiSquare(const double parameters[]); // returns a chi squared value 
TVector3 GetHitOnStrip(TVector3 Beamspot, int StripNumber, double angle); // calculate impact position on strip middle line
double GetParticleAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate angle from position assuming beam spot
double GetParticlePhiAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate angle from position assuming beam spot
TVector3 GetNormalOnDetector(double hyperstrip);
TVector3 GetPointOnStrip(double hyperstrip, double pos=0); // pos = [-1;+1]

// MAIN
void MinimizeBeamPosition(double Angle=-1, TString data="BarrelAngles_Mg25.txt"){

  if(Angle==-1) {
    cout << " ERROR ---- Provide Angle ! ----- "<< endl; 
    cout << "execute:\n .x MinimizeBeamPosition.C++( <angle-degree>, <input-ascii-file> ) "<< endl; 
    exit(-1); 
  }
  else {
    cout << " Angle provided " << Angle << endl;
    }
  
  gAngle = Angle*TMath::DegToRad();
// Read the data file and store values in c-vectors
  int barrel,strip, hstrip;
  double angle, sangle;
  string line;
  ifstream myfile (data.Data());
  if (myfile.is_open()) {
    while ( myfile>>barrel){
    if(barrel==0) {
      getline (myfile,line);
      continue ; // skip comments line starting with zero
      }
    myfile>>strip>>angle>>sangle;
    if(!(angle>0)) continue ; 
      gBarrelNumber.push_back(barrel);
      gStripNumber.push_back(strip);
      gHyperStrip.push_back((barrel-1)*4+strip);
      gDataAngle.push_back(angle);
      gDataAngleSigma.push_back(sangle);
      gNull.push_back(0); // for plotting purposes
    }
    myfile.close();
  }
  else { cout << "Unable to open file"; exit(-1);} 
  
// Construct the absolute positions of the central line along every strip
for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
  TVector3 temp = GetPointOnStrip(gHyperStrip[i],0);
  gPointOnStrip.push_back(temp);
  gPhiStrip.push_back(temp.Phi()*TMath::RadToDeg());
  //cout << gHyperStrip[i] << "  " <<  temp.Phi()*TMath::RadToDeg()<< endl ; 
}

// Construct the plane of every strip (1 plane for every 4 strip)
for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
  gNormalOnDetector.push_back(GetNormalOnDetector(gHyperStrip[i]));
}


// perform the numerical minimisation
  int result = 0 ;
  result = Minimise(); 
  cout << "Minimizing Result " << result << endl;
  
  //Show results
  TCanvas* cantest= new TCanvas("cantest","cantest",650,650);
  cantest->Divide(2,2);
  
  TGraphErrors* grStripvsData = new TGraphErrors(gDataAngle.size(),&gHyperStrip[0],&gDataAngle[0],&gNull[0],&gDataAngleSigma[0]);
  grStripvsData->SetTitle(" minTheta vs HyperStrip");
  grStripvsData->SetMarkerStyle(24);
  
  TGraphErrors* grStripvsAssumedData = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gAssumedAngle[0],&gNull[0],&gNull[0]);
  grStripvsAssumedData->SetMarkerColor(2);
  grStripvsAssumedData->SetMarkerStyle(20);
  
  TGraph* grThetavsPhi = new TGraphErrors(gAssumedAngle.size(),&gAssumedPhiAngle[0],&gAssumedAngle[0]);
  grThetavsPhi->SetTitle(" minTheta vs minPhi");
  grThetavsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsPhi = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPhiStrip[0]);
  grStripvsPhi->SetTitle(" Phi vs HyperStrip");
  grStripvsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsZ = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripZ[0]);
  grStripvsZ->SetTitle(" HitPosition(Z) vs HyperStrip");
  grStripvsZ->SetMarkerStyle(20);

  cantest->cd(1);
  grStripvsData->Draw("alp"); 
  grStripvsAssumedData->Draw("p same"); 
  cantest->cd(2);
  grThetavsPhi->Draw("ap");
  cantest->cd(3);
  grStripvsPhi->Draw("ap");
  cantest->cd(4);
  grStripvsZ->Draw("alp");
  cantest->Draw();

  TCanvas* can2= new TCanvas("can2","can2",650,650);
  grStripvsData->Draw("alp"); 
  grStripvsAssumedData->Draw("p same"); 
  can2->Draw();

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
   cout << "BeamSpot:    " << gFinalMinimPosition.X() << " " << gFinalMinimPosition.Y() << " " << gFinalMinimPosition.Z() << "\n" ;
   cout << "BeamSpotErr: " << gFinalMinimPositionError.X() << " " << gFinalMinimPositionError.Y() << " " << gFinalMinimPositionError.Z() << "\n" ;
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
  
  //create a beam spot and another point along the beam
  TVector3 BeamSpot(parameters[0],parameters[1],parameters[2]);
  // this is another point along the beam to define plane of action, 1 is arbitrary
  TVector3 AlongBeam(parameters[0],parameters[1],parameters[2]+1);
 
 // Construct the absolute positions of the central line along every strip
  for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
    //construct plane of action and get it's normal vector
    clsPlane3D*  ActionPlane = new clsPlane3D(gPointOnStrip[i], BeamSpot, AlongBeam);
    TVector3 norm = ActionPlane->fNorm;
    //construct line of action
    TVector3 direction(0,0,1);
    direction.Rotate(gAngle, norm);
    clsLine3D* ActionLine = new clsLine3D(BeamSpot, direction);
    //Construc the plane of the strip
    TVector3 NormalOnStrip = GetNormalOnDetector(gHyperStrip[i]); 
    clsPlane3D* StripPlane = new clsPlane3D(gPointOnStrip[i], gNormalOnDetector[i]);
    //Get Intersection
    TVector3 HitPosition = GetLinePlaneIntersect(ActionLine,StripPlane);
    gHitOnStrip.push_back(HitPosition);
    gHitOnStripZ.push_back(HitPosition.Z());
    //Calculate assumed angle
    double angle = GetParticleAngle(HitPosition,TVector3(0,0,0));
    gAssumedAngle.push_back(angle*TMath::RadToDeg());
    angle = GetParticlePhiAngle(HitPosition,TVector3(0,0,0));
    gAssumedPhiAngle.push_back(angle*TMath::RadToDeg());
  } 
  
  for(unsigned int i=0; i<gHyperStrip.size(); i++){
    //cout << " Assumed " << gAssumedAngle[i] << endl; 
    diff = (gAssumedAngle[i]-gDataAngle[i]);
    //ChiSquare += diff*diff; // all have the same weight
    ChiSquare += diff*diff/(gDataAngleSigma[i]*gDataAngleSigma[i]);
  }
  ChiSquare = ChiSquare/(gHyperStrip.size()-3); // 3 is the parameters size 
  
  cout << "   Chi2 is " << ChiSquare << endl; // used for testing
  
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
	// All in mm
  double INNERBARREL_PCB_Width  = 27.76;
  double INNERBARREL_ActiveWafer_Length = 94.80;
  double INNERBARREL_ActiveWafer_Width = 24.0;
  double StripPitch = INNERBARREL_ActiveWafer_Width/4.0;
  //Calculate the hit position as if it hits detector 3 (at 12 o'clock i.e. perpendiculart on the positive y-axis)
  double Z = (0.5*INNERBARREL_ActiveWafer_Length) * (pos);
  double Y = INNERBARREL_PCB_Width*(0.5+sin(45*deg));
  double X = (strip*StripPitch-0.5*INNERBARREL_ActiveWafer_Width)-(0.5*StripPitch);
  TVector3 aPos(X,Y,Z);        // since RowPos = (U-D)/(U+D) => Downstream hit (i.e. Z>0) has RowPos<0, thus the sign
  aPos.RotateZ((3-det)*45*deg);// looking downstream Detector 1 is at 3 o'clock (negative x-axis)
 
  return( aPos ) ;
}


/*****************************************************************************************************************/
double GetParticleAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate angle from position assuming beam spot
	
	TVector3 zaxis(0,0,1);
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(zaxis) ) ;
}

/*****************************************************************************************************************/
double GetParticlePhiAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate angle from position assuming beam spot
	
	TVector3 zaxis(1,0,0);
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(zaxis) ) ;
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






