/*
This code is used to creat mock data of the inelastic channel hits 
on the barrel with a user-defined emission angle

input: 
1. Angle of emission (typically from an inelastic kinematic line calculation 
at the right beam energy)
2. x y z of the real beam spot
3. x y z of the assumed beam spot


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
#include "TVector3.h"
#include "TMath.h"
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
  vector <int> gBarrelNumber;
  vector <int> gStripNumber;
  vector <double> gHyperStrip;
  vector <double> gPhiStrip; // Phi of a strip
  vector <double> gNull;
  vector <double> gUnit;
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
 
  
//Functions
TVector3 GetLinePlaneIntersect(clsLine3D* l,clsPlane3D* p);
TVector3 GetHitOnStrip(TVector3 Beamspot, int StripNumber, double angle); // calculate impact position on strip middle line
double GetParticleAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate angle from position assuming beam spot
double GetParticlePhiAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate phi angle from position assuming beam spot
TVector3 GetNormalOnDetector(double hyperstrip);
TVector3 GetPointOnStrip(double hyperstrip, double pos=0); // pos = [-1;+1]

// MAIN
void CreateMockData(double Angle=-1, double sx=0, double sy=0, double sz=3, double ax=0, double ay=0, double az=0){

  if(Angle==-1) {
    cout << " ERROR ---- Provide Angle ! ----- "<< endl; 
    cout << " To execute:\n\n .x CreateMockData.C++( <angle-degree>,   <x>,<y>,<z> (real pos. mm),   <x>,<y>,<z> (assumed pos. mm)) \n"<< endl; 
    exit(-1); 
  }
  else {
    cout << " Angle input (deg.)                 : " << Angle << endl;
    cout << " Real Source position input (mm)    : " << sx << " " << sy << " " << sz << endl;
    cout << " Assumed Source position input (mm) : " << ax << " " << ay << " " << az << endl;
    }
  
  //Convert angle
  double AngleRad = Angle*TMath::DegToRad();
  //General Ntuple to stor values
  TNtuple *gHitScatterPlot = new TNtuple("gHitScatterPlot", "gHitScatterPlot", "x:y:z:hyperstrip");

// Create the strip and barrel vectors
  int barrel,strip;
   for (unsigned iBarrel=1 ; iBarrel <=8; iBarrel++ ){ // barrel
      for (unsigned iStrip=1 ; iStrip <=4; iStrip++ ){ // strip
      gBarrelNumber.push_back(iBarrel);
      gStripNumber.push_back(iStrip);
      gHyperStrip.push_back((iBarrel-1)*4+iStrip);
      gNull.push_back(0); // for plotting purposes
      gUnit.push_back(1); // for plotting purposes
      cout << " Filling " << iBarrel << " " << iStrip << " " << (iBarrel-1)*4+iStrip  << endl; 
      }
    }

    
// Construct the absolute positions (with respect to zero) of the central line along every strip
for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
  TVector3 temp = GetPointOnStrip(gHyperStrip[i],0);
  gPointOnStrip.push_back(temp);
  gPointOnStripY.push_back(temp.Y());
  gPointOnStripX.push_back(temp.X());
  gPointOnStripZ.push_back(temp.Z());
  gPhiStrip.push_back(temp.Phi()*TMath::RadToDeg());
  //Fill this point in the Ntuple 
  //gHitScatterPlot->Fill(temp.X(), temp.Y(), temp.Z(),gHyperStrip[i]);
  //inspect
  //cout << gHyperStrip[i] << "  " <<  temp.Phi()*TMath::RadToDeg()<< endl ; 
}
  
  
//Create mock data 
  
  //initiate
  gAssumedAngle.clear();
  gAssumedPhiAngle.clear();
  gHitOnStrip.clear();
  gHitOnStripZ.clear();
  gHitOnStripY.clear();
  gHitOnStripX.clear(); 
  
  //create a beam spot and another point along the beam
  TVector3 BeamSpot(sx,sy,sz);
  // this is another point along the beam to define plane of action
  TVector3 AlongBeam(sx,sy,sz+1); // 1 is arbitrary
 
 // Construct the absolute positions of the central line along every strip
  for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
    //construct plane of action and get it's normal vector
    clsPlane3D*  ActionPlane = new clsPlane3D(gPointOnStrip[i], BeamSpot, AlongBeam); // order count
    TVector3 norm = ActionPlane->fNorm;
    //construct line of action
    TVector3 direction(0,0,1); // initiate direction of hit
    direction.Rotate(AngleRad, norm); // rotate in the plane of action
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
    //Calculate assumed angle
    double angle = GetParticleAngle(HitPosition,TVector3(ax,ay,az));
    gAssumedAngle.push_back(angle*TMath::RadToDeg());
    angle = GetParticlePhiAngle(HitPosition,TVector3(ax,ay,az));
    gAssumedPhiAngle.push_back(angle*TMath::RadToDeg());
  } 
  
  //Fill in the 3D hist
  gHitScatterPlot->Fill(sx, sy, sz,0);
  for (unsigned i=0; i < gHitOnStrip.size() ; i++)
    gHitScatterPlot->Fill(gHitOnStrip[i].X(), gHitOnStrip[i].Y(), gHitOnStrip[i].Z(), gHyperStrip[i]);
    
  //Show results
  TCanvas* cantest= new TCanvas("cantest","cantest",650,650);
  cantest->Divide(3,3);
  
  TGraphErrors* grStripvsAssumedData = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gAssumedAngle[0],&gNull[0],&gNull[0]);
  grStripvsAssumedData->SetTitle(" (Apparent) Theta vs HyperStrip");
  grStripvsAssumedData->SetMarkerColor(2);
  grStripvsAssumedData->SetMarkerStyle(20);
  
  TGraph* grThetavsPhi = new TGraphErrors(gAssumedAngle.size(),&gAssumedPhiAngle[0],&gAssumedAngle[0]);
  grThetavsPhi->SetTitle(" (Apparent) Theta vs (Apparent) Phi");
  grThetavsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsPhi = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gPhiStrip[0]);
  grStripvsPhi->SetTitle(" Phi vs HyperStrip");
  grStripvsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsZ = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripZ[0]);
  grStripvsZ->SetTitle(" HitPosition(Z) vs HyperStrip");
  grStripvsZ->SetMarkerStyle(20);

  TGraph* grStripvsY = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripY[0]);
  grStripvsY->SetTitle(" HitPosition(Y) vs HyperStrip");
  grStripvsY->SetMarkerStyle(20);
  
  TGraph* grStripvsX = new TGraphErrors(gAssumedAngle.size(),&gHyperStrip[0],&gHitOnStripX[0]);
  grStripvsX->SetTitle(" HitPosition(X) vs HyperStrip");
  grStripvsX->SetMarkerStyle(20);
  
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
  
  cantest->cd(1);
  grStripvsAssumedData->Draw("alp"); 
  cantest->cd(2);
  grThetavsPhi->Draw("ap");
  cantest->cd(3);
  grStripvsPhi->Draw("ap");
  cantest->cd(4);
  grStripvsZ->Draw("alp");
  cantest->Draw();
  cantest->cd(5);
  grStripvsY->Draw("alp");
  cantest->cd(6);
  grStripvsX->Draw("alp");
  cantest->cd(7);
  grStripvsPointZ->Draw("alp");
  cantest->cd(8);
  grStripvsPointY->Draw("alp");
  cantest->cd(9);
  grStripvsPointX->Draw("alp");
  grStripvsPointY->Draw("same p");
  cantest->Draw();
    

  gHitScatterPlot->SetMarkerStyle(20);
  TCanvas* can3= new TCanvas("scatter3D","scatter3D",650,650);
  can3->Divide(2,3);
  can3->cd(1);
  gHitScatterPlot->Draw("y:x"); 
  can3->cd(2);
  gHitScatterPlot->Draw("y:x:z"); 
  can3->cd(3);
  gHitScatterPlot->Draw("x:hyperstrip"); 
  can3->cd(4);
  gHitScatterPlot->Draw("y:hyperstrip"); 
  can3->cd(5);
  gHitScatterPlot->Draw("y:z"); 
  can3->cd(6);
  gHitScatterPlot->Draw("x:z"); 
  
  can3->Draw();
  
  
  //Print out mock data 
  ofstream writedata;
  TString filename = Form("mockData_A%.1f_Real_%.1f_%.1f_%.1f_Assumed_%.1f_%.1f_%.1f.txt",Angle,sx,sy,sz,ax,ay,az);
  writedata.open(filename.Data());
  writedata << " 000   angle in degrees, position in mm " << endl  ;
  for (unsigned i = 0; i < gBarrelNumber.size() ; i++ ){
      cout << gBarrelNumber[i] << " " << gStripNumber[i] << " " << gAssumedAngle[i] << " " << gUnit[i] << endl  ;
      writedata << gBarrelNumber[i] << " " << gStripNumber[i] << " " << gAssumedAngle[i] << " " << gUnit[i] << endl  ;
      }
      
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
double GetParticleAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate theta angle (emission) from position assuming beam spot
	
	TVector3 zaxis(0,0,1);
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(zaxis) ) ;
}

/*****************************************************************************************************************/
double GetParticlePhiAngle(TVector3 HitPosition, TVector3 BeamSpot){ // Calculate angle from position assuming beam spot
	
	TVector3 xaxis(1,0,0);
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(xaxis) ) ;
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






