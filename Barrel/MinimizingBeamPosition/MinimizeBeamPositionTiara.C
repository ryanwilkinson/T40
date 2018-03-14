/*
This code is used to minimise the position of the beam spot
from observing the inelastic channel in the barrel at a specific energy

input: 
1. Angle of emission 
(typically from an inelastic kinematic line calculation 
at the right beam energy)
2. text file of the observed
strip angle fwhm_angle

3. Set the reaction parameters (global variable)

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
#include "TEllipse.h"

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
#include "NPReaction.h"
#include "NPEnergyLoss.h"


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




//
//Global variables
//
  //Barrel
  double gAngle;
  vector <int> gBarrelNumber;
  vector <int> gStripNumber;
  vector <double> gHyperStrip;
  vector <double> gPhiOfStrip;// Phi of a strip [0,360]
  vector <double> gDataBarrel; // Z position
  vector <double> gDataBarrelErr;
  vector <double> gDataHitPhi;      //Calculated from gDataBarrel and assumed beam position
  vector <double> gNull;
  vector <double> gAssumedThetaAngle;
  vector <double> gAssumedPhiAngle;
  vector <TVector3> gZeroPointOnBarrelStrip;
  vector <TVector3> gHitOnStrip;
  vector <double> gHitOnStripZ;
  vector <double> gHitOnStripY;
  vector <double> gHitOnStripX;
  vector <double> gZeroPointOnBarrelStripZ;
  vector <double> gZeroPointOnBarrelStripY;
  vector <double> gZeroPointOnBarrelStripX;
  TNtuple *gMinimisationTree;
 
  //
  //Hyball
  //
  NPL::Reaction gReaction("19F(d,p)20F@189.805");  //@189.805 //187.264
  double gHyballZ = -147;     
  vector <double> gExcitation;
  vector <double> gRing;   //NB: rings go from [0 to 15]
  vector <double> ghyPhiData; // Local phi of the segments in the Hyball
  vector < vector <double> > ghyAssumedEnergy;  //proton energy in the hyball from assumed position
  vector < vector <double> > ghyEnergyData;
  vector < vector <double> > ghyEnergyDataErr;
  NPL::EnergyLoss ElossTarget;
  NPL::EnergyLoss ElossSilicon;
  double gTargetThick = 0*0.94*micrometer; // turn to zero in case the data are deadlayer corrected
  double gSiDeadLayerThick = 0*0.4*micrometer; // idem
      
  //fixed
  TVector3 gUserBeamSpot; 
  TVector3 gUserBeamDir; 
  double gUserV; //vertical deviation angle
  double gUserH; //horizontal deviation angle 
  
  //final result
	TVector3 gFinalMinimPosition;
  TVector3 gFinalMinimPositionError;
  double gFinalVBeamDeviat;
  double gFinalVBeamDeviatErr;
  double gFinalHBeamDeviat;
  double gFinalHBeamDeviatErr;
  TVector3 gFinalBeamDirection;

  //Functions
  TVector3 GetLinePlaneIntersect(clsLine3D* l,clsPlane3D* p);
  int Minimization(void); // numerical minimisation function which produces the final calibration parameters
  double GetChiSquareBarrelZ(const double parameters[]); // returns a chi squared value, using the Z of the hit 
  double GetChiSquareBarrelAngle(const double parameters[]); // returns a chi squared value, using the angle of the hit
  double GetChiSquareHyball(const double parameters[]); // returns a chi squared value 
  double GetChiSquare(const double parameters[]); // returns a chi squared value 
  TVector3 GetHitOnStrip(TVector3 Beamspot, int StripNumber, double angle); // calculate impact position on strip middle line
  double GetHitThetaAngle(TVector3 HitPosition, TVector3 BeamSpot, TVector3 BeamDir); // Calculate angle from position assuming beam spot and the
  double GetHitPhiAngle(TVector3 HitPosition, TVector3 BeamSpot); // Calculate phi angle [0,360] from position assuming beam spot
  TVector3 GetNormalOnDetector(double hyperstrip);
  TVector3 GetPointOnStrip(double hyperstrip, double pos=0); // pos = [-1;+1]
  double GetHitIndexOnLine(clsLine3D* line,TVector3 dir, TVector3 spot, double angle);
  double EnergyHyball(double localPhi, int iRing, double totalbeamE, double excit, double Z, TVector3 beamspot, TVector3 beamdir);



//
//
// MAIN
//
//
void MinimizeBeamPositionTiara(double Angle=-1, // angle in degree
  double ax=0, double ay=0, double az=0, //assumed beam position
  double bx=0, double by=0, double bz=1){ //assumed beam angle


//Input file
  TString BarrelData="Data/Barrel_000.txt";
  //TString BarrelData="Data/Barrel_shifted_000.txt"; //subtracting the Z shadow position 
  TString hyData    ="Data/Hyball_000_2states.txt";

  if(Angle==-1) {
    cout << " ERROR ---- Provide Angle ! ----- "<< endl; 
    cout << " To execute:\n\n .x MinimizeBeamPosition.C++( <angle-degree>, spotx,spoty,spotz,  dirx,diry,dirz  ) \n"<< endl; 
    exit(-1); 
  }
  else {
    cout << endl; 
    cout << "Angle provided: " << Angle << endl;
    cout << "File provided: " << BarrelData << endl;
    cout << "Assumed Source position input (mm) : " << ax << " " << ay << " " << az << endl;
    cout << "NB0: The Assumed Source position input \"MUST\" be the position of the beam spot used\nduring the BarrelData analysis \"AND\" calculated relative to the geometrical barrel centroid.\n" <<endl; 
    }
  

  TString name = BarrelData;
  name.ReplaceAll("Data/","");
  if(name.EndsWith(".txt")) name.ReplaceAll(".txt",".root");
  else name+=".root";
  TFile* Rootfile = new TFile(name,"RECREATE");
  
  //crunch the input data
  gAngle = Angle*TMath::DegToRad();
  gUserBeamSpot.SetXYZ(ax,ay,az);
  gUserBeamDir.SetXYZ(bx,by,bz);
  //Transform into aV and aH the deviation angles with respect 
  // to the horizontal (y=0) plane, and the vertical (x=0) plane
  gUserV=gUserBeamDir.Angle(TVector3(bx,0,bz));
  gUserH=gUserBeamDir.Angle(TVector3(0,by,bz));
  
  
// 
//
//HYBALL
//
//  
//
  //Excitation energies
  gExcitation.push_back(0.656*MeV);
  gExcitation.push_back(2.044*MeV);
  //Rings
  gRing.push_back(5); // the values will be averaged on three rings iring-1,iring,iring+1
  gRing.push_back(5);
  //target And silicon eenrgy losses     
  ElossTarget= NPL::EnergyLoss("proton_CD2.SRIM","SRIM",10);
  ElossSilicon= NPL::EnergyLoss("proton_Si.SRIM","SRIM",10);

  ghyEnergyData.resize(gExcitation.size());
  ghyEnergyDataErr.resize(gExcitation.size());

  vector <double> hyWedgeData;
  vector <double> hySectorData;
  vector <double> ghyNull;

  // Read the hyball data file and store values in c-vectors
  int hySector, hyWedge;
  double hyPhi;
  double hyE, hyErr;
  string line;

  ifstream hyfile(hyData.Data());
  if (hyfile.is_open()) {
    while ( hyfile>>hySector){
      if(hySector==0) {
        getline (hyfile,line);
        cout << line << endl;
        continue ; // skip "comments" line starting with zero
        }
      hyfile>>hyWedge;
      hyfile>>hyPhi;

      hySectorData.push_back(hySector);      
      hyWedgeData.push_back(hyWedge);
      if(hyPhi<0) hyPhi+=360;
      ghyPhiData.push_back(hyPhi);
      
      //printf("%d %d %f ",hySector,hyWedge, hyPhi) ;      
      for(unsigned i = 0 ; i < gExcitation.size() ; i++ ){
        hyfile>>hyE>>hyErr;
        ghyEnergyData[i].push_back(hyE);
        ghyEnergyDataErr[i].push_back(hyErr);
        //printf(" %f %f ",hyE,hyErr) ;      
        }
        //printf("\n");
        
      ghyNull.push_back(0); // for plotting purposes      
      }
    hyfile.close();
    }
  else { cout << "Unable to open file"; exit(-1);} 
 


// 
//
//Barrel
//
//  
// Read the data file of the barrel and store values in c-vectors
  int barrel,strip, hstrip;
  double bdata, bdataerr;
  ifstream myfile (BarrelData.Data());
  if (myfile.is_open()) {
    while ( myfile>>barrel){
    if(barrel==0) {
      getline (myfile,line);
      continue ; // skip "comments" line starting with zero
      }
    myfile>>strip>>bdata>>bdataerr;
    if(!(bdata>0)) continue ; 
      gBarrelNumber.push_back(barrel);
      gStripNumber.push_back(strip);
      gHyperStrip.push_back((barrel-1)*4+strip);
      gDataBarrel.push_back(bdata-4.84);
      gDataBarrelErr.push_back(bdataerr);
      gNull.push_back(0); // for plotting purposes
      //cout << " Reading " << barrel << " " << strip << " " << (barrel-1)*4+strip << " " <<  bdata << " " << bdataerr << endl; 
    }
    myfile.close();
  }
  else { cout << "Unable to open file"; exit(-1);} 
 
  TNtuple *BarrelHitScatterPlot = new TNtuple("BarrelHitScatterPlot", "BarrelHitScatterPlot", "x:y:z:hyperstrip");
  BarrelHitScatterPlot->Fill(gUserBeamSpot.X(), gUserBeamSpot.Y(), gUserBeamSpot.Z(),-1);
   
// Construct the absolute positions of the (z=0) point lying on the central line along every strip
for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
  TVector3 temp = GetPointOnStrip(gHyperStrip[i],0);
  gZeroPointOnBarrelStrip.push_back(temp);
  gZeroPointOnBarrelStripZ.push_back(temp.Z());
  gZeroPointOnBarrelStripY.push_back(temp.Y());
  gZeroPointOnBarrelStripX.push_back(temp.X()); 
  double phi_360 = temp.Phi()*TMath::RadToDeg();
  if(temp.Y()<0)  phi_360= 360+phi_360;  
  gPhiOfStrip.push_back(phi_360);
  //Fill this point in the Ntuple 
  BarrelHitScatterPlot->Fill(temp.X(), temp.Y(), temp.Z(),gHyperStrip[i]);
  //cout << gHyperStrip[i] << "  " <<  gDataBarrel.back()<< endl ; 
  }


 //
 //
 // Minimization
 //
 //

  gMinimisationTree = new TNtuple("gMinimisationTree", "gMinimisationTree", "x:y:z:v:h:chi2");  
  int result = 0 ;
  result = Minimization(); 
  cout << "Minimizing Result " << result << endl << endl;
  //Fill in the 3D after minimisation
  BarrelHitScatterPlot->Fill(gFinalMinimPosition.X(), gFinalMinimPosition.Y(), gFinalMinimPosition.Z(),0);
  for (unsigned i=0; i < gHitOnStrip.size() ; i++)
    BarrelHitScatterPlot->Fill(gHitOnStrip[i].X(), gHitOnStrip[i].Y(), gHitOnStrip[i].Z(), gHyperStrip[i]);
 
 //
 //
 // Show results
 //
 //
    
  TCanvas* canInspection= new TCanvas("canInspection","canInspection",650,650);
  canInspection->Divide(2,2);
  
  TGraphErrors* grStripvsData = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gDataBarrel[0],&gNull[0],&gDataBarrelErr[0]);
  grStripvsData->SetTitle("Z-Data vs HyperStrip");
  grStripvsData->SetMarkerColor(kBlack);
  grStripvsData->SetMarkerStyle(33);
  grStripvsData->SetMarkerSize(2.3);
  
  TGraphErrors* grStripvsAssumedData = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gAssumedThetaAngle[0],&gNull[0],&gNull[0]);
  grStripvsAssumedData->SetMarkerColor(kRed);
  grStripvsAssumedData->SetMarkerStyle(20);
  
  TGraph* grThetavsPhi = new TGraphErrors(gDataBarrel.size(),&gAssumedPhiAngle[0],&gAssumedThetaAngle[0]);
  grThetavsPhi->SetTitle("Theta (minimised) vs Phi (minimised)");
  grThetavsPhi->SetMarkerColor(kRed);
  grThetavsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsPhi = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gPhiOfStrip[0]);
  grStripvsPhi->SetTitle(" Phi (of Strip) vs HyperStrip");
  grStripvsPhi->SetMarkerStyle(20);
  
  TGraph* grStripvsZ = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gHitOnStripZ[0]);
  grStripvsZ->SetTitle(" HitPosition(Z) vs HyperStrip");
  grStripvsZ->SetMarkerStyle(20);
  grStripvsZ->SetMarkerColor(kRed);
  
  TGraphAsymmErrors* grStripvsDataZ = new TGraphAsymmErrors(gDataBarrel.size(),&gHyperStrip[0],&gDataBarrel[0],&gNull[0],&gNull[0],&gDataBarrelErr[0],&gDataBarrelErr[0]);
  grStripvsDataZ->SetTitle(" HitPosition(Z) (Data and minimized) vs HyperStrip");
  grStripvsDataZ->SetMarkerStyle(33);
  grStripvsDataZ->SetMarkerColor(kBlack);
  grStripvsDataZ->SetMarkerSize(2.3);

  TGraph* grStripvsMinimizedPhi = new TGraphErrors(gDataBarrel.size(),&gPhiOfStrip[0],&gAssumedPhiAngle[0]);
  grStripvsMinimizedPhi->SetTitle(" Phi (minimised) vs Phi (of Strip)");
  grStripvsMinimizedPhi->SetMarkerStyle(20);
  grStripvsMinimizedPhi->SetMarkerColor(kRed);

   TGraph* grHyperStripvsMinimizedPhi = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gAssumedPhiAngle[0]);
  grHyperStripvsMinimizedPhi->SetTitle(" Phi (minimised) vs Phi (of Strip)");
  grHyperStripvsMinimizedPhi->SetMarkerStyle(20);
  grHyperStripvsMinimizedPhi->SetMarkerColor(kRed);
   
  TGraph* grStripvsY = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gHitOnStripY[0]);
  grStripvsY->SetTitle(" HitPosition(Y) vs HyperStrip");
  grStripvsY->SetMarkerStyle(20);
  grStripvsY->SetMarkerColor(kRed);
  
  TGraph* grStripvsX = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gHitOnStripX[0]);
  grStripvsX->SetTitle(" HitPosition(X) vs HyperStrip");
  grStripvsX->SetMarkerStyle(20);
  grStripvsX->SetMarkerColor(kRed);
  
  TGraph* grStripvsPointZ = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gZeroPointOnBarrelStripZ[0]);
  grStripvsPointZ->SetTitle(" Point on Strip Position(Z) vs HyperStrip");
  grStripvsPointZ->SetMarkerStyle(28);
  grStripvsPointZ->SetMarkerColor(3);

  TGraph* grStripvsPointY = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gZeroPointOnBarrelStripY[0]);
  grStripvsPointY->SetTitle(" Point on Strip Position(Y) vs HyperStrip");
  grStripvsPointY->SetMarkerStyle(20);
  grStripvsPointY->SetMarkerColor(3);
  
  TGraph* grStripvsPointX = new TGraphErrors(gDataBarrel.size(),&gHyperStrip[0],&gZeroPointOnBarrelStripX[0]);
  grStripvsPointX->SetTitle(" Point on Strip Position(X) vs HyperStrip");
  grStripvsPointX->SetMarkerStyle(20);
  
  TGraphAsymmErrors* grDataZvsPhi = new TGraphAsymmErrors(gDataBarrel.size(),&gPhiOfStrip[0],&gDataBarrel[0],&gNull[0],&gNull[0],&gDataBarrelErr[0],&gDataBarrelErr[0]);
  grDataZvsPhi->SetTitle("HitPosition(Z) (Data and minimised) vs Phi (of Strip)");
  grDataZvsPhi->SetMarkerStyle(33);
  grDataZvsPhi->SetMarkerColor(kBlack);
  grDataZvsPhi->SetMarkerSize(2.3);
  
  TGraph* grMinimZvsPhi = new TGraphErrors(gDataBarrel.size(),&gPhiOfStrip[0],&gHitOnStripZ[0]);
  grMinimZvsPhi->SetTitle(" HitPosition(Z) vs HyperStrip");
  grMinimZvsPhi->SetMarkerStyle(20);
  grMinimZvsPhi->SetMarkerColor(kRed);
  
  canInspection->cd(1);
  grHyperStripvsMinimizedPhi->Draw("ap");
  canInspection->cd(2);
  grStripvsPhi->Draw("ap");
  canInspection->cd(3);
  grStripvsY->Draw("alp");
  canInspection->cd(4);
  grStripvsX->Draw("alp");
  
  //canInspection->cd(7);
  //grStripvsPointZ->Draw("alp");
  //canInspection->cd(8);
  //grStripvsPointY->Draw("alp");
  //canInspection->cd(9);
  //grStripvsPointX->Draw("alp");
  //grStripvsPointY->Draw("same p");
  
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
  BarrelHitScatterPlot->SetMarkerStyle(20);
  canInspection3D->cd(1);
  BarrelHitScatterPlot->Draw("y:x>>hisxy(160,-40,40,160,-40,40)","hyperstrip>0","");
  BarrelHitScatterPlot->GetHistogram()->SetTitle("GREEN(1) -> RED(2)");
  BarrelHitScatterPlot->SetMarkerColor(kGreen);
  BarrelHitScatterPlot->Draw("y:x","z==0 && hyperstrip>=1 && hyperstrip<=4","same"); 
  BarrelHitScatterPlot->SetMarkerColor(kRed);
  BarrelHitScatterPlot->Draw("y:x","z==0 && hyperstrip>=5 && hyperstrip<=8","same"); 
  BarrelHitScatterPlot->SetMarkerStyle(47);
  BarrelHitScatterPlot->SetMarkerColor(kMagenta); BarrelHitScatterPlot->Draw("y:x","hyperstrip==-1","same");
  BarrelHitScatterPlot->SetMarkerStyle(34);
  BarrelHitScatterPlot->SetMarkerColor(kCyan); BarrelHitScatterPlot->Draw("y:x","hyperstrip==0","same");
  

  canInspection3D->cd(2);
  BarrelHitScatterPlot->SetMarkerStyle(1);
  BarrelHitScatterPlot->Draw("y:x:z>>hisyxz","",""); //>>hxyz(100,-40,40,100,-40,40,100,-40,40) 
  BarrelHitScatterPlot->SetMarkerStyle(7);
  BarrelHitScatterPlot->SetMarkerColor(kBlue);
  BarrelHitScatterPlot->Draw("y:x:z","z==0","same"); 
  BarrelHitScatterPlot->SetMarkerColor(kGreen);
  BarrelHitScatterPlot->Draw("y:x:z","z==0 && hyperstrip>=1 && hyperstrip<=4","same"); 
  BarrelHitScatterPlot->SetMarkerColor(kRed);
  BarrelHitScatterPlot->Draw("y:x:z","z==0 && hyperstrip>=5 && hyperstrip<=8","same"); 
  BarrelHitScatterPlot->SetMarkerStyle(20);
  BarrelHitScatterPlot->SetMarkerColor(kMagenta); 
  BarrelHitScatterPlot->Draw("y:x:z","hyperstrip==-1","same");
  BarrelHitScatterPlot->SetMarkerColor(kCyan); 
  BarrelHitScatterPlot->Draw("y:x:z","hyperstrip==0","same");
  BarrelHitScatterPlot->SetMarkerColor(kBlack);
  BarrelHitScatterPlot->Draw("y:x:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");

  
  canInspection3D->cd(3);
  BarrelHitScatterPlot->SetMarkerStyle(1);
  BarrelHitScatterPlot->Draw("y:z>>hisyz","",""); 
  BarrelHitScatterPlot->SetMarkerStyle(20);
  BarrelHitScatterPlot->SetMarkerColor(kMagenta); BarrelHitScatterPlot->Draw("y:z","hyperstrip==-1","same");
  BarrelHitScatterPlot->SetMarkerColor(kCyan); BarrelHitScatterPlot->Draw("y:z","hyperstrip==0","same");
  BarrelHitScatterPlot->SetMarkerColor(kBlack);
  BarrelHitScatterPlot->Draw("y:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");

  
  canInspection3D->cd(4);
  BarrelHitScatterPlot->SetMarkerStyle(1);
  BarrelHitScatterPlot->Draw("x:z>>hisxz","",""); 
  BarrelHitScatterPlot->SetMarkerStyle(20);
  BarrelHitScatterPlot->SetMarkerColor(kMagenta); BarrelHitScatterPlot->Draw("x:z","hyperstrip==-1","same");
  BarrelHitScatterPlot->SetMarkerColor(kCyan); BarrelHitScatterPlot->Draw("x:z","hyperstrip==0","same");
  BarrelHitScatterPlot->SetMarkerColor(kBlack);
  BarrelHitScatterPlot->Draw("x:z","z>0 && z<94./2 && (hyperstrip!=-1 && hyperstrip!=0) ","same");



  //Hyball Graphs
  vector <TGraphErrors*> grHyExvsPhi;
  grHyExvsPhi.resize(gExcitation.size());
  for (unsigned i = 0 ; i < gExcitation.size() ; i++){
    grHyExvsPhi[i] = new TGraphErrors(ghyPhiData.size(),&ghyPhiData[0],&ghyEnergyData[i][0], &ghyNull[0], &ghyEnergyDataErr[i][0]);
    grHyExvsPhi[i]->SetTitle(Form("Ex %.2f vs hyPhi (data)",gExcitation[i]));
    grHyExvsPhi[i]->SetMarkerColor(kRed);
    grHyExvsPhi[i]->SetMarkerStyle(20);
    }
  
  vector <TGraph*> grHyExvsPhi_m;
  grHyExvsPhi_m.resize(gExcitation.size());
  for (unsigned i = 0 ; i < gExcitation.size() ; i++){  
    grHyExvsPhi_m[i] = new TGraph(ghyPhiData.size(), &ghyPhiData[0] ,&ghyAssumedEnergy[i][0]);
    grHyExvsPhi_m[i]->SetTitle(Form("Ex %.2f (minimised) vs hyPhi (data)",gExcitation[i]));
    grHyExvsPhi_m[i]->SetMarkerColor(kBlue);
    grHyExvsPhi_m[i]->SetMarkerStyle(20);
    }
  
  //plot 
  TCanvas* hyCan = new TCanvas("hyCan","Hyball",900,900);
  hyCan->Divide(2,2); 

/*for (unsigned i = 0 ; i < gExcitation.size() ; i++)
  if(i==0)
    grHyExvsPhi[i]->Draw("ap");
  else grHyExvsPhi[i]->Draw("p same");
  for (unsigned i = 0 ; i < gExcitation.size() ; i++) */
 
  hyCan->cd(1); 
  grHyExvsPhi[0]->GetYaxis()->SetRangeUser(ghyEnergyData[0][0]-0.1*MeV, ghyEnergyData[0][0]+0.1*MeV);
  grHyExvsPhi[0]->Draw("ap");
  grHyExvsPhi_m[0]->Draw("p same");

  hyCan->cd(3);
  grHyExvsPhi[1]->GetYaxis()->SetRangeUser(ghyEnergyData[1][0]-0.1*MeV, ghyEnergyData[1][0]+0.1*MeV);
  grHyExvsPhi[1]->Draw("ap");
  grHyExvsPhi_m[1]->Draw("p same");
   
  hyCan->cd(2); 
  for (unsigned i = 0 ; i < gExcitation.size() ; i++){  
    gReaction.SetExcitation4(gExcitation[i]);
    gReaction.GetKinematicLine3()->SetLineColor(i+1);
    if(i==0) gReaction.GetKinematicLine3()->Draw("alp");
    else gReaction.GetKinematicLine3()->Draw("lp same");
    }
     
  hyCan->cd(4);
  TGraph* grHyPolar = new TGraph((int)ghyPhiData.size());
  for(unsigned i = 0 ; i < ghyPhiData.size() ; i++)
    grHyPolar->SetPoint(i,cos(ghyPhiData[i]*TMath::DegToRad()),sin(ghyPhiData[i]*TMath::DegToRad()));
  grHyPolar->Set(ghyPhiData.size());  
  grHyPolar->SetTitle("hyPhi (data) in Unit circle");
  grHyPolar->SetMarkerColor(kRed);
  grHyPolar->SetMarkerStyle(20);
  grHyPolar->Draw("ap");
  TEllipse* unitcircle = new TEllipse(0,0,1,1);
  unitcircle->SetFillStyle(0);
  unitcircle->Draw("same");

  //
  //Wrap up 
  //  
  gMinimisationTree->Write();
  
  Rootfile->Write();
  Rootfile->Close();
  
//Sanity check, get intersection Beam wit hyball Plane
clsLine3D* aBeam       = new clsLine3D(gFinalMinimPosition,gFinalBeamDirection); 
clsPlane3D* HyballPlan  = new clsPlane3D(TVector3(0,0,-147), TVector3(0,0,1)); //
TVector3 Intersect = GetLinePlaneIntersect(aBeam,HyballPlan);
cout << " Intersect (beam, Hyball) : " << Intersect.X() << " "<< Intersect.Y() << " " << Intersect.Z() << ")\n" ;

  
}

/*****************************************************************************************************************/
int Minimization(void){

  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  //ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");

	min->SetMaxFunctionCalls(1000000000);
	min->SetMaxIterations(1000000000);
	min->SetTolerance(0.0001);

  unsigned dim = 9 ;
	ROOT::Math::Functor f(&GetChiSquare,dim);
  //parameters: X,Y,Z of beam spot, Vertical and horizontal deviations, Calibration offsets, HyballZ, beam energy
	//double variable[] = {0,0,0,0,0,0,0,gHyballZ,gReaction.GetBeamEnergy()}; 
	double variable[] = {2.9,4.5,-2.3,0,0,0,0,gHyballZ, gReaction.GetBeamEnergy()};  
  double step[]     = {0.1,0.1,0.1,0.01,0.01,0.1,0.1,0.1,0.1};

	min->SetFunction(f);
	// Set the free variables to be minimized
  for(unsigned int i = 0 ; i < dim ; i++)
	  min->SetVariable(i,Form("Par%i",i),variable[i], step[i]);
	
	min->SetVariableLimits(0,variable[0]-3, variable[0]+8); //X mm
	min->SetVariableLimits(1,variable[1]-3, variable[1]+8); //Y mm
	min->SetVariableLimits(2,variable[2]-8, variable[2]+8); //Z mm
	min->SetVariableLimits(3,variable[3]-2, variable[3]+2); //V degrees
	min->SetVariableLimits(4,variable[4]-2, variable[4]+2); //H degrees
  min->SetVariableLimits(5,variable[5]-2, variable[5]+2); //Calibration offset
	min->SetVariableLimits(6,variable[6]-2, variable[6]+2); //Calibration offset
  min->SetVariableLimits(7,variable[7]-3, variable[7]+3); // Hyball Z in mm
  min->SetVariableLimits(8,variable[8]-3, variable[8]+3); // BeamEnergy in MeV

//Uncommet to fix the variable value
	//min->FixVariable(0);  //X mm
	//min->FixVariable(1);  //Y mm
	//min->FixVariable(2);  //Z mm
  //min->FixVariable(3);  //V degrees
  //min->FixVariable(4);  //H degrees
	//min->FixVariable(5);    //Calib offset state 1
	//min->FixVariable(6);    //Calib offset state 2
  min->FixVariable(7);  //Hyball Z shift
  min->FixVariable(8);  //Beam Energy Shift

	gFinalMinimPosition.SetXYZ(-10, -10, -10);
  gFinalMinimPositionError.SetXYZ(-10, -10, -10);
	gFinalVBeamDeviat=-1;
  gFinalHBeamDeviat=-1;
	gFinalVBeamDeviatErr=-1;
  gFinalHBeamDeviatErr=-1;
  gFinalBeamDirection.SetXYZ(-10,-10,-10);

  int result = min->Minimize();

  const double * minPos;
  const double * minPosErr;
    
 // if (result){
    minPos = min->X();
    minPosErr = min->Errors();
    gFinalMinimPosition.SetXYZ(minPos[0],minPos[1],minPos[2]);
    gFinalMinimPositionError.SetXYZ(minPosErr[0],minPosErr[1],minPosErr[2]);
    gFinalVBeamDeviat=minPos[3];
    gFinalHBeamDeviat=minPos[4];
    gFinalVBeamDeviatErr=minPosErr[3];
    gFinalHBeamDeviatErr=minPosErr[4];
    gFinalBeamDirection.SetXYZ(
      cos(gFinalVBeamDeviat*TMath::DegToRad())*sin(gFinalHBeamDeviat*TMath::DegToRad()),
      sin(gFinalVBeamDeviat*TMath::DegToRad())*cos(gFinalHBeamDeviat*TMath::DegToRad()),
      cos(gFinalVBeamDeviat*TMath::DegToRad())*cos(gFinalHBeamDeviat*TMath::DegToRad()));
  //}
  
   cout  << "RESULT " << result << "  \n";
   cout << "BeamSpot (mm)       : " << gFinalMinimPosition.X() << " " << gFinalMinimPosition.Y() << " " << gFinalMinimPosition.Z() << "\n" ;
   cout << "BeamSpotErr (mm)    : " << gFinalMinimPositionError.X() << " " << gFinalMinimPositionError.Y() << " " << gFinalMinimPositionError.Z() << "\n\n" ;
   cout << "Beam V angle (deg)  : " << gFinalVBeamDeviat << " +/- " << gFinalVBeamDeviatErr << "\n" ;
   cout << "Beam H angle (deg)  : " << gFinalHBeamDeviat << " +/- " << gFinalHBeamDeviatErr << "\n" ;
   cout << " (=> Beam Direction : " << gFinalBeamDirection.X() << " "
                                    << gFinalBeamDirection.Y() << " " 
                                    << gFinalBeamDirection.Z() << ")\n\n" ;
   cout << "Shift Proton E1 (MeV)  : " << minPos[5] << " +/- " << minPosErr[5] << "\n" ;
   cout << "Shift Proton E2 (MeV)  : " << minPos[6] << " +/- "  << minPosErr[6] << "\n\n" ;
   cout << "Shift in Z (mm)        : " << minPos[7] << " +/- "  << minPosErr[7] << "\n" ;
   cout << "Beam Energy in (MeV)   : " << minPos[8] << " +/- "  << minPosErr[8] << "\n\n" ;
   cout << "NB1: Above are the \"ABSOLUTE\" positions, calculated relative to the geometrical centroid of the barrel." <<endl;
   cout << "NB2: The geometrical centroid of the barrel is \"NOT\" necessary at the center of the experiment.\n" <<endl; 
  return result;
}

/*****************************************************************************************************************/ 
double GetChiSquare(const double parameters[]){

  double chi2barrel = GetChiSquareBarrelZ(parameters); // beam spot and beam direction
  double chi2hyball = GetChiSquareHyball(parameters);// beam spot and beam direction
  //cout << " ChiSquare:  Hyball ; Barrel       " <<  chi2hyball << " " << chi2barrel 
  //     << "               ratio: " << chi2hyball/chi2barrel<< endl; 
  //Do some operations
  double chi2 = (1*chi2barrel) + (0.2*chi2hyball);

return chi2;

}




/*****************************************************************************************************************/
double GetChiSquareBarrelZ(const double parameters[]){ 

  double ChiSquare=0;
  double diff = 0 ;
  gAssumedPhiAngle.clear();
  gHitOnStrip.clear();
  gHitOnStripZ.clear(); // Z calculated from the Assumed position
  gHitOnStripY.clear();
  gHitOnStripX.clear(); 
  //create a beam spot and another point along the beam
  TVector3 BeamSpot(parameters[0],parameters[1],parameters[2]);
  //create the beam direction for this set of parameters
  double DeviatV = parameters[3]*TMath::DegToRad();
  double DeviatH = parameters[4]*TMath::DegToRad();
  TVector3 BeamDirection(cos(DeviatV)*sin(DeviatH), cos(DeviatH)*sin(DeviatV), cos(DeviatV)*cos(DeviatH)); //calculated at the center
    
 // Construct the absolute positions of the central line along every strip
  for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
    //Construct the backbone of the strip
    clsLine3D* StripCentralLine = new clsLine3D(gZeroPointOnBarrelStrip[i], TVector3(0,0,1)); // the Barrel is fixed and parallel to z-axis always
   
    //Get Intersection
    double t = GetHitIndexOnLine(StripCentralLine, BeamDirection, BeamSpot, gAngle);
    TVector3 HitPosition = StripCentralLine->GetPoint(t);
    gHitOnStrip.push_back(HitPosition);
    gHitOnStripZ.push_back(HitPosition.Z());
    gHitOnStripY.push_back(HitPosition.Y());
    gHitOnStripX.push_back(HitPosition.X());

    //Calculate assumed angle (or the angle generated by the user analysis)
    gAssumedThetaAngle.push_back(GetHitThetaAngle(HitPosition,gUserBeamSpot,gUserBeamDir));
    gAssumedPhiAngle.push_back(GetHitPhiAngle(HitPosition,gUserBeamSpot));
  } 
  
  for(unsigned int i=0; i<gHyperStrip.size(); i++){
    //cout << " Assumed " << gAssumedThetaAngle[i] << endl; 
    diff = (gHitOnStripZ[i]-gDataBarrel[i]);
    if(abs(diff)>10) cout << " Minimizer(Z) and Data(Z) " << gHitOnStripZ[i] << " - " << gDataBarrel[i] << " = " << diff << endl ; 
    ChiSquare += diff*diff/(gDataBarrelErr[i]*gDataBarrelErr[i]);
  }
  ChiSquare = ChiSquare/(gHyperStrip.size()-5); // mius the number of free parameters 
  
  gMinimisationTree->Fill(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],ChiSquare);
  //cout << "   Chi2 is " << ChiSquare << endl; // used for testing
  
  return ChiSquare;
}


/*****************************************************************************************************************/
double GetChiSquareBarrelAngle(const double parameters[]){

  double ChiSquare=0;
  double diff = 0 ;
  gAssumedThetaAngle.clear();
  gAssumedPhiAngle.clear();
  gHitOnStrip.clear();
  gHitOnStripZ.clear();
  gHitOnStripY.clear();
  gHitOnStripX.clear(); 
  //create a beam spot and another point along the beam
  TVector3 BeamSpot(parameters[0],parameters[1],parameters[2]);
  //create the beam direction for this set of parameters
  double DeviatV = parameters[3]*TMath::DegToRad();
  double DeviatH = parameters[4]*TMath::DegToRad();
  TVector3 BeamDirection(sin(DeviatH), sin(DeviatV), cos(DeviatV)*cos(DeviatH));
    
 // Construct the absolute positions of the central line along every strip
  for (unsigned i = 0 ; i < gHyperStrip.size() ; i++ ){
    //Construct the backbone of the strip
    clsLine3D* StripCentralLine = new clsLine3D(gZeroPointOnBarrelStrip[i], TVector3(0,0,1)); // the Barrel is fixed and parallel to z-axis always
   
    //Get Intersection
    double t = GetHitIndexOnLine(StripCentralLine, BeamDirection, BeamSpot, gAngle);
    TVector3 HitPosition = StripCentralLine->GetPoint(t);
    gHitOnStrip.push_back(HitPosition);
    gHitOnStripZ.push_back(HitPosition.Z());
    gHitOnStripY.push_back(HitPosition.Y());
    gHitOnStripX.push_back(HitPosition.X());
        
    //Calculate assumed angle (or the angle generated by the analysis)
    double angle = GetHitThetaAngle(HitPosition,gUserBeamSpot,gUserBeamDir);
    gAssumedThetaAngle.push_back(angle*TMath::RadToDeg());
    gAssumedPhiAngle.push_back(GetHitPhiAngle(HitPosition,gUserBeamSpot));
  } 
  
  for(unsigned int i=0; i<gHyperStrip.size(); i++){
    //cout << " Assumed " << gAssumedThetaAngle[i] << endl; 
    diff = (gAssumedThetaAngle[i]-gDataBarrel[i]);
    ChiSquare += diff*diff/(gDataBarrelErr[i]*gDataBarrelErr[i]);
  }
  ChiSquare = ChiSquare/(gHyperStrip.size()-5); // mius the number of free parameters 
  
  gMinimisationTree->Fill(parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],ChiSquare);
  //cout << "   Chi2 is " << ChiSquare << endl; // used for testing
  
  return ChiSquare;
}

/*****************************************************************************************************************/ 
double GetChiSquareHyball(const double parameters[]){

  double ChiSquare=0;
  double diff = 0 ;
  double AverageE = 0; //Average Energy on 3 rings
  
  ghyAssumedEnergy.clear();
  ghyAssumedEnergy.resize(gExcitation.size());

  //create a beam spot and another point along the beam
  TVector3 BeamSpot(parameters[0],parameters[1],parameters[2]);
  //create the beam direction for this set of parameters
  double DeviatV = parameters[3]*TMath::DegToRad();
  double DeviatH = parameters[4]*TMath::DegToRad();
  TVector3 BeamDirection(sin(DeviatH), sin(DeviatV), cos(DeviatV)*cos(DeviatH));
  
  //Extra hyball parameters
  double shift[2]= {parameters[5],parameters[6]};
  double hyZ= parameters[7]; 
  double beamEnergy = parameters[8];

    
  //Reference: EnergyHyball(double localPhi, int iRing, double excit=0, double Z=-147, TVector3 beamspot(0,0,0),TVector3 beamdir(0,0,1));
  for (unsigned i = 0 ; i < gExcitation.size() ; i++)
    for(unsigned j=0; j < ghyPhiData.size() ; j++){
      AverageE  = EnergyHyball(ghyPhiData[j], gRing[i]-1, beamEnergy, gExcitation[i], hyZ, BeamSpot,BeamDirection);
      AverageE += EnergyHyball(ghyPhiData[j], gRing[i]  , beamEnergy, gExcitation[i], hyZ, BeamSpot,BeamDirection);
      AverageE += EnergyHyball(ghyPhiData[j], gRing[i]+1, beamEnergy, gExcitation[i], hyZ, BeamSpot,BeamDirection);
      AverageE /=3.0;
      ghyAssumedEnergy[i].push_back(shift[i]+AverageE); //average of the three rings
    }
   for (unsigned i = 0 ; i < gExcitation.size() ; i++)
    for(unsigned j=0; j < ghyPhiData.size() ; j++){
      diff = (ghyAssumedEnergy[i][j] - ghyEnergyData[i][j]);
      ChiSquare += diff*diff/(ghyEnergyDataErr[i][j]*ghyEnergyDataErr[i][j]);
    }
  ChiSquare = ChiSquare/(gExcitation.size()*ghyPhiData.size()-6); //  minus the number of free parameters 
  //cout << "   Chi2 is " << ChiSquare << endl; // used for testing
  
  return ChiSquare;
}

/*****************************************************************************************************************/
TVector3 GetNormalOnDetector(double hyperstrip){
  // get the detector number
  int det = (((int)hyperstrip-1)/4) + 1;
  //Calculate the hit position as if it hits detector 3 (at 12 o'clock i.e. perpendicular on the positive y-axis)
  TVector3 NormalOnDet(0,1,0);    // initiate with the normal on detector 3 at 12 o'clock (positive y-axis)
  NormalOnDet.RotateZ((3-det)*45*TMath::DegToRad());
 
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
  double Y = INNERBARREL_PCB_Width*(0.5+sin(45*TMath::DegToRad()));
  double X = ((strip-2.5)*StripPitch);
  TVector3 aPos(X,Y,Z);        
  aPos.RotateZ((3-det)*45*TMath::DegToRad());// looking downstream, Detector 1 is at 3 o'clock (negative x-axis)
 
  return( aPos ) ;
}


/*****************************************************************************************************************/
double GetHitThetaAngle(TVector3 HitPosition, TVector3 BeamSpot, TVector3 BeamDir){ 
// Calculate angle from position assuming beam spot AND beam direction
	
	TVector3 particle = HitPosition - BeamSpot;
	
  return( particle.Angle(BeamDir) ) ;
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


double GetHitIndexOnLine(clsLine3D* line, TVector3 dir, TVector3 spot, double angle){
//This function calculates the index of the point from the a 3d-line
// forming 'angle' with a vector of direction 'dir' at vertex 'spot'
// Demonstration is simpl'ish'

  TVector3 linedir = line->fDir;
  TVector3 linepoint = line->fPoint;
   
  double A = dir.X() * (linepoint.X()-spot.X()) + dir.Y() * (linepoint.Y()-spot.Y());
  double B = pow((linepoint.X()-spot.X()),2) + pow((linepoint.Y()-spot.Y()),2);
  double bz = dir.Z();  
  //Solve equation alpha x^2() + beta x + gamma, where;
  double alpha = bz*bz-cos(angle)*cos(angle);
  double beta = 2*A*bz;
  double gamma = A*A-B*cos(angle)*cos(angle);
  
  /* 
  double A = dir.Dot(linepoint-spot);
  double B = (linepoint-spot).Mag2();
  double bz = dir.Z(); 
  double z0 = linepoint.Z();
  double zs = spot.Z();
  double c  = linedir.Z();
  //Solve equation alpha x^2() + beta x + gamma, where;
  double alpha = 1-(cos(angle)*cos(angle)/(bz*bz));
  double beta = 2*A - 2/bz*(z0-zs)*cos(angle)*cos(angle);
  double gamma = A*A-B*B*cos(angle)*cos(angle); 
   */
   
  double delta = beta*beta - 4*alpha*gamma;
  if (delta<0) {
   cout << " No real solutions, something is wrong!\n";
   exit(-1);
  }

  // we are interested in the positive solution downstream
  double z = (-beta + sqrt(delta))/(2*alpha);
  if (z<0) z = (-beta - sqrt(delta))/(2*alpha); 
  if (z<0) { // if it's still < than zero than there's a problem
    cout << " Both solutions are negative, something is wrong!\n";
    exit(-1);
  }

  //extract the index of this point onthe line
  double t = (z + spot.Z() - linepoint.Z())/linedir.Z();
  //double t = z/(bz*c);

//cout << t << endl; 
//cin.get();

  return t;

}


double EnergyHyball(double localPhi, int iRing, double totalbeamE, double excit, double Z, TVector3 beamspot,TVector3 beamdir){

//Returns hyball energy at a specific postion
//The position is calculated with respect to the beam spot
//The angle is calculated with respect to the beam direction
//Averages on 3 rings centered at 6 

//Get position
//rings go from [0-15]
  if(iRing < 0) return -1;
  if(iRing > 15) return -1;
  
  double r_min = 32.6;
  double r_max = 135.1;
  int NumberOfRings   = 16 ;   // 16
  int NumberOfSectors = 8 ; // 8
  double ring_pitch   = (r_max-r_min)/NumberOfRings  ;
  TVector3 StripCenter(0,0,0);
  StripCenter.SetX(r_min + (iRing+0.5)*ring_pitch); // build the detector at angle phi=0, then rotate
  StripCenter.SetY(0);
  StripCenter.SetZ(Z);
  StripCenter.RotateZ(localPhi*deg); //https://static.miraheze.org/t40wiki/5/55/TIARA_Detector_Map.png
  

  //Recalculate the position with respect to the current beam spot position
  TVector3 StripCenterAlt = StripCenter-beamspot;
  double angle = beamdir.Angle(StripCenterAlt) ; // angle with beam direction
  double ElossAngle = StripCenterAlt.Angle(TVector3(0,0,1)); //supposing Target and S1 orthogonal to (0,0,1)
  
  //Set excitation energy
  gReaction.SetBeamEnergy(totalbeamE); //MeV
  gReaction.SetExcitation4(excit); //MeV
  double E = gReaction.GetKinematicLine3()->Eval(angle/deg);
  E = ElossTarget.Slow( E ,gTargetThick/2., ElossAngle); // this angle is with respect to the target or detector normal
  E = ElossSilicon.Slow( E ,gSiDeadLayerThick, ElossAngle);
      
  return E;

}    

