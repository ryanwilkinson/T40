// This macro is used for offline analysis to load files/scripts, load cuts and set relevant aliases.
// Run this script in an 'nptoolroot' session
#include "TCutG.h"
#include "TFile.h"
#include <iostream>

void nptoolrootsetup(){

/// load up file to work on

// load up relevent macros within the 'nptool/Projects/T40/' subdirectory
  //gROOT->ProcessLine(".L /home/rw00227/nptool/Projects/T40/");

// defining aliases for the plotting of certain quantities (such as RF TOF)
  cout << "\nSetting aliases..." << endl;
  cout << "...done!" << endl;

// loading cuts
  cout << "\nLoading cuts..." << endl;
  TFile* EdEcutFile = new TFile("/home/rw00227/nptool/Projects/T40/cuts/EdEcut.root");
  TCutG* EdEcut = (TCutG*) EdEcutFile->Get("EdEcut");
  TFile* AwcutFile = new TFile("/home/rw00227/nptool/Projects/T40/cuts/Awcut.root");
  TCutG* Awcut = (TCutG*) AwcutFile->Get("Awcut");
  cout << "...done!" << endl;

// set decent colour PlastLeftTime
  gStyle->SetPalette(57);
}
