#include <cstdio>
#include <fstream>
#include <TDirectory.h>
#include <TTree.h>
#include <TH1.h>

//
// Routine to convert ROOT histograms into Radware *.spe files
// Shamelessly stolen from the GRIFFIN Collaboration repository
// https://github.com/GRIFFINCollaboration/GRSISort
//
// Modified by GAC to make into a stand-alone package
//
// To use:
//    root.exe
//    .L Root2Rad.cxx+
//    Root2Rad(histo, "spectrum.spe");  // create from existing ROOT TH1
//    // or //
//    Root2Rad(T40Tree, "EGamma", "EGamma > 500", 4096, 0, 4096, "spectrum.spe");
//    // ^^ Make a new histogram from a TTree, then write it ^^ //

struct SpeHeader {
	Int_t buffsize;          /*fortran file, each record starts with record size */            // 14
	Char_t label[8];
	Int_t size;
	Int_t junk1;
	Int_t junk2;
	Int_t junk3;
	Int_t buffcheck;         /*fortran file, record ends with record size :) */                // 14
} __attribute__((packed));


// Write from existing histogram //
void Root2Rad(TH1 *hist, const char* output_file_name){
	std::ofstream outfile(output_file_name);
	SpeHeader spehead;
	spehead.buffsize = 24;
	strncpy(spehead.label,hist->GetName(),8);

	if(hist->GetRMS() > 16384/2) {
		while(hist->GetNbinsX()>16384){
			hist = hist->Rebin(2);
			fprintf(stderr, "\t!!  WARNING: %s has been compressed by 2.\n",hist->GetName());
		}
		spehead.size = hist->GetNbinsX();
	} else if(hist->GetNbinsX()>16384) {
		spehead.size = 16384;
	} else {
		spehead.size = hist->GetNbinsX();
	}

	spehead.junk1 = 1;
	spehead.junk2 = 1;
	spehead.junk3 = 1;
	spehead.buffcheck = 24;         /*fortran file, record ends with record size :) */                // 14

	outfile.write((Char_t*)(&spehead),sizeof(SpeHeader));

	Int_t histsizeinbytes = spehead.size *4;

	outfile.write((Char_t*)&histsizeinbytes,sizeof(Int_t));
	Float_t bin = 0.0;
	for(int x=1;x<=spehead.size;x++){
		if(x<=hist->GetNbinsX()){
			bin = (Float_t)hist->GetBinContent(x);
			outfile.write((Char_t*)&bin,sizeof(Int_t));
		}
		else {
			bin = 0.0;
			outfile.write((Char_t*)&bin,sizeof(Int_t));
		}
	}

	outfile.write((Char_t*)&histsizeinbytes,sizeof(Int_t));

	return;
}

// Write directly from TTree //
void Root2Rad(TTree* t, const char* param, const char* cut, Int_t bins, Float_t low, Float_t high, const char* output_file_name)
{
	t->Draw(Form("%s>>RADWARE_DUMMY_HIST(%i, %f, %f)", param, bins, low, high), cut, "goff");
	TH1 *hist = dynamic_cast<TH1*>(gDirectory->Get("RADWARE_DUMMY_HIST"));
	if(!hist) {
		fprintf(stderr,"\t!!  ERROR: Couldn't create histogram.\n");
		return;
	}
	Root2Rad(hist, output_file_name);
	hist->Delete();
}
