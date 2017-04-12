/// \file fitEuPeaks.cxx
/// \brief Routines to find the centroids of strong peaks in the 152Eu gamma source spectrum
///
#include <vector>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TRandom.h>

// ebennett additions (necessary?)
#include <TMath.h>
#include <TTree.h>
#include <TSpectrum.h>
#include <TMarker.h>
#include <fstream>

/// Random function
Double_t randomBin()
{
    return gRandom->Uniform(-0.5, 0.5);
}

/// Fit function - Gaussian Peak shape on top of polynomial background (up to 3rd order)
Double_t peakFitFunction(Double_t* x, Double_t* p)
{
	/**	\param x Length one array, variable
	 ** \param p Array of fit parameters:
	 ** p[0] = scale of gaussian
	 ** p[1] = centroid of gaussian
	 ** p[2] = sigma of gaussian
     ** p[3-6] = 0th - 3rd order parameters of polynomial
	 */
	
	return p[0]*TMath::Gaus(x[0], p[1], p[2]) + 
		p[3] + p[4]*x[0] + p[5]*TMath::Power(x[0], 2) + p[6]*TMath::Power(x[0], 3);
}


/// Perform peak fit in a limited range, return fit params
std::vector<Double_t> fitPeak(TH1* hist, Double_t nominalEnergy, const char* saveDir = 0)
{
	gStyle->SetOptFit(1);
	
	// MIN, MAX of fit region
	const Double_t MIN = nominalEnergy - 10;
	const Double_t MAX = nominalEnergy + 10;
  // Nominal starting width
	const Double_t peakWidth = 1.5;
	
	// Clone histogram to not modify the original
	TH1* hnew = (TH1*)hist->Clone("");

	// Set new histogram range to MIN, MAX
	hnew->GetXaxis()->SetRangeUser(MIN, MAX);

	// Get Peak height and center of new histogram (in range MIN, MAX)
	Double_t peakY = hnew->GetMaximum();
	Double_t peakX = hnew->GetBinCenter(hnew->GetMaximumBin());
  // Get nominal height of flat BG
	Double_t flatHeight = hnew->GetBinContent(hnew->FindBin(MIN));
	
	// Create TF1 w/ 7 parameters
	TF1 f("", peakFitFunction, MIN, MAX, 7);
	// Set initial parameters for gaussian and flat BG
	f.SetParameters(peakY-flatHeight, peakX, peakWidth, flatHeight);
	// Fix 1st-3rd order polynomial params to zero
	for (Int_t i=4; i< 7; ++i) { f.FixParameter(i, 0); }

	// Now perform the fit
	TFitResultPtr fitResult = hnew->Fit(&f, "qns", "", MIN, MAX);
	// Collect fit parameters to return
	std::vector<Double_t> fitpar(7);
	for(size_t i=0; i< 7; ++i) {
		fitpar[i] = fitResult.Get()->GetParams()[i];
	}

	// Now draw and save the fit result, if requested
	if(saveDir != 0) {
		TString fname = Form("%s/%s.png", saveDir, hist->GetName());
		TCanvas* csave = new TCanvas(Form("c_%s", hist->GetName()), hist->GetName());
//		hnew->SetMarkerStyle(20);
		hnew->Draw("");

		TF1 fgaus("","[0]*TMath::Gaus(x,[1],[2])", MIN, MAX);
		TF1 fbg("", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", MIN, MAX);
		for(int i=0; i< 3; ++i) { fgaus.SetParameter(i, f.GetParameter(i)); }
		for(int i=0; i< 4; ++i) { fbg.SetParameter(i, f.GetParameter(i+3)); }
		fgaus.SetLineColor(3);
		fbg.SetLineColor(6);

		fgaus.Draw("SAME");
		fbg.Draw("SAME");
		f.Draw("SAME");
		hnew->Draw("EHSAME");
		
		csave->SaveAs(fname);
		csave->Close();
	}

	// Now clean up memory items on the heap (pointers)
	hnew->Delete();

	return fitpar;
}

//// Create function that takes clover # and crystal # and generates fit parameters
//// Note:  assumes file is loaded
//// We're also going to ignore error handling, assuming this will be an internal functiond
std::vector<Double_t> coFit (TTree* physTree, int coClove, int coCryst)
{
    // Define variables (because you're not using python, Eames...)
    Double_t coSlop;
    Double_t coInpt;
    Double_t coLowP = 1173.2;
    Double_t coHigP = 1332.5;
    Double_t coPeak = coHigP - coLowP;

    // Create canvas
    TCanvas* coCanv = new TCanvas("coCanv","coCanv",600,1000);
    coCanv->Divide(1,2);

    // Create hist based on clover and crystal
    //      Make title and gating arguments
    TString tempTitl = Form("RAW 60Co Spectrum for Clover %d Crystal %d",coClove,coCryst);
    TString tempGate = Form("fGeTAMU_Core_CloverNbr_E==%d && fGeTAMU_Core_CrystalNbr_E==%d",coClove,coCryst);
    //      Init and fill hist
    coCanv->cd(1);
    TH1D* coTemp = new TH1D("coTemp",tempTitl,4096,0,4096);
    physTree->Project("coTemp","fGeTAMU_Core_Energy",tempGate,"");
    coTemp -> GetXaxis() -> SetRangeUser(0,4000); // this gets rid of overflow bin
    
    // Create TSpec and fit
    //      Init TSpec and find peaks
    TSpectrum coSpec;
    coSpec.Search(coTemp,2,"",0.1);
    //      Create ordered list of peaks
    //      Note:  ordered low to high
    Double_t* xPos = coSpec.GetPositionX();
    std::vector<int> coIndx(coSpec.GetNPeaks());
    TMath::Sort(coSpec.GetNPeaks(),xPos,&coIndx[0],false);
    
    // Calculate linear fit from peaks
    coSlop = coPeak/(xPos[coIndx.at(1)] - xPos[coIndx.at(0)]);
//    cout << coSlop << "\n";
    coInpt = coLowP - coSlop*xPos[coIndx.at(0)]; // sum math needed
//    cout << coInpt << "\n";
    
    // Create hist based on clover and crystal
    //      Make title and argument
    TString caliTitl = Form("Calibrated 60Co Spectrum for Clover %d Crystal %d",coClove,coCryst);
    TString caliArgu = Form("(fGeTAMU_Core_Energy+%f)*%f+%f",randomBin(),coSlop,coInpt);
    //      Make hist
    coCanv->cd(2);
    TH1D* coCali = new TH1D("coCali",caliTitl,4096,0,4096);
    physTree->Project("coCali",caliArgu,tempGate,"");
    coCali->Draw();
    //      Make markers (for visual check)
    TMarker m1(coLowP, 1000, 20);
    m1.Draw();
    TMarker m2(coHigP, 1000, 20);
    m2.Draw();
    //      Save Canvas
    TString fileTitl = Form("coFit/coFit_cl%d_cr%d.png",coClove,coCryst);
    coCanv->SaveAs(fileTitl);

    //  Create vector
    //      return [x^0, x^1, x^2, ...]
    std::vector<Double_t> coFitParam = {coInpt, coSlop};

    // Clean up
    coCanv->Close();
    coTemp->Delete();
    coCali->Delete();
    
    // Return vector
    return coFitParam;
}

void npFit (TTree* physTree)
{
    // Open file for writing
    std::ofstream outputFile;
    outputFile.open("GeCoCalParam.txt");
    // Create vectors to store params
    std::vector<vector<Double_t>> coCalParam;
    std::vector<TString> coCalTitle;
    // Loop over clovers
    for (int cloverNbr=1; cloverNbr<=4; cloverNbr++){
        // Loop over crystals
        for (int crystalNbr=1; crystalNbr<=4; crystalNbr++){
            // Call fit function
            coCalParam.push_back(coFit(physTree, cloverNbr, crystalNbr));
            // Create title
            TString tempToken = Form("GETAMU_D%d_CRY%d_E", cloverNbr, crystalNbr);
            coCalTitle.push_back(tempToken);
        }
    }
    // Loop over lists to output file
    // I know this is inefficient, but I want to stay flexible and have these lists ready for other uses in the future
    for (int totalToken=0; totalToken<static_cast<int>(coCalParam.size()); totalToken++){
        outputFile << coCalTitle.at(totalToken) << " " << coCalParam.at(totalToken).at(0) << " " << coCalParam.at(totalToken).at(1) << endl;
    }
    outputFile.close();
}
//  EAMES!  GET Eu POINTS FROM SHUYA RADWARE FIT

