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
#include <TROOT.h>
#include <TGraph.h>

// ebennett additions (necessary?)
#include <TMath.h>
#include <TTree.h>
#include <TSpectrum.h>
#include <TMarker.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

/// Random function
/*Double_t randomBin()
{
    return gRandom->Uniform(-0.5, 0.5);
}
*/

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
	// EAB:  Added 170420 to deal with empty fitpar bug
    if (peakY == 0){
        hnew->Delete();
        std::vector<Double_t> fitpar = {0, 0, 0, 0, 0, 0, 0};
        return fitpar;
    }
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
    // Make sure random file is loaded
    if (gROOT->GetClass("randomBinDummy") == 0)
        gROOT->ProcessLine(".L randomBin.cxx+");

    // Define variables (because you're not using python, Eames...)
    Double_t coSlop;
    Double_t coInpt;
    Double_t coLowP = 1173.238;
    Double_t coHigP = 1332.502;
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
    coSpec.Search(coTemp,2,"",0.3);
    //      Create ordered list of peaks
    //      Note:  ordered low to high
    Double_t* xPos = coSpec.GetPositionX();
    std::vector<int> coIndx(coSpec.GetNPeaks());
    TMath::Sort(coSpec.GetNPeaks(),xPos,&coIndx[0],false);
    
    // Calculate linear fit from peaks
    coSlop = coPeak/(xPos[coIndx.at(1)] - xPos[coIndx.at(0)]);
    coInpt = coLowP - coSlop*xPos[coIndx.at(0)]; // sum math needed
    
    // Create hist based on clover and crystal
    //      Make title and argument
    TString caliTitl = Form("Calibrated 60Co Spectrum for Clover %d Crystal %d",coClove,coCryst);
    TString caliArgu = Form("(fGeTAMU_Core_Energy+randomBin())*%f+%f",coSlop,coInpt);
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

//// This is a preliminary function to generate calibration files for nptool
//// v1 -- 170329
int npFit (TTree* physTree)
{
    // Open file for writing
    std::ofstream outputFile;
    outputFile.open("coFit/GeCoCalParam.txt");
    // Create vectors to store params
    std::vector<std::vector<Double_t>> coCalParam;
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
        outputFile << coCalTitle.at(totalToken) << " " << coCalParam.at(totalToken).at(0) << " " << coCalParam.at(totalToken).at(1) << std::endl;
    }
    outputFile.close();
    return 0;
}

//// Rough program for reading data out of nptool calibration file
std::vector<std::string> getFitParam (std::string fitParamFile, int cloverNbr, int crystalNbr)
{
    // create token
    std::string tempToken = "GETAMU_D"+ std::to_string(cloverNbr)+"_CRY"+ std::to_string(crystalNbr) +"_E";
    
    // Start by reading in the correct fit parameters from file
    //      We want a string var for the line we read in and a buffer for the string stream
    std::string line;
    std::string buff;
    std::vector<std::string> lineData;
    //      Now we open the file for reading
    std::ifstream coParamFile (fitParamFile);
    if (coParamFile.is_open()){
        // This is going to iterate over each line in the file
        while (std::getline(coParamFile,line)){
            // Insert our line into a string stream and then push to a vector
            std::stringstream ss(line);
            while (ss >> buff)
                lineData.push_back(buff);
            // Check to see if token is what we're looking for
            if(lineData.at(0) == tempToken)
                break; // If so end the loop
            else
                lineData.clear(); // If not, clear vector and try next line
        }
        // Close our file and return
        coParamFile.close();
        return lineData;
    }
    else{
        // In the case we can't open the file, we'll throw this error (can't read as sub-function, though)
        std::cout << "Unable to find calibration parameters." << std::endl;
        return lineData; // This should just return an empty list...  Need better error handling here
    }
}

//// fitEuPeaks finally!  This *should* take the calibration parameters from coFit and then find peaks in Eu
int fitEuPeak (TTree* physTree, std::string fitFile, int cloverNbr, int crystalNbr, std::vector<Double_t>& retFit, std::vector<Double_t>& retChn)
{
    // Make sure random file is loaded
    if (gROOT->GetClass("randomBinDummy") == 0)
        gROOT->ProcessLine(".L randomBin.cxx+");

    // Start by defining the Eu peaks we want to fit (thanks, Shuya)
    Double_t euP1 = 121.7824;
    Double_t euP2 = 244.6989;
    Double_t euP3 = 344.2811;
    Double_t euP4 = 778.903;
    Double_t euP5 = 964.055;
    Double_t euP6 = 1112.087;
    Double_t euP7 = 1408.022;
    std::vector<Double_t> euPeaks = {euP1,euP2,euP3,euP4,euP5,euP6,euP7};
    retFit.clear();
    retChn.clear();

    // Read in parameters 
    std::vector<std::string> fitParam = getFitParam(fitFile, cloverNbr, crystalNbr);
    if (fitParam.empty())
        return 1;
    // Convert strings to doubles
    Double_t coSlop = std::stod(fitParam.at(2));
    Double_t coInpt = std::stod(fitParam.at(1));
    
    // Now we can create a histogram with the calibrated Eu data
    TCanvas* euCanv = new TCanvas("euCanv","euCanv",600,600);
    //      Make title and gating arguments
    TString euTitl = Form("152Eu Spectrum for Clover %d Crystal %d",cloverNbr,crystalNbr);
    TString euArgu = Form("(fGeTAMU_Core_Energy+randomBin())*%f+%f",coSlop,coInpt);
    TString euGate = Form("fGeTAMU_Core_CloverNbr_E==%d && fGeTAMU_Core_CrystalNbr_E==%d",cloverNbr,crystalNbr);
    //      Init and fill hist
    TH1D* euHist = new TH1D("euHist", euTitl,4096,0,4096);
    physTree->Project("euHist",euArgu,euGate,"");
    euHist->Draw();
    //      Save hist
    TString fileTitl = Form("euFit/euFit_cl%d_cr%d.png",cloverNbr,crystalNbr);
    euCanv->SaveAs(fileTitl);


    // Okay, so now let's call Greg's fitting function
    //      We'll have to loop over each peak
    for(size_t i=0; i<7; ++i){
        Double_t peakValue = fitPeak(euHist, euPeaks.at(i), 0).at(1);
        if (peakValue == 0){
            continue;
        }
        else{
            retFit.push_back(euPeaks.at(i));
            // Back-convert to channel number
            Double_t chanValue = (peakValue - coInpt)/coSlop;
            retChn.push_back(chanValue);
        }
    }

    // Clean up our mess
    euHist->Delete();
    euCanv->Close();

    return 0;
}

//// script to find peaks in background data
int fitBkPeak(TTree* physTree, std::string fitFile, int cloverNbr, int crystalNbr, std::vector<Double_t>& retFit, std::vector<Double_t>& retChn)
{
    // Make sure random file is loaded
    if (gROOT->GetClass("randomBinDummy") == 0)
        gROOT->ProcessLine(".L randomBin.cxx+");

    // Start by defining the Eu peaks we want to fit (thanks, Shuya)
    Double_t bkP1 = 1764.515;
    Double_t bkP2 = 2614.522;
    std::vector<Double_t> bkPeaks = {bkP1,bkP2};

    // Read in parameters 
    std::vector<std::string> fitParam = getFitParam(fitFile, cloverNbr, crystalNbr);
    if (fitParam.empty())
        return 1;
    // Convert strings to doubles
    Double_t coSlop = std::stod(fitParam.at(2));
    Double_t coInpt = std::stod(fitParam.at(1));
    
    // Now we can create a histogram with the calibrated Eu data
    TCanvas* bkCanv = new TCanvas("bkCanv","bkCanv",600,600);
    //      Make title and gating arguments
    TString bkTitl = Form("Background Spectrum for Clover %d Crystal %d",cloverNbr,crystalNbr);
    TString bkArgu = Form("(fGeTAMU_Core_Energy+randomBin())*%f+%f",coSlop,coInpt);
    TString bkGate = Form("fGeTAMU_Core_CloverNbr_E==%d && fGeTAMU_Core_CrystalNbr_E==%d",cloverNbr,crystalNbr);
    //      Init and fill hist
    TH1D* bkHist = new TH1D("bkHist", bkTitl,4096,0,4096);
    physTree->Project("bkHist",bkArgu,bkGate,"");
    bkHist->Draw();

    // Okay, so now let's call Greg's fitting function
    //      We'll have to loop over each peak
    for(size_t i=0; i<2; ++i){
        Double_t peakValue = fitPeak(bkHist, bkPeaks.at(i), 0).at(1);
        if (peakValue == 0){
            continue;
        }
        else{
            retFit.push_back(bkPeaks.at(i));
            // Back-convert to channel number
            Double_t chanValue = (peakValue - coInpt)/coSlop;
            retChn.push_back(chanValue);
        }
    }

    // Clean up our mess
    bkHist->Delete();
    bkCanv->Close();

    return 0;
}

//// Create function that iterates over all clovers and crystals
//// Fits Eu peaks and then back-converts to channel number
//// Results written to .txt file
int fitEuPeaks(TTree* physTree, TTree* bkTree = 0)
{
    // Make my goddamn life easier because I can't reliably type 16 characters correctly
    std::string fitFile = "coFit/GeCoCalParam.txt";
    // Create vectors to hold peak and channel values 
    std::vector<Double_t> peakVec;
    std::vector<Double_t> chanVec;

    // Now we want to loop over clovers and crystals and fit the chosen Eu peaks
    // Loop over clovers
    for (int cloverNbr=1; cloverNbr<=4; cloverNbr++){
        // Loop over crystals
        for (int crystalNbr=1; crystalNbr<=4; crystalNbr++){
            // Call fit function
            fitEuPeak(physTree, fitFile, cloverNbr, crystalNbr, peakVec, chanVec);
            if (bkTree != 0){
                // Handle TFile argu to load correct files (?)
                // Call fit function
                //      Should pass back vectors with appended peak and channel number
                fitBkPeak(bkTree, fitFile, cloverNbr, crystalNbr, peakVec, chanVec);
            }
            // Write out energies and channels to text file
            if (bkTree == 0){
                std::string tempTitl = "euFit/euFit_DET"+ std::to_string(cloverNbr)+"_CRY"+ std::to_string(crystalNbr)+ ".txt";
                std::ofstream outputFile;
                outputFile.open(tempTitl);
                for (int totalToken=0; totalToken<static_cast<int>(chanVec.size()); totalToken++){
                    outputFile << chanVec.at(totalToken) << " " << peakVec.at(totalToken) << std::endl;
                }
                outputFile.close();
            }
            else{
                std::string tempTitl = "bkFit/bkFit_DET"+ std::to_string(cloverNbr)+"_CRY"+ std::to_string(crystalNbr)+ ".txt";
                std::ofstream outputFile;
                outputFile.open(tempTitl);
                for (int totalToken=0; totalToken<static_cast<int>(chanVec.size()); totalToken++){
                    outputFile << chanVec.at(totalToken) << " " << peakVec.at(totalToken) << std::endl;
                }
                outputFile.close();
            }
        }
    }
    return 0;
}

//// Just a script to generate the graph for each fit on each crystal
int fitEuFit(int cloverNbr, int crystalNbr, std::vector<Double_t>& linParVec, std::vector<Double_t>& qudParVec, const char* bkFind = 0)
{
    // This might work with TString, but I'm doing the constant char conversion just for the hell of it
    //      NOTE:  Talked with Greg, TString is better in the future
    std::string tempPath;
    TString tempTitl;
    TString tempTRes;
    TString tempSave;
    TString tempSubt = "Red = linear fit // Blue = quadratic fit";
    if (bkFind == 0){
        tempPath = "euFit/euFit_DET"+ std::to_string(cloverNbr)+"_CRY"+ std::to_string(crystalNbr)+ ".txt";
        tempTitl = Form("152Eu Fit for Clover %d Crystal %d",cloverNbr,crystalNbr);
        tempTRes = Form("152Eu Residuals for Clover %d Crystal %d",cloverNbr,crystalNbr);
        tempSave = Form("euFit/euFit_DET%d_CRY%d.png",cloverNbr,crystalNbr);
    }
    else {
        tempPath = "bkFit/bkFit_DET"+ std::to_string(cloverNbr)+"_CRY"+ std::to_string(crystalNbr)+ ".txt";
        tempTitl = Form("152Eu and Background Fit for Clover %d Crystal %d",cloverNbr,crystalNbr);
        tempTRes = Form("152Eu and Background Residuals for Clover %d Crystal %d",cloverNbr,crystalNbr);
        tempSave = Form("bkFit/bkFit_DET%d_CRY%d.png",cloverNbr,crystalNbr);
    }
    char const *charPath = tempPath.c_str();
    // clear the passed vectors
    qudParVec.clear();
    linParVec.clear();
            
    // now let's make the actual graph
    //      Create a canvas
    TCanvas* fitCanv = new TCanvas("fitCanv","fitCanv",600,1000);
    fitCanv->Divide(1,2);
    fitCanv->cd(1);
    //      Create our chan vs energy graph
    TGraph* texFit = new TGraph(charPath);
    texFit->SetMarkerStyle(20);
    texFit->SetTitle(tempTitl);
    texFit->GetXaxis()->SetTitle("Channel");
    texFit->GetYaxis()->SetTitle("Energy (keV)");
    texFit->Draw();

    //      Create our residuals graphs
    fitCanv->cd(2);
    TGraph* rsLFit = new TGraph();
    TGraph* rsQFit = new TGraph();
    //      Fit 4 lyfe
    //          Start with linear fit
    texFit->Fit("pol1");
    Int_t numPoints = texFit->GetN();
    for(int i=0; i<numPoints;++i){
        // pull x and y values
        Double_t x = texFit->GetX()[i]; 
        Double_t y = texFit->GetY()[i]; 
        // evaluate for fit functions
        Double_t lFit = texFit->GetFunction("pol1")->Eval(x); 
        // create point in TGraph
        rsLFit->SetPoint(i,x,lFit-y);
    }
    //          get linear fit parameters
    linParVec.push_back(texFit->GetFunction("pol1")->GetParameter(0));
    linParVec.push_back(texFit->GetFunction("pol1")->GetParameter(1));

    //          Now we do our quadratic fit
    texFit->Fit("pol2");
    for(int i=0; i<numPoints;++i){
        // pull x and y values
        Double_t x = texFit->GetX()[i]; 
        Double_t y = texFit->GetY()[i]; 
        // evaluate for fit functions
        Double_t qFit = texFit->GetFunction("pol2")->Eval(x);
        // create point in TGraph
        rsQFit->SetPoint(i,x,qFit-y);
    }
    //          get quad fit parameters
    qudParVec.push_back(texFit->GetFunction("pol2")->GetParameter(0));
    qudParVec.push_back(texFit->GetFunction("pol2")->GetParameter(1));
    qudParVec.push_back(texFit->GetFunction("pol2")->GetParameter(2));

    //          Make them look purrrty.png
    rsLFit->SetTitle(tempTRes);
    rsLFit->SetMarkerStyle(20);
    rsQFit->SetMarkerStyle(21);
    rsLFit->SetMarkerColor(2);
    rsQFit->SetMarkerColor(4);
    rsLFit->Draw("PA");
    rsQFit->Draw("P");

    // Save canvas
    fitCanv->SaveAs(tempSave);

    // Clean up
    fitCanv->Close();
    rsLFit->Delete();
    rsQFit->Delete();
    texFit->Delete();

    return 0;
}

////    YOU NEED TO CREATE CONTAINERS FOR FIT PARAMETERS AND WRITE TO A FILE
////    OTHERWISE PEOPLE HAVE TO READ OFF THE TERMINAL AND WRITE THINGS DOWN BY HAND
////    and that's not nice
////    ain't nobody got time for that

//// And now let's do it for everything
void fitEuFits()
{
    // vectors to pass
    std::vector<TString> toknVec;
    std::vector<Double_t> linParVec;
    std::vector<Double_t> qudParVec;
    std::vector<std::vector<Double_t>> totLinrVec;
    std::vector<std::vector<Double_t>> totQuadVec;
    // Loop over clovers
    for (int cloverNbr=1; cloverNbr<=4; cloverNbr++){
        // Loop over crystals
        for (int crystalNbr=1; crystalNbr<=4; crystalNbr++){
            // Call fit function
            fitEuFit(cloverNbr,crystalNbr,linParVec,qudParVec);
            totLinrVec.push_back(linParVec);
            totQuadVec.push_back(qudParVec);
            TString tempToken = Form("GETAMU_D%d_CRY%d_E", cloverNbr, crystalNbr);
            toknVec.push_back(tempToken);
        }
    }
    // Write calibration parameters to file
    //      Open linear file for writing
    std::ofstream outputFile;
    outputFile.open("euFit/linrCalParam.txt");

    //      Loop over vecs to create lin fit file
    for (int totalToken=0; totalToken<static_cast<int>(toknVec.size()); totalToken++){
        outputFile << toknVec.at(totalToken) << " " << totLinrVec.at(totalToken).at(0) << " " << totLinrVec.at(totalToken).at(1) << std::endl;
    }
    outputFile.close();
    //      Open quad file for writing
    outputFile.open("euFit/quadCalParam.txt");
    //      Loop over vecs to create qud fit file
    for (int totalToken=0; totalToken<static_cast<int>(toknVec.size()); totalToken++){
        outputFile << toknVec.at(totalToken) << " " << totQuadVec.at(totalToken).at(0) << " " << totQuadVec.at(totalToken).at(1) << " " << totQuadVec.at(totalToken).at(2) << std::endl;
    }
    outputFile.close();
    // Done?
}

//// Also for the background fits
void fitBkFits()
{
    // vectors to pass
    std::vector<TString> toknVec;
    std::vector<Double_t> linParVec;
    std::vector<Double_t> qudParVec;
    std::vector<std::vector<Double_t>> totLinrVec;
    std::vector<std::vector<Double_t>> totQuadVec;
    // Loop over clovers
    for (int cloverNbr=1; cloverNbr<=4; cloverNbr++){
        // Loop over crystals
        for (int crystalNbr=1; crystalNbr<=4; crystalNbr++){
            // Call fit function
            fitEuFit(cloverNbr,crystalNbr,linParVec,qudParVec,"back");
            totLinrVec.push_back(linParVec);
            totQuadVec.push_back(qudParVec);
            TString tempToken = Form("GETAMU_D%d_CRY%d_E", cloverNbr, crystalNbr);
            toknVec.push_back(tempToken);
        }
    }
    // Write calibration parameters to file
    //      Open linear file for writing
    std::ofstream outputFile;
    outputFile.open("bkFit/linrCalParam.txt");

    //      Loop over vecs to create lin fit file
    for (int totalToken=0; totalToken<static_cast<int>(toknVec.size()); totalToken++){
        outputFile << toknVec.at(totalToken) << " " << totLinrVec.at(totalToken).at(0) << " " << totLinrVec.at(totalToken).at(1) << std::endl;
    }
    outputFile.close();
    //      Open quad file for writing
    outputFile.open("bkFit/quadCalParam.txt");
    //      Loop over vecs to create qud fit file
    for (int totalToken=0; totalToken<static_cast<int>(toknVec.size()); totalToken++){
        outputFile << toknVec.at(totalToken) << " " << totQuadVec.at(totalToken).at(0) << " " << totQuadVec.at(totalToken).at(1) << " " << totQuadVec.at(totalToken).at(2) << std::endl;
    }
    outputFile.close();
    // Done.
}

//// I'm lazy and want a way to quickly draw calibrated spectra
int quikCalSpec(TTree* physTree, string fitFile, int cloverNbr, int crystalNbr)
{
    // Make sure random file is loaded
    if (gROOT->GetClass("randomBinDummy") == 0)
        gROOT->ProcessLine(".L randomBin.cxx+");

    // Read in parameters 
    std::vector<string> fitParam = getFitParam(fitFile, cloverNbr, crystalNbr);
    if (fitParam.empty())
        return 1;
    // Convert strings to doubles
    Double_t coSlop = std::stod(fitParam.at(2));
    Double_t coInpt = std::stod(fitParam.at(1));
    
    // Now we can create a histogram with the calibrated data
    TCanvas* quikCanv = new TCanvas("quikCanv","quikCanv",600,600);
    //      Make title and gating arguments
    TString quikArgu = Form("(fGeTAMU_Core_Energy+randomBin())*%f+%f",coSlop,coInpt);
    TString quikGate = Form("fGeTAMU_Core_CloverNbr_E==%d && fGeTAMU_Core_CrystalNbr_E==%d",cloverNbr,crystalNbr);
    //      Init and fill hist
    TH1D* quikHist = new TH1D("quikHist","",4096,0,4096);
    physTree->Project("quikHist",quikArgu,quikGate,"");
    quikHist->Draw();

    return 0;
}

//// meta function - ties everything together
void geCal(const char* coName, const char* euName, const char* bkName = 0)
{
    // open sesame
    TFile* coFile = TFile::Open(coName);
    TTree* coTree = static_cast<TTree*>(coFile->Get("T40Tree"));
    gSystem->Exec("mkdir -p ${PWD}/coFit/");
    TFile* euFile = TFile::Open(euName);
    TTree* euTree = static_cast<TTree*>(euFile->Get("T40Tree"));
    gSystem->Exec("mkdir -p ${PWD}/euFit/");
    TFile* bkFile = 0;
    TTree* bkTree = 0;
    if (bkName != 0){
        bkFile = TFile::Open(bkName);
        bkTree = static_cast<TTree*>(bkFile->Get("T40Tree"));
        gSystem->Exec("mkdir -p ${PWD}/bkFit/");
    }

    // execute coFit
    npFit(coTree);

    // execute euFit (and/or bkFit)
    fitEuPeaks(euTree,bkTree);
    if (bkName != 0){
        fitBkFits();
    }
    else{
        fitEuFits();
    }

    // clean up clean up everybody do your share
    coFile->Close();
    coFile->Delete();
    euFile->Close();
    euFile->Delete();
    if (bkName != 0){
        bkFile->Close();
        bkFile->Delete();
    }
}
