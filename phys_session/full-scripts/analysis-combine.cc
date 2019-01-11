// This program reads the *.root files of jet and hadron spectra from different pTHatBins
// and combine them together into the final spectra
//
// This program is very similar to analysis-spectra, details can be found there.

//C++ header
#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"

#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "JetScapeBanner.h"
#include "fjcore.hh"
#include "Pythia8/Pythia.h"

#include <GTL/dfs.h>
//ROOT headers
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"

using namespace std;
using namespace Jetscape;

//using namespace Pythia8;
int main(int argc, char* argv[]){
    
    int nListJets =1;
    int StartTime = time(NULL);
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    
    // Create a file on which histogram(s) can be saved.
    char inFileName[1000];
    char outFileName[1000];
    sprintf(outFileName,"CombinedSpectra.root");
    TFile* outFile = new TFile( outFileName, "RECREATE");
    
    //List here all pTHatBins for which you have generated "test_out.dat" and "SigmaHard.dat"
    int pTHatMin[48] = {1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90,  100,  110, 120, 130, 140, 150, 160, 170, 180, 190,200,                       210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};
    int pTHatMax[48] = {2, 3, 4, 5, 10,15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110,  120, 130, 140, 150, 160, 170, 180, 190, 200,210,                       220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 1380};
    
    ////By default we want to look at one pTHatBin spectrum
    //int pTHatMin[1]={110};
    //int pTHatMax[1]={120};
    int NpTHardBin = sizeof(pTHatMin)/sizeof(pTHatMin[0]);//sizeof(pTHatMin[]); or sizeof(pTHatMax[]);
    int Events[NpTHardBin];
    double HardCrossSection[NpTHardBin];
    double HardCrossSectionError[NpTHardBin];
    double DetectorEtaCut= 2.6;
    
    //Variables for jet spectrum
    //(For CMS at 2760 GeV with jet radius R=2,0.3, or 0.4, |eta_jet|<2.0)
    double JetpTBin[13]  ={70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300}; //in GeV
    int NpTJetBin = 13-1;     double JetpTMin = 10; //in GeV
    
    double JetEtaCut =2.0;
    double JetRadius = 0.3;
    double dNdpTCountJet[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]
    double pTHardBinJetBinError[NpTHardBin][NpTJetBin];
    double TotalDifferentialJetCrossSection[NpTJetBin];
    double TotalDifferentialJetCrossSectionError[NpTJetBin];
    
    //Variable for single hadron spectrum
    double SingleHadronpTBin[11] = {6.0, 7.2, 10.8, 14.4, 21.6, 28.8, 38.4, 48.0, 67.2, 86.4, 112.2};
    int NpTSingleHadronBin = 10;
    
    double SingleHadronEtaCut =1.0;
    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYield[NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin];
    
    double JetpTCut[NpTHardBin];
    int TriggeredJetNumber[NpTHardBin];
    //Variables for jet shape
    double JetShaperBin[7] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
    int NrJetShapeBin = 6;
    TH1D *CountJetShape[NpTHardBin];
    //Variables for  jet fragmentation function
    double JetFFzBin[11] = {0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1};
    int NzJetFFBin = 10;
    TH1D *CountJetFF[NpTHardBin];
    //Variables for jet mass
    double JetMassBin[13] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
    int NJetMassBin = 12;
    TH1D *CountJetMass[NpTHardBin];
    
    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    
    std::vector <fjcore::PseudoJet> fjInputs;
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius);
    
    Pythia8::Pythia pythia;//("",false);
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files
    // Instead of get histogram from test_out.dat directly, read from different root files for different pTHat bins
    //
    for (int k = 0; k<NpTHardBin; ++k)
    {
        //sprintf(inFileName,"SpectraBin.root");
        //sprintf(inFileName,"SpectraBin%i_%i.root", pTHatMin[k], pTHatMax[k]);
        sprintf(inFileName,"/Users/sscao/Work/JETSCAPE-public/combine/Root_Files_SpectraPP/SpectraBin%i_%i.root", pTHatMin[k], pTHatMax[k]);
        
        TFile *file = new TFile(inFileName);
        
        // Reset for a new pTHatBin
        char HistName[100];
        
        HistTempJet->Reset();
        sprintf(HistName,"CountVspTJetSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJet->SetName(HistName);
        
        HistTempSingleHadron->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadron->SetName(HistName);
        
        char NameOfHistogram[400]; // for reading in files
        sprintf(NameOfHistogram,"CountVspTJetSpectrumBin%d_%d",pTHatMin[k],pTHatMax[k]);
        HistTempJet = (TH1D *)file->FindObjectAny(NameOfHistogram);
        sprintf(NameOfHistogram,"CountVspTSingleHadronSpectrumBin%d_%d",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadron = (TH1D *)file->FindObjectAny(NameOfHistogram);
        
        char NameOfVector[400]; // for reading in files
        sprintf(NameOfVector,"EventInfo");
        TVector * EventInfo = (TVector*)file->FindObjectAny(NameOfVector);
        HardCrossSection[k] = (*EventInfo)[0];
        HardCrossSectionError[k]= (*EventInfo)[1];
        Events[k]= (*EventInfo)[2];
        cout << "Hard cross section: " << scientific << HardCrossSection[k] << " " << HardCrossSectionError[k]<< endl;
        
        sprintf(NameOfVector,"TriggeredJetInfo");
        TVector * TriggeredJetInfo = (TVector*)file->FindObjectAny(NameOfVector);
        JetpTCut[k] = (*TriggeredJetInfo)[0];
        TriggeredJetNumber[k] = (*TriggeredJetInfo)[1];
        cout << "Number of Jets with pT>" << int((*TriggeredJetInfo)[0])
        << " GeV/c: " << TriggeredJetNumber[k] << endl;
        
        //For jet shape
        sprintf(NameOfHistogram,"CountVsrJetShapeBin%d_%d",pTHatMin[k],pTHatMax[k]);
        CountJetShape[k] = (TH1D *)file->FindObjectAny(NameOfHistogram);
        CountJetShape[k]->SetDirectory(0);
        CountJetShape[k]->Sumw2();
        CountJetShape[k]->Scale( (1.0/Events[k]) , "width" );
        
        //For jet fragmentation function
        sprintf(NameOfHistogram,"CountVszJetFFBin%i_%i",pTHatMin[k],pTHatMax[k]);
        CountJetFF[k] = (TH1D *)file->FindObjectAny(NameOfHistogram);
        CountJetFF[k]->SetDirectory(0);
        CountJetFF[k]->Sumw2();
        CountJetFF[k]->Scale( (1.0/Events[k]) , "width" );

        //For jet mass
        sprintf(NameOfHistogram,"CountVsJetMassBin%i_%i",pTHatMin[k],pTHatMax[k]);
        CountJetMass[k] = (TH1D *)file->FindObjectAny(NameOfHistogram);
        CountJetMass[k]->SetDirectory(0);
        CountJetMass[k]->Sumw2();
        CountJetMass[k]->Scale( (1.0/Events[k]) , "width" );
        
        //For jet spectrum
        for(int j=0;j<NpTJetBin;j++)
        {
            dNdpTCountJet[k][j]= HistTempJet->GetBinContent(j+1);
            if(dNdpTCountJet[k][j]> 0.0)
            {
                pTHardBinJetBinError[k][j] = (dNdpTCountJet[k][j]*HardCrossSection[k]/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut))*TMath::Sqrt( (1/dNdpTCountJet[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
                cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountJet[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut)<<endl;
            }
            else
            {
                pTHardBinJetBinError[k][j] = 0.0;
                cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
            }
        }
        
        //For single Hadron spectrum
        for(int j=0;j<NpTSingleHadronBin;j++)
        {
            dNdpTCountSingleHadron[k][j]= HistTempSingleHadron->GetBinContent(j+1);
            if(dNdpTCountSingleHadron[k][j]> 0.0)
            {
                pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k]/(Events[k]*68*2*M_PI*((SingleHadronpTBin[j+1] + SingleHadronpTBin[j])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError[k]/HardCrossSection[k],2.0));
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*68*2*M_PI*((SingleHadronpTBin[j+1]+SingleHadronpTBin[j])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
            }
            else
            {
                pTHardBinSingleHadronBinError[k][j] = 0.0;
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
            }
        }
        

        //Write histogram into a root file
        
        // change pointer to outFile
        outFile->cd();
        
        HistTempJet->Write();
        HistTempSingleHadron->Write();
        CountJetShape[k]->Write();
        CountJetFF[k]->Write();
        //file->Close();
        
    }// End file for k-loop
    
    // now add contribution from different pTHat bins together
    //Plot of pTHatcrosssection vs pTHatBin
    double pTHatError[NpTHardBin], pTHatAvg[NpTHardBin];
    for(int k=0; k<NpTHardBin; k++)
    {
        pTHatError[k]=0.0;
        pTHatAvg[k]= ((pTHatMin[k] + pTHatMax[k])/2.0);
    }
    TGraphErrors *GHardCrossSection = new TGraphErrors(NpTHardBin,pTHatAvg,HardCrossSection,pTHatError,HardCrossSectionError);
    GHardCrossSection->SetName("HardCrossSection");
    GHardCrossSection->Write();
    
    //Plots for jet spectrum
    double DifferentialJetCrossSection[NpTJetBin],DifferentialJetCrossSectionError[NpTJetBin],JetpT[NpTJetBin],JetpTError[NpTJetBin];
    TGraphErrors * GEJet;
    for(int j=0; j<NpTJetBin; j++)
    {
        TotalDifferentialJetCrossSection[j] = 0.0;
        TotalDifferentialJetCrossSectionError[j] = 0.0;
    }
    
    for(int k=0;k<NpTHardBin;k++)
    {cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
        for(int j=0; j<NpTJetBin;j++)
        {
            DifferentialJetCrossSection[j] = (dNdpTCountJet[k][j]*HardCrossSection[k])/(Events[k]*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut);
            DifferentialJetCrossSectionError[j] = pTHardBinJetBinError[k][j];
            TotalDifferentialJetCrossSection[j] = TotalDifferentialJetCrossSection[j] + DifferentialJetCrossSection[j];
            TotalDifferentialJetCrossSectionError[j] = TotalDifferentialJetCrossSectionError[j] + TMath::Power( pTHardBinJetBinError[k][j], 2.0);
            JetpT[j] = (JetpTBin[j]+JetpTBin[j+1])/2.0;
            JetpTError[j] = (JetpTBin[j+1]-JetpTBin[j])/2.0;
            cout<<JetpT[j]<<"\t"<<DifferentialJetCrossSection[j]<<"\t"<<DifferentialJetCrossSectionError[j]<<endl;
        }
        
        GEJet = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetCrossSection,JetpTError,DifferentialJetCrossSectionError);
        char MyGraphName[100];
        sprintf(MyGraphName,"DifferentialJetCrossSectionBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJet->SetName(MyGraphName);
        GEJet->Write();
    }
    
    //Final Plot TotalJetCrossSection with errorbars
    cout<<"Final results for total crossSection"<<endl;
    cout<<"JetpT \t"<<"Total Differential CrossSection \t"<<"Error"<<endl;
    for(int j=0;j<NpTJetBin;j++)
    {
        TotalDifferentialJetCrossSectionError[j] =TMath::Sqrt(TotalDifferentialJetCrossSectionError[j]);
        cout<<JetpT[j]<<"\t"<<TotalDifferentialJetCrossSection[j]<<" +/- \t"<<TotalDifferentialJetCrossSectionError[j]<<endl;
    }
    cout<<"Create Final TGraphError plot "<<endl;
    //Final Plot and save also
    TGraphErrors *GEJetTotal = new TGraphErrors(NpTJetBin,JetpT,TotalDifferentialJetCrossSection,JetpTError,TotalDifferentialJetCrossSectionError);
    GEJetTotal->SetName("TotalDifferentialJetCrossSection");
    GEJetTotal->Write();
    
    // For charged Hadron spectrum
    double DifferentialSingleHadronYield[NpTSingleHadronBin],DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
    TGraphErrors * GESingleHadron;
    for(int j=0; j<NpTSingleHadronBin;j++)
    {
        TotalDifferentialSingleHadronYield[j] = 0.0;
        TotalDifferentialSingleHadronYieldError[j] = 0.0;
    }
    
    for(int k=0;k<NpTHardBin;k++)
    {cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
        for(int j=0; j<NpTSingleHadronBin;j++)
        {
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection[k])/(Events[k]*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*68*2*M_PI*((SingleHadronpTBin[j+1]+SingleHadronpTBin[j])/2.0)*2.0*SingleHadronEtaCut);
            DifferentialSingleHadronYieldError[j] = pTHardBinSingleHadronBinError[k][j];
            TotalDifferentialSingleHadronYield[j] = TotalDifferentialSingleHadronYield[j] + DifferentialSingleHadronYield[j];
            TotalDifferentialSingleHadronYieldError[j] = TotalDifferentialSingleHadronYieldError[j] + TMath::Power( pTHardBinSingleHadronBinError[k][j], 2.0);
            SingleHadronpT[j] = (SingleHadronpTBin[j+1] + SingleHadronpTBin[j])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1] - SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYield[j]<<"\t"<<DifferentialSingleHadronYieldError[j]<<endl;
        }
        
        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GESingleHadron->SetName(MyGraphName2);
        GESingleHadron->Write();
    }
    //Final Plot TotalDifferentialSingleHadronYield with errorbars
    cout<<"Final results for total differential single hadron yield"<<endl;
    cout<<"SingleHadronpT \t"<<"total differential single hadron yield \t"<<"Error"<<endl;
    for(int j=0;j<NpTSingleHadronBin;j++)
    {
        TotalDifferentialSingleHadronYieldError[j] =TMath::Sqrt(TotalDifferentialSingleHadronYieldError[j]);
        cout<<SingleHadronpT[j]<<"\t"<<TotalDifferentialSingleHadronYield[j]<<" +/- \t"<<TotalDifferentialSingleHadronYieldError[j]<<endl;
    }
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    TGraphErrors *GESingleHadronTotal
    = new TGraphErrors(NpTSingleHadronBin, SingleHadronpT, TotalDifferentialSingleHadronYield, SingleHadronpTError, TotalDifferentialSingleHadronYieldError);
    GESingleHadronTotal->SetName("TotalDifferentialSingleHadronYield");
    GESingleHadronTotal->Write();
    
    // For Jet Shape
    double n_jetshape = 0;
    TH1D *HistJetShape = new TH1D("TotlJetShape", "Totl Jet Shape", NrJetShapeBin, JetShaperBin);
    HistJetShape->Sumw2();
    for(int k=0;k<NpTHardBin;k++){
        n_jetshape += TriggeredJetNumber[k]*HardCrossSection[k]/Events[k];
        // *Usually Errors of Hard Cross Sections are Much Smaller than Statistical Errors*
        HistJetShape->Add( CountJetShape[k], HardCrossSection[k] ); // accumulate distribution weighted by cross section
    }
    HistJetShape->Scale( 1/n_jetshape );
    double ErrorNormJetShape;
    double NormJetShape
    = HistJetShape->IntegralAndError( 1, NrJetShapeBin,
                                      ErrorNormJetShape,
                                      "width" );
    TH1D *Norm = (TH1D*)HistJetShape->Clone("Normalization");
    Norm->Sumw2();Norm->SetTitle("Normalization");
    int nbins = Norm->GetSize();
    for( int i=1; i < nbins-1; i++){
        Norm->SetBinContent(i, NormJetShape);
        Norm->SetBinError(i, ErrorNormJetShape);
    }
    HistJetShape->Divide(Norm);
    delete Norm;

    double JetShape[NrJetShapeBin], JetShapeError[NrJetShapeBin], JetShapeR[NrJetShapeBin], JetShapeRError[NrJetShapeBin];
    cout<<"Final results for  Jet Shape is Below ( pTjet >"<< int(JetpTCut[0]) << " GeV/c)" <<endl;
    cout<<"r \t rho(r) \t"<<"Error"<<endl;
    for(int j=0;j<NrJetShapeBin;j++)
    {
        JetShape[j] = HistJetShape->GetBinContent(j+1);
        JetShapeError[j] = HistJetShape->GetBinError(j+1);
        JetShapeR[j] = (JetShaperBin[j]+JetShaperBin[j+1])/2.0;
        JetShapeRError[j] = (JetShaperBin[j+1]-JetShaperBin[j])/2.0;
        cout<<JetShapeR[j]<<"\t"<<JetShape[j]<<"\t"<<JetShapeError[j]<<endl;
    }
    delete HistJetShape;
    for(int k=0;k<NpTHardBin;k++){delete CountJetShape[k];}
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    TGraphErrors * GEJetShape
    = new TGraphErrors( NrJetShapeBin, JetShapeR, JetShape, JetShapeRError, JetShapeError );
    GEJetShape->SetName("TotalJetShape");
    GEJetShape->Write();
    
    // For Jet Fragmentation Function
    double n_ff = 0;
    TH1D *HistJetFF = new TH1D("TotlJetFF", "Totl Jet FF", NzJetFFBin, JetFFzBin);
    HistJetFF->Sumw2();
    for(int k=0;k<NpTHardBin;k++){
        n_ff += TriggeredJetNumber[k]*HardCrossSection[k]/Events[k];
        // *Usually Errors of Hard Cross Sections are Much Smaller than Statistical Errors*
        HistJetFF->Add( CountJetFF[k], HardCrossSection[k] );
    }
    HistJetFF->Scale( 1/n_ff );
    double JetFF[NzJetFFBin], JetFFError[NzJetFFBin], JetZ[NzJetFFBin], JetZError[NzJetFFBin];
    TGraphErrors * GEJetFF;
    cout<<"Final results for  Jet Fragmentation Function is Below ( pTjet >"<< int(JetpTCut[0]) << " GeV/c)" <<endl;
    cout<<"z \t (1/N_jet)dN/dz \t"<<"Error"<<endl;
    for(int j=0;j<NzJetFFBin;j++)
    {
        JetFF[j] = HistJetFF->GetBinContent(j+1);
        JetFFError[j] = HistJetFF->GetBinError(j+1);
        JetZ[j] = (JetFFzBin[j]+JetFFzBin[j+1])/2.0;
        JetZError[j] = (JetFFzBin[j+1]-JetFFzBin[j])/2.0;
        cout<<JetZ[j]<<"\t"<<JetFF[j]<<"\t"<<JetFFError[j]<<endl;
    }
//    delete HistJetFF;
//    for(int k=0;k<NpTHardBin;k++){delete CountJetFF[k];}
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    GEJetFF
    = new TGraphErrors(NzJetFFBin,JetZ,JetFF,JetZError,JetFFError);
    GEJetFF->SetName("TotalJetFragmentationFunction");
    GEJetFF->Write();
    
    
    // For Jet Mass
    double n_jetmass = 0;
    TH1D *HistJetMass = new TH1D("TotlJetMass", "Totl Jet Mass", NJetMassBin, JetMassBin);
    HistJetMass->Sumw2();
    for(int k=0;k<NpTHardBin;k++){
        n_jetmass += TriggeredJetNumber[k]*HardCrossSection[k]/Events[k];
        // *Usually Errors of Hard Cross Sections are Much Smaller than Statistical Errors*
        HistJetMass->Add( CountJetMass[k], HardCrossSection[k] );
    }
    HistJetMass->Scale( 1/n_jetmass );

    double JetMassDist[NJetMassBin], JetMassDistError[NJetMassBin], JetMass[NJetMassBin], JetMassError[NJetMassBin];
    TGraphErrors * GEJetMass;
    cout<<"Final results for  Jet Mass is Below ( pTjet >"<< int(JetpTCut[0]) << " GeV/c)" <<endl;
    cout<<"m \t (1/N_jet)dN_jet/dm \t"<<"Error"<<endl;
    for(int j=0; j<NJetMassBin;j++){
        JetMassDist[j] = HistJetMass->GetBinContent(j+1);
        JetMassDistError[j] = HistJetMass->GetBinError(j+1);
        JetMass[j] = (JetMassBin[j]+JetMassBin[j+1])/2.0;
        JetMassError[j] = (JetMassBin[j+1]-JetMassBin[j])/2.0;
        cout<<JetMass[j]<<"\t"<<JetMassDist[j]<<"\t"<<JetMassDistError[j]<<endl;
    }
    delete HistJetMass;
    for(int k=0;k<NpTHardBin;k++){delete CountJetMass[k];}
    cout<<"Create Final TGraphError plot "<<endl;
    // Final Plot and save also
    GEJetMass
    = new TGraphErrors(NJetMassBin,JetMass,JetMassDist,JetMassError,JetMassDistError);
    GEJetMass->SetName("TotalJetMass");
    GEJetMass->Write();

    delete GEJet; delete GESingleHadron;
    delete GEJetTotal; delete GESingleHadronTotal;
    delete GEJetShape;
    delete GEJetFF;
    delete GEJetMass;
    delete GHardCrossSection;
    
    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    
    return 0;
}
