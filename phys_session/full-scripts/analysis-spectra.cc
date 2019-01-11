// This program reads "test_out.dat" produced from the JETSCAPE and Generates spectrums for jets and charged hadron yield.
// JETSCAPE by default gives pTHatCrossSection in milibarn and particle's momentum in GeV.
// With slight change, this program also allows reading multiple pTHatBins data and produce weighted spectra. 
// Output of this program is a ROOT file which contains following plots
// 1. Number of jets vs pT of the jet. (Graph name is CountVspTJetSpectrumBinpTHatMin_pTHatMax)
// 2. Number of charged hadron vs pT of the hadron. (Graph name is CountVspTSingleHadronSpectrumBinpTHatMin_pTHatMax)
// 3. Weighted differential jet cross section vs pT of the jet. Here Weighted differential crosssection means = sigmapTHat*dN/(dpTJet*dEta*TotalEvents)
//    (Graph name is DifferentialJetCrossSectionBinpTHatMin_pTHatMax)
//
// 4. Weighted charged hadron differential yield vs pT of hadron. Here Weighted differential Yield = sigmapTHat*dN/(TotalInelasticCrossSection*2*PI*dpTHadron*dEta*TotalEvents)  
//    (Graph name is DifferentialSingleHadronYieldBinpTHatMin_pTHatMax)
//
// 5. Total differential jet cross section (i.e. summed over all pTHat bins) vs pT of the jet. (if you have multiple pTHatBins data)
//    (Graph name is TotalDifferentialJetCrossSection)
//
// 6. Total differential charged hadron yield (i.e. summed over all pTHat bins) vs pT of the hadron.
//    (Graph name is TotalDifferentialSingleHadronYield)
//
// 7. pTHat crosssection i.e. hard scattering crosssection vs pTHat bin.
//    (Graph name is HardCrossSection) 
// 
// Note, plots are saved in two ROOT format TH1D and TGraphErrors.
// For jet spectrum, we use anti-kT algorithm with jet radius, R=0.3 and |eta_jet|<2.0. Inside jet cone, we include all particle except neutrinos (CMS def).
// For charged hadron spectrum, we use |eta_hadron|<1.0. Only charged hadrons are included in the spectrum.

// Authorship: written by Amit Kumar, revised by Shanshan Cao.

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
    
    int nListJets =1;  // for output control
    int StartTime = time(NULL);
    // Create the ROOT application environment.
    TApplication theApp("hist", &argc, argv);
    
    ////List here all pTHatBins for which you have generated "test_out.dat"
    //int pTHatMin[48] = {1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90,  100,  110, 120, 130, 140, 150, 160, 170, 180, 190,200,                       210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600};
    //int pTHatMax[48] = {2, 3, 4, 5, 10,15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110,  120, 130, 140, 150, 160, 170, 180, 190, 200,210,                       220, 230, 240, 250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 1380};
    
    //By default we want to look at only one pTHatBin spectrum
    int pTHatMin[1]={110};
    int pTHatMax[1]={120};
    int NpTHardBin = sizeof(pTHatMin)/sizeof(pTHatMin[0]); //sizeof(pTHatMin[]); or sizeof(pTHatMax[]); // # of bin we analyze
    double DetectorEtaCut= 2.6;
    
    //Variables for jet spectrum
    //(For CMS at 2760 GeV with jet radius R=2,0.3, or 0.4, |eta_jet|<2.0)
    double JetpTBin[13] = {70, 80, 90, 100, 110, 130, 150, 170, 190, 210, 240, 270, 300}; //in GeV 
    int NpTJetBin = 13-1;     double JetpTMin = 10; //in GeV
    
    double JetEtaCut = 2.0;
    double JetRadius = 0.3;
    double dNdpTCountJet[NpTHardBin][NpTJetBin];  //[ptHatBin] [Regular pt]  // count for jet in each bin
    double pTHardBinJetBinError[NpTHardBin][NpTJetBin];        // error
    double TotalDifferentialJetCrossSection[NpTJetBin];        // weighted cross section 
    double TotalDifferentialJetCrossSectionError[NpTJetBin];   // weighted error
    
    //Variable for single hadron spectrum
    double SingleHadronpTBin[11] = {6.0, 7.2, 10.8, 14.4, 21.6, 28.8, 38.4, 48.0, 67.2, 86.4, 112.2};
    int NpTSingleHadronBin = 10;
    
    double SingleHadronEtaCut =1.0;
    double dNdpTCountSingleHadron[NpTHardBin][NpTSingleHadronBin];  //[ptHatBin] [Regular pt]
    double pTHardBinSingleHadronBinError[NpTHardBin][NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYield[NpTSingleHadronBin];
    double TotalDifferentialSingleHadronYieldError[NpTSingleHadronBin];
   
    // for jet substructure
    double JetpTCut = 100.0;
    //Variables for jet shape
    double JetShaperBin[7] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
    int NrJetShapeBin = 6;
    //Variables for  jet fragmentation function
    double JetFFzBin[11] = {0.01, 0.016, 0.025, 0.04, 0.063, 0.1, 0.16, 0.25, 0.4, 0.63, 1};
    int NzJetFFBin = 10;
    //Variables for jet mass
    double JetMassBin[13] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};
    int NJetMassBin = 12;
    
    // Histograms.1. Number of jets vs pT of the jet. 2. Number of charged hadron vs pT of the hadron
    TH1D *HistTempJet = new TH1D("JetSpectrumBin", "Jet Spectrum pT", NpTJetBin, JetpTBin); //CountVspT for jets
    TH1D *HistTempSingleHadron = new TH1D("SingleHadronSpectrumBin", "Single Hadron Spectrum pT", NpTSingleHadronBin, SingleHadronpTBin); //CountVspT for single-hadron
    TH1D *HistTempJetShape = new TH1D("JetShapeBin", "Jet Shape", NrJetShapeBin, JetShaperBin); //CountVsr for Jet Shape
    TH1D *HistTempJetFF = new TH1D("JetFragmentaionFunctionBin", "Jet Fragmentaion Function", NzJetFFBin, JetFFzBin);
    TH1D *HistTempJetMass = new TH1D("JetMassBin", "Jet Mass", NJetMassBin, JetMassBin);
    
    std::vector <fjcore::PseudoJet> fjInputs;
    std::vector <int> chargeList;
    fjcore::JetDefinition jetDef(fjcore::antikt_algorithm, JetRadius); // define jet finding algorithm
    
    Pythia8::Pythia pythia;//("",false);
    
    cout<<"These are pTHat loops "<<endl;
    // For loop to open different pTHat bin files, in this default example, only one file
    for (int k = 0; k<NpTHardBin; ++k)
    {
        char HadronFile[300], pTBinString[100];
        
        //sprintf(HadronFile,"test_outBin%i_%i.dat", pTHatMin[k], pTHatMax[k]);
        sprintf(HadronFile,"test_out.dat");  // name of the input file
        
        auto myfile  = make_shared<JetScapeReaderAscii>(HadronFile);
        sprintf(pTBinString,"Current pTHatBin is %i (%i,%i) GeV",k,pTHatMin[k],pTHatMax[k]);
        
        int  SN=0,PID=0,Status=0;
        double Px, Py, Pz, E, Eta, Phi;
        int Events =0;
        int TriggeredJetNumber=0;
        
        // Create a file on which histogram(s) can be saved.
        char outFileName[1000];
        sprintf(outFileName,"SpectraBin%i_%i.root",pTHatMin[k],pTHatMax[k]); // name of the output root file
        TFile* outFile = new TFile( outFileName, "RECREATE");
        // Reset for each pTHardBin
        char HistName[100];
        
        HistTempJet->Reset();
        sprintf(HistName,"CountVspTJetSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJet->SetName(HistName);
        
        HistTempSingleHadron->Reset();
        sprintf(HistName,"CountVspTSingleHadronSpectrumBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempSingleHadron->SetName(HistName);
        
        HistTempJetShape->Reset();
        sprintf(HistName,"CountVsrJetShapeBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetShape->SetName(HistName);
        HistTempJetShape->Sumw2();
        
        HistTempJetFF->Reset();
        sprintf(HistName,"CountVszJetFFBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetFF->SetName(HistName);
        HistTempJetFF->Sumw2();
        
        HistTempJetMass->Reset();
        sprintf(HistName,"CountVsJetMassBin%i_%i",pTHatMin[k],pTHatMax[k]);
        HistTempJetMass->SetName(HistName);
        HistTempJetMass->Sumw2();
        
        fjInputs.resize(0);
        chargeList.resize(0);
        
        // Read file
        vector<shared_ptr<Hadron>> hadrons;
        
        while (!myfile->Finished())
        {
            myfile->Next(); // default JETSCAPE function that takes one event from test_out.dat
            hadrons = myfile->GetHadrons(); // default JETSCAPE function that takes hadron information from the event
            //cout<<"Number of hadrons is: " << hadrons.size() << endl;
            Events++;
            for(unsigned int i=0; i<hadrons.size(); i++)
            {
                SN = i;
                PID= hadrons[i].get()->pid();
                E  = hadrons[i].get()->e();
                Px = hadrons[i].get()->px();
                Py = hadrons[i].get()->py();
                Pz = hadrons[i].get()->pz();
                Eta = hadrons[i].get()->eta();
                Phi = hadrons[i].get()->phi();
		Status = hadrons[i].get()->pstat();
                double PT = TMath::Sqrt( (Px*Px) + (Py*Py));
                
                // Add this particle for Jet spectrum
                if( SN==0 && fjInputs.size()>0 )  // construct jet for the previous event
                {
                    // Run Fastjet algorithm and sort jets in pT order.
                    vector <fjcore::PseudoJet> UnSortedJets, SortedJets;
                    fjcore::ClusterSequence clustSeq(fjInputs, jetDef);
                    UnSortedJets = clustSeq.inclusive_jets(JetpTMin);
                    SortedJets    = sorted_by_pt(UnSortedJets);
                    
                    // List first few FastJet jets and some info about them.
                    if (nListJets)
                    {
                        cout << "\n --------  FastJet jets, R = " << JetRadius << "anti-kt for pTHatBin = "<<k
                        << "  --------------------------------------------------\n\n "
                        << "  i         pT        y     eta      phi  " << endl;
                        for (int i = 0; i < int(SortedJets.size()); ++i)
                        {
                            vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();
                            cout << setw(4) << i << fixed << setprecision(3) << setw(11)
                            << SortedJets[i].perp() << setw(9)  << SortedJets[i].rap()
                            << setw(9) << SortedJets[i].eta() << setw(9)  << SortedJets[i].phi_std() << endl;
                        }
                        cout << "\n --------  End FastJet Listing  ------------------"
                        << "---------------------------------" << endl;
                    }
                    
                    int pFast = SortedJets.size();
                    for (int i = 0; i < pFast; ++i)  // fill jets from this event into observable bins
                    {
                        if(-JetEtaCut < SortedJets[i].eta() && SortedJets[i].eta()< JetEtaCut )
                        {
                            HistTempJet->Fill( SortedJets[i].perp() ); // fill the pT spectrum
                            if( SortedJets[i].perp() >= JetpTCut ){ // now for substructure

                                TriggeredJetNumber ++;
                                vector<fjcore::PseudoJet> constituents = SortedJets[i].constituents();

                                //Jet Shape---
                                for( int j = 0; j < constituents.size(); j++ ){
                                    
                                    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    HistTempJetShape->Fill( delta_r, constituents[j].perp() ); // fill the pT vs. r, second variable is the weight to fill in
                                }
                                
                                //Jet Fragmentation Function---
                                for( int j = 0; j < fjInputs.size(); j++ ){
                                    double delta_eta = fjInputs[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(fjInputs[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    if( fabs(chargeList[j]) > 0.01 && delta_r <= JetRadius){//charged particle in jet cone
                                        double z_jet = fjInputs[j].perp()/SortedJets[i].perp();
                                        HistTempJetFF->Fill( z_jet ); // fill in z distribution
                                    }
                                }

                                //Jet Mass---
                                double jet_e = 0.0, jet_px = 0.0, jet_py = 0.0, jet_pz = 0.0;
                                for( int j = 0; j < constituents.size(); j++ ){
                                    double delta_eta = constituents[j].eta() - SortedJets[i].eta();
                                    double delta_phi = SortedJets[i].delta_phi_to(constituents[j]);
                                    double delta_r = TMath::Sqrt( delta_eta*delta_eta + delta_phi*delta_phi);
                                    if( delta_r <= JetRadius){// "ALL" particle in jet cone
                                        jet_e += constituents[j].e();
                                        jet_px += constituents[j].px();
                                        jet_py += constituents[j].py();
                                        jet_pz += constituents[j].pz();
                                    }
                                }
                                double jet_mass
                                = TMath::Sqrt(jet_e*jet_e - jet_px*jet_px - jet_py*jet_py - jet_pz*jet_pz);
                                HistTempJetMass->Fill( jet_mass ); // fill in mass distribution

                            }
                        }
                    }
                    fjInputs.resize(0);  // reset jet information and put the first particle into a new jet list
                    //cout<<"Found a Jet \t "<<pTBinString<<"\t NetJetevents is \t"<<NetJetEvents<<endl;
                    if(Status>=0 && fabs(Eta) < DetectorEtaCut &&  PT>0.01  && PID!=12 && PID!=14 && PID!=16 && PID!=18){
                        fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                        chargeList.push_back( pythia.particleData.charge( PID ) );
                    }
		    // note for Status>=0: here we only analyze positive (jet + recoil) partons in the particle list in Pb-Pb collisions, in principle, one should subtract negative (back-reaction) particle contribution in realistic calcualtions
                }
                else
                {
                    //cout<<EventLabel << " "<<PID <<" 1"<< E <<" "<< Px<<" " << Py<<" " << Pz<<" " << Eta<<" " << Phi<<endl;
                    if(Status>=0 && fabs(Eta) < DetectorEtaCut && PT>0.01  &&  PID!=12 && PID!=14 && PID!=16 && PID!=18 ){
                        fjInputs.push_back(fjcore::PseudoJet(Px,Py,Pz,E));
                        chargeList.push_back( pythia.particleData.charge( PID ) );
                    }
                }
                
                
                // Add this particle into SingleHadron spectrum
                if(Status>=0 && fabs(Eta) < SingleHadronEtaCut && PT>0.01  && fabs(PID) >100 &&  pythia.particleData.charge( PID )!=0   )
                {
                    //cout<<"PID "<<PID<<"\t charge = "<<pythia.particleData.charge( PID)<<endl;
                    HistTempSingleHadron->Fill(PT); // fill pT distribution of single hadron
                }
                
            } // end of one event
        } // end of reading one file
        
        //Extract cross section for this pTHatBin from the last few lines of the test_out.dat file
        string command_line = string("tail -20 ") + HadronFile + string(" > tail_tempFile.dat");
        system(command_line.c_str());
        
        ifstream input("tail_tempFile.dat");
        
        string line;
        
        bool FoundVal = false;
        bool FoundError = false;
        
        double sigma_val, sigma_err;
        
        while(getline(input,line)) {
            
            char str[1024];
            strncpy(str, line.c_str(), sizeof(str));
            
            if( str[0] != '#' ) continue;
            
            char* divideChar[5];
            divideChar[0] = strtok(str," ");
            
            int nChar=0;
            while(divideChar[nChar]!=NULL) {
                if(nChar==4) break;
                nChar++;
                divideChar[nChar] = strtok(NULL," ");
            }
            
            if(nChar==4) {
                string firstStr(divideChar[0]);
                string secondStr(divideChar[1]);
                string thirdStr(divideChar[2]);
                string fourthStr(divideChar[3]);
                if(thirdStr=="=") {
                    if(secondStr=="sigmaGen") {
                        sigma_val=atof(divideChar[3]);
                        FoundVal=true;
                    } else if (secondStr=="sigmaErr") {
                        sigma_err=atof(divideChar[3]);
                        FoundError=true;
                    }
                }
            }
            
            if(FoundVal && FoundError) break;
            
        }
        
        if(FoundVal && FoundError) {
            cout << scientific;
            cout << "sigmaGen:  " << sigma_val << endl;
            cout << "sigmaErr:  " << sigma_err << endl;
        } else {
            cout << "cannot found sigma information" << endl;
            exit(EXIT_FAILURE);
        }
        
        input.close();
        
        system("rm tail_tempFile.dat");
        
        double HardCrossSection = sigma_val;
        double HardCrossSectionError = sigma_err;
        
        // end of reading cross section
        
        //For jet spectrum, weighted by cross section and combined through multiple pTHatBins
        for(int j=0;j<NpTJetBin;j++)
        {
            dNdpTCountJet[k][j]= HistTempJet->GetBinContent(j+1);
            if(dNdpTCountJet[k][j] > 0.0)
            {
                pTHardBinJetBinError[k][j] = (dNdpTCountJet[k][j]*HardCrossSection/(Events*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut))*TMath::Sqrt( (1/dNdpTCountJet[k][j]) + TMath::Power(HardCrossSectionError/HardCrossSection,2.0));
                cout<<"For JetBin j = "<<j<<" \t BinContent = \t"<<HistTempJet->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountJet[k][j]*HardCrossSection)/(Events*(JetpTBin[j+1]-JetpTBin[j])*2.0*JetEtaCut)<<endl;
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
                pTHardBinSingleHadronBinError[k][j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection/(Events*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut))*TMath::Sqrt( (1/dNdpTCountSingleHadron[k][j]) + TMath::Power(HardCrossSectionError/HardCrossSection,2.0));
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<(dNdpTCountSingleHadron[k][j]*HardCrossSection)/(Events*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*2.0*SingleHadronEtaCut)<<endl;
            }
            else
            {
                pTHardBinSingleHadronBinError[k][j] = 0.0;
                cout<<"For SingleHadronBin j = "<<j<<" \t BinContent = \t"<<HistTempSingleHadron->GetBinContent(j+1)<<"\t Scaled Value = "<<0.0<<endl;
            }
        }
        
        //Write histogram into a root file
        HistTempJet->Write();
        HistTempJetShape->Write();
        HistTempSingleHadron->Write();
        HistTempJetFF->Write();
        HistTempJetMass->Write();
        
        myfile->Close();
        
        TVector EventInfo(3);
        EventInfo[0] = HardCrossSection;
        EventInfo[1] = HardCrossSectionError;
        EventInfo[2] = Events;
        EventInfo.Write("EventInfo");
        
        TVector TriggeredJetInfo(2);
        TriggeredJetInfo[0] = JetpTCut;
        TriggeredJetInfo[1] = TriggeredJetNumber;
        TriggeredJetInfo.Write("TriggeredJetInfo");
        
        //Plots for jet spectrum
        double DifferentialJetCrossSection[NpTJetBin],DifferentialJetCrossSectionError[NpTJetBin],JetpT[NpTJetBin],JetpTError[NpTJetBin];
        TGraphErrors * GEJet;
        
        cout<<"For ptHardBin = "<<k+1<<"\t CrossSection is Below "<<endl;
        for(int j=0; j<NpTJetBin;j++)  // fill in x, y observables together with error
        {
            DifferentialJetCrossSection[j] = (dNdpTCountJet[k][j]*HardCrossSection)/(Events*((JetpTBin[j]+JetpTBin[j+1])/2.0)*2.0*JetEtaCut);
            DifferentialJetCrossSectionError[j] = pTHardBinJetBinError[k][j];
            JetpT[j] = (JetpTBin[j] + JetpTBin[j+1])/2.0;
            JetpTError[j] = (JetpTBin[j+1] - JetpTBin[j])/2.0;
            cout<<JetpT[j]<<"\t"<<DifferentialJetCrossSection[j]<<"\t"<<DifferentialJetCrossSectionError[j]<<endl;
        }
        
        // put observables and error bars into a root graph
	GEJet = new TGraphErrors(NpTJetBin,JetpT,DifferentialJetCrossSection,JetpTError,DifferentialJetCrossSectionError);
        char MyGraphName[100];
        sprintf(MyGraphName,"DifferentialJetCrossSectionBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJet->SetName(MyGraphName);
        GEJet->Write();
        
        // For charged Hadron spectrum
        double DifferentialSingleHadronYield[NpTSingleHadronBin],DifferentialSingleHadronYieldError[NpTSingleHadronBin],SingleHadronpT[NpTSingleHadronBin],SingleHadronpTError[NpTSingleHadronBin];
        TGraphErrors * GESingleHadron;
        
        cout<<"For ptHardBin = "<<k+1<<"\t SingleHadron differential yield is Below "<<endl;
        for(int j=0; j<NpTSingleHadronBin;j++)
        {
            DifferentialSingleHadronYield[j] = (dNdpTCountSingleHadron[k][j]*HardCrossSection)/(Events*(SingleHadronpTBin[j+1]-SingleHadronpTBin[j])*68*2*M_PI*((SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0)*2.0*SingleHadronEtaCut);
            DifferentialSingleHadronYieldError[j] = pTHardBinSingleHadronBinError[k][j];
            SingleHadronpT[j] = (SingleHadronpTBin[j]+SingleHadronpTBin[j+1])/2.0;
            SingleHadronpTError[j] = (SingleHadronpTBin[j+1]-SingleHadronpTBin[j])/2.0;
            cout<<SingleHadronpT[j]<<"\t"<<DifferentialSingleHadronYield[j]<<"\t"<<DifferentialSingleHadronYieldError[j]<<endl;
        }
        
        GESingleHadron = new TGraphErrors(NpTSingleHadronBin,SingleHadronpT,DifferentialSingleHadronYield,SingleHadronpTError,DifferentialSingleHadronYieldError);
        char MyGraphName2[100];
        sprintf(MyGraphName2,"DifferentialSingleHadronYieldBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GESingleHadron->SetName(MyGraphName2);
        GESingleHadron->Write();

        
        // for jet shape
        HistTempJetShape->Scale( (1.0/TriggeredJetNumber), "width" ); // average over number of triggered jets
        double ErrorNormJetShape;
        double NormJetShape
        = HistTempJetShape->IntegralAndError( 1, NrJetShapeBin,
                                              ErrorNormJetShape,
                                              "width" );
        TH1D *Norm = (TH1D*)HistTempJetShape->Clone("Normalization");
        Norm->Sumw2();Norm->SetTitle("Normalization");
        int nbins = Norm->GetSize();
        for( int i=1; i < nbins-1; i++){
            Norm->SetBinContent(i, NormJetShape);
            Norm->SetBinError(i, ErrorNormJetShape);
        }
        HistTempJetShape->Divide(Norm); // jet shape is defined as a self-normalized observable
        delete Norm;

        double JetShape[NrJetShapeBin],JetShapeError[NrJetShapeBin],JetShapeR[NrJetShapeBin],JetShapeRError[NrJetShapeBin];
        TGraphErrors * GEJetShape;
        
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Shape is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NrJetShapeBin;j++)
        {
            JetShape[j] = HistTempJetShape->GetBinContent(j+1);
            JetShapeError[j] = HistTempJetShape->GetBinError(j+1);
            JetShapeR[j] = (JetShaperBin[j]+JetShaperBin[j+1])/2.0;
            JetShapeRError[j] = (JetShaperBin[j+1]-JetShaperBin[j])/2.0;
            cout<<JetShapeR[j]<<"\t"<<JetShape[j]<<"\t"<<JetShapeError[j]<<endl;
        }
        GEJetShape
        = new TGraphErrors(NrJetShapeBin,JetShapeR,JetShape,JetShapeRError,JetShapeError);
        char MyGraphName3[100];
        sprintf(MyGraphName3,"JetShapeBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetShape->SetName(MyGraphName3);
        GEJetShape->Write();
        
        
        
        // for jet fragmentation function
        HistTempJetFF->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetFF[NzJetFFBin], JetFFError[NzJetFFBin], JetZ[NzJetFFBin], JetZError[NzJetFFBin];
        TGraphErrors * GEJetFF;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Fragmentation Function is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NzJetFFBin;j++){
            JetFF[j] = HistTempJetFF->GetBinContent(j+1);
            JetFFError[j] = HistTempJetFF->GetBinError(j+1);
            JetZ[j] = (JetFFzBin[j]+JetFFzBin[j+1])/2.0;
            JetZError[j] = (JetFFzBin[j+1]-JetFFzBin[j])/2.0;
            cout<<JetZ[j]<<"\t"<<JetFF[j]<<"\t"<<JetFFError[j]<<endl;
        }
        GEJetFF
        = new TGraphErrors(NzJetFFBin,JetZ,JetFF,JetZError,JetFFError);
        char MyGraphName4[100];
        sprintf(MyGraphName4,"JetFFBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetFF->SetName(MyGraphName4);
        GEJetFF->Write();
        
        // for jet mass
        HistTempJetMass->Scale( (1.0/TriggeredJetNumber), "width" );
        double JetMassDist[NJetMassBin], JetMassDistError[NJetMassBin], JetMass[NJetMassBin], JetMassError[NJetMassBin];
        TGraphErrors * GEJetMass;
        cout<<"For ptHardBin = "<<k+1<<"\t Jet Mass Distribution is Below ( pTjet >"<< int(JetpTCut) << " GeV/c)" <<endl;
        for(int j=0; j<NJetMassBin;j++){
            JetMassDist[j] = HistTempJetMass->GetBinContent(j+1);
            JetMassDistError[j] = HistTempJetMass->GetBinError(j+1);
            JetMass[j] = (JetMassBin[j]+JetMassBin[j+1])/2.0;
            JetMassError[j] = (JetMassBin[j+1]-JetMassBin[j])/2.0;
            cout<<JetMass[j]<<"\t"<<JetMassDist[j]<<"\t"<<JetMassDistError[j]<<endl;
        }
        GEJetMass
        = new TGraphErrors(NJetMassBin,JetMass,JetMassDist,JetMassError,JetMassDistError);
        char MyGraphName5[100];
        sprintf(MyGraphName5,"JetMassDistBin%i_%i",pTHatMin[k],pTHatMax[k]);
        GEJetMass->SetName(MyGraphName5);
        GEJetMass->Write();

        delete GEJet; delete GESingleHadron; delete GEJetShape; delete GEJetFF; delete GEJetMass;
        delete outFile;
    } //k-loop ends here (pTHatBin loop)
    
    delete HistTempJet; delete HistTempSingleHadron;
    delete HistTempJetShape; delete HistTempJetFF; delete HistTempJetMass;
    
    //Done. Script run time
    int EndTime = time(NULL);
    int Hour = (EndTime-StartTime)/3600;
    int Minute = ((EndTime-StartTime)/60)-Hour*60;
    int Second = (EndTime-StartTime)-Hour*60*60 - Minute*60;
    cout<<"Programme run time = "<<Hour<<"::"<<Minute<<"::"<<Second<<endl;
    
    return 0;
}
