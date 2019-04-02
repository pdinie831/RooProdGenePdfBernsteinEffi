//*******************************************************************************
//*                                                                             *
//*  test code to compare efficiency projection and fit parametrization         *
//*  using RooBernsteinClass.                                                   *
//*                                                                             *
//*  P.Dini fecit, Anno Domini MMXVIII                                          *
//*                                                                             *
//*******************************************************************************


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
#include <vector>
#include <math.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH1F.h> 
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>   
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TChain.h>
#include "TBranch.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <TBranch.h>
#include <TApplication.h>
#include <TFile.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAbsReal.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsCategory.h>
#include <RooEffProd.h>
#include <RooProdPdf.h>
#include <RooProduct.h>
#include <RooGenericPdf.h>
#include "RooBernsteinEffi.h"
#include "RooProdGenePdfBernsteinEffi.h"
#include "RooGenePdf.h"
#include <RooMinuit.h>
#include "RooNumIntConfig.h"
#include <TStopwatch.h>
#include <RooWorkspace.h>

using namespace std;
void CreateInputHistoFile(char * OutFileNameInputHisto);
std::map<std::string, std::string>  ReadNamelist(int argc, char** argv);

char InputRecoB0TreeName[10]	= "ntuple";
char OutputRecoB0TreeName[10]	= "ntuple";
char InputFileNameMCReco[300]   = "2016MC_LMNR_addmumukMass.root";
char eosRecoDir[100] = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/";
char eosGeneDir[100] = "/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN/";
// char RecoDir[100]    =  ".";
// char GeneDir[100]    =  ".";
char RecoDir[100]    =  "finalNtuples";
char GeneDir[100]    =  "GEN_NoFilter";


double SetIntPrecisionAbs = 1e-7;
double SetIntPrecisionRel = 1e-7;
int NumNormIntegrals = 11;
//
double Q2Min = 0.; 
double Q2Max = 0.; 
int    Q2Bin = -1;
int    xCosLHBin =  0;
int    xCosKHBin =  0;
int    xPhiHBin  =  0;
double XMinCosThetaL = -1.;
double XMaxCosThetaL = 1.;
double XMinCosThetaK = -1.;
double XMaxCosThetaK = 1.;
double XMinPhi       =-TMath::Pi();
double XMaxPhi       = TMath::Pi();
bool   wrongTagged = false;

char OutFileNameInputHisto[300] =  "";
char OutFileName[300]       =  "";
char SaveFileName[300]      =  "";
char ListParName[300] 	    =  "";
char ListGenValues[300]     =  "";
double xMinQMuMu = 1.;
double xMaxQMuMu = 19.;
double XMinSign = 5.0;
double XMaxSign = 5.6;
double B0Mass	= 5.27962;
double B0Sigma  = 0.030;
double NSigma  = 3.;
double XMinSignW = B0Mass - NSigma*B0Sigma;
double XMaxSignW = B0Mass + NSigma*B0Sigma;
double XStepSign = 0.001;
float xMassHBin  = (XMaxSign -XMinSign)/XStepSign;
float xQ2HBin	 = (xMaxQMuMu -xMinQMuMu)/0.1;
float xMassHBin2  =  xMassHBin /5; // plot only!
int   NCPU = 8;
//
TFile*OutFileInputHisto;
TTree* RecoB0TreeOut=0;
// RooRealVar x("x", "CosL",  XMinCosThetaL,XMaxCosThetaL);
// RooRealVar y("y", "CosK",  XMinCosThetaK,XMaxCosThetaK);
// RooRealVar z("z", "Phi" ,  XMinPhi,XMaxPhi);
// RooDataSet  recoData("recoData", "recoData", RooArgSet(x,y,z));
// RooRealVar* ctL = new RooRealVar("ctL", "ctL",  XMinCosThetaK,XMaxCosThetaK);
// RooRealVar* ctK = new RooRealVar("ctK", "ctK",  XMinCosThetaL,XMaxCosThetaL);
// RooRealVar* phi = new RooRealVar("phi", "phi",  XMinPhi,XMaxPhi);
RooRealVar* ctL=0;
RooRealVar* ctK=0;
RooRealVar* phi=0;
RooDataSet* recoData=0;


// RooRealVar F_L("F_L","F_L", 0.6 , 0.,1.);
// RooRealVar P_1("P_1","P_1",-0.2 ,-1.,1.);
// RooRealVar P_2("P_2","P_2", 0.4 ,-1.,1.);
// RooRealVar P_3("P_3","P_3", 0.0,-1.,1.);
// RooRealVar P4p("P4p","P4p", 0.0,-1.,1.);
// RooRealVar P5p("P5p","P5p", 0.0,-1.,1.);
// RooRealVar P6p("P6p","P6p", 0.0,-1.,1.);
// RooRealVar P8p("P8p","P8p", 0.0,-1.,1.);

// RooRealVar* FlS = new  RooRealVar( "FlS" ,"FlS" , 0.6 , 0.,1.);
// RooRealVar P1S ( "P1S" ,"P1S" ,-0.2 ,-1.,1.);
// RooRealVar P2S ( "P2S" ,"P2S" , 0.4 ,-1.,1.);
// RooRealVar P3S ( "P3S" ,"P3S" , 0.0 ,-1.,1.);
// RooRealVar P4pS("P4pS","P4pS", 0.0 ,-1.,1.);
// RooRealVar P5pS("P5pS","P5pS", 0.0 ,-1.,1.);
// RooRealVar P6pS("P6pS","P6pS", 0.0 ,-1.,1.);
// RooRealVar P8pS("P8pS","P8pS", 0.0 ,-1.,1.);

int main (int argc, char** argv) {
 gROOT ->Reset();
//
 TStopwatch TimeWatch;
 TimeWatch.Start();
//  

//============================
// maxDegree START
// now defined in NAMELIST 
//============================
  int maxDegree1 =0;
  int maxDegree2 =0;
  int maxDegree3 =0;
//============================
// now defined in NAMELIST 
// maxDegree END
//============================

if (argc<=1 ){
    cout<<"Q2Bin not set"<<endl;
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8]\n"<<endl;
    exit(1);
}   


 
switch ( *argv[1] ) {

  case '0' : 
   Q2Min = 1.; 
   Q2Max = 2.;
   Q2Bin = 0;
      xCosLHBin =  25;
      xCosKHBin =  25;
      xPhiHBin  =  25;
    break;
  case '1' : 
   Q2Min = 2.; 
   Q2Max = 4.3; 
   Q2Bin = 1;
   xCosLHBin =  25;
   xCosKHBin =  25;
   xPhiHBin  =  25;
    break;
  case '2' : 
   Q2Min = 4.3; 
   Q2Max = 6.; 
   Q2Bin = 2;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '3' : 
   Q2Min = 6.; 
   Q2Max = 8.68; 
   Q2Bin = 3;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '4' : 
   Q2Min = 8.68; 
   Q2Max = 10.09; 
   Q2Bin = 4;
   xCosLHBin =  25;
   xCosKHBin =  25;
   xPhiHBin  =  25;
    break;
  case '5' : 
   Q2Min = 10.09; 
   Q2Max = 12.86; 
   Q2Bin = 5;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '6' : 
   Q2Min = 12.86; 
   Q2Max = 14.18; 
   Q2Bin = 6;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '7' : 
   Q2Min = 14.18; 
   Q2Max = 16.; 
   Q2Bin = 7;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;
  case '8' : 
   Q2Min = 16; 
   Q2Max = 19.; 
   Q2Bin = 8;
   xCosLHBin = 25;
   xCosKHBin = 25;
   xPhiHBin  = 25;
    break;

  default : 
    // Process for all other cases.
    cout<<"Q2Bin not set correctly!!!"<<endl;
    cout<<"Usage: "<<argv[0]<< " QBin2 [where QBin2=0,1,2,3,4,5,6,7,8]\n"<<endl;
    exit(1);

}
 char NameList[300];;
  
  
  if(wrongTagged){
   sprintf(NameList,"namelist-Effi3DB0-2016-Q2Bin-%d-wrongTagged.lis", Q2Bin);
  }else{
   sprintf(NameList,"namelist-Effi3DB0-2016-Q2Bin-%d-correTagged.lis", Q2Bin);
  }  
  
  char*argn[]={NameList};
  
  std::map<std::string, std::string> mappa = ReadNamelist(1,argn );
//
  maxDegree1     =    atoi (mappa["maxDegree1"].c_str() ) ;
  maxDegree2     =    atoi (mappa["maxDegree2"].c_str() ) ;
  maxDegree3     =    atoi (mappa["maxDegree3"].c_str() ) ;
  xCosLHBin      =    atof (mappa["xCosLHBin" ].c_str() ) ;
  xCosKHBin      =    atof (mappa["xCosKHBin" ].c_str() ) ;
  xPhiHBin       =    atof (mappa["xPhiHBin"  ].c_str() ) ;

  std::cout<<" Num Param Bernstein polynomial CosL:  "<<maxDegree1<<std::endl;
  std::cout<<" Num Param Bernstein polynomial CosK:  "<<maxDegree2<<std::endl;
  std::cout<<" Num Param Bernstein polynomial Phi :  "<<maxDegree3<<std::endl;
  std::cout<<" Binning choice for CosL		:  "<<xCosLHBin<<std::endl;
  std::cout<<" Binning choice for CosK		:  "<<xCosKHBin<<std::endl;
  std::cout<<" Binning choice for Phi		:  "<<xPhiHBin <<std::endl;
  sprintf(OutFileNameInputHisto,"Ntupla-GeneReFit-2016-InputHisto-Q2Bin-%d.root", Q2Bin); 
  sprintf(ListGenValues,"ListGenValues-2016-Q2Bin-%d.txt",Q2Bin);
  std::cout<<"--------------------------------------------\n"<<endl;
  std::cout<<" Setting selection for Q^2 bin: "<<*argv[1]<<" ==> "<<Q2Min<<"<Q^2<"<<Q2Max<<std::endl;
  std::cout<<"--------------------------------------------\n"<<endl;
if (argc>2 && (strcmp(argv[2],"w") == 0) ){
  std::cout<<"====================================="<<endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<" Setting Wrong tagged analysis"<<std::endl;
  std::cout<<"====================================="<<endl;
  sprintf(SaveFileName,"saveEffi-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(OutFileName,"saveEffiGeneReFit-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(ListParName,"ListParValues-2016-strategy-global-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d-wrongTagged.plo", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  wrongTagged = true;
}else{
  std::cout<<"======================================="<<endl;
  std::cout<<"======== DEFAULT: TAGGED =============="<<std::endl;
  std::cout<<"======================================="<<endl;
  sprintf(SaveFileName,"saveEffi-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(OutFileName,"saveEffiGeneReFit-2016-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.root", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  sprintf(ListParName,"ListParValues-2016-strategy-global-Q2Bin-%d-Bins-%d-%d-%d-BernDeg-%d-%d-%d.plo", Q2Bin,xCosLHBin,xCosKHBin,xPhiHBin,maxDegree1,maxDegree2,maxDegree3); 
  wrongTagged = false;
}
if (argc>3 && (strcmp(argv[3],"c") == 0) ){
  std::cout<<"=========================================================================="<<endl;
  std::cout<<"Setting the option: Efficiency function evaluated in the CENTER of the bin"<<std::endl;
  std::cout<<"=========================================================================="<<endl;
//  
  std::string sav(SaveFileName) ;
  sav.replace(sav.end()-5,sav.end(),"-centerBin.root");
  sprintf(SaveFileName,sav.c_str());
//  
  std::string str(OutFileName) ;
  str.replace(str.end()-5,str.end(),"-centerBin.root");
  sprintf(OutFileName,str.c_str());
//  
  std::string plo(ListParName) ;
  plo.replace(plo.end()-4,plo.end(),"-centerBin.plo");
  sprintf(ListParName,plo.c_str());
}else{
  std::cout<<"======================================================================="<<endl;
  std::cout<<"======== DEFAULT: INTEGRAL of the Efficiency in the bin  =============="<<std::endl;
  std::cout<<"======================================================================="<<endl;
//  
  std::string sav(SaveFileName) ;
  sav.replace(sav.end()-5,sav.end(),"-integraBin.root");
  sprintf(SaveFileName,sav.c_str());
//  
  std::string str(OutFileName) ;
  str.replace(str.end()-5,str.end(),"-integraBin.root");
  sprintf(OutFileName,str.c_str());
//  
  std::string plo(ListParName) ;
  plo.replace(plo.end()-4,plo.end(),"-integraBin.plo");
  sprintf(ListParName,plo.c_str());
}
  

  TCanvas* ProjEffiPlots = new TCanvas("ProjEffiPlots","Efficiency projections",200,10,900,780);
  ProjEffiPlots->Divide(2,2);  

  TCanvas* ClosurePlots = new TCanvas("ClosurePlots","Closure test plots ",200,10,900,780);
  ClosurePlots->Divide(2,2);  

  TCanvas* TestEffiPlots = new TCanvas("TestEffiPlots","Efficiency test plots",200,10,900,780);
  TestEffiPlots->Divide(2,2);  

  TCanvas* ClosurePlotsPdf = new TCanvas("ClosurePlotsPdf","Closure test Pdf plots",200,10,900,780);
  ClosurePlotsPdf->Divide(2,2);  

  TCanvas* GenePlotsPdf = new TCanvas("GenePlotsPdf","Gene test Pdf plots",200,10,900,780);
  GenePlotsPdf->Divide(2,2);  
//
  TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
//

  int numParameters = (maxDegree1+1)*(maxDegree2+1)*(maxDegree3+1);
//-------------------------------------------------------

  TFile * SaveFile = new TFile(SaveFileName,"READ");
  if (SaveFile->IsOpen())
  {
    std::cout<<"File:"<<SaveFileName<<" not found!!! Exit..."<<std::endl;
  }else{
  } 
  
  RooWorkspace* ws = (RooWorkspace*)SaveFile->Get("w");
  if ( !ws || ws->IsZombie() ) {
    cout<<Form("Workspace not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } else {
  }
  
  RooBernsteinEffi* Effi = (RooBernsteinEffi*) ws->function("Effi");
  if ( !Effi || Effi->IsZombie() ) {
    cout<<Form("Efficiency not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  ctL = (RooRealVar*) ws->var("ctL");
  if ( !ctL || ctL->IsZombie() ) {
    cout<<Form("ctL not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  ctK = (RooRealVar*) ws->var("ctK");
  if ( !ctK || ctK->IsZombie() ) {
    cout<<Form("ctK not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  phi = (RooRealVar*) ws->var("phi");
  if ( !phi || phi->IsZombie() ) {
    cout<<Form("phi not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* Fl = (RooRealVar*) ws->var("Fl");
  if ( !Fl || Fl->IsZombie() ) {
    cout<<Form("Fl not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P1 = (RooRealVar*) ws->var("P1");
  if ( !P1 || P1->IsZombie() ) {
    cout<<Form("P1 not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P2 = (RooRealVar*) ws->var("P2");
  if ( !P2 || P2->IsZombie() ) {
    cout<<Form("P2S not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P3 = (RooRealVar*) ws->var("P3");
  if ( !P3 || P3->IsZombie() ) {
    cout<<Form("P2S not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P4p = (RooRealVar*) ws->var("P4p");
  if ( !P4p || P4p->IsZombie() ) {
    cout<<Form("P4p not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P5p = (RooRealVar*) ws->var("P5p");
  if ( !P5p || P5p->IsZombie() ) {
    cout<<Form("P5p not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P6p = (RooRealVar*) ws->var("P6p");
  if ( !P6p || P6p->IsZombie() ) {
    cout<<Form("P6p not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  RooRealVar* P8p = (RooRealVar*) ws->var("P8p");
  if ( !P8p || P8p->IsZombie() ) {
    cout<<Form("P8p not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 
  
  recoData = new RooDataSet("recoData", "recoData", RooArgSet(*ctL,*ctK,*phi));
//-------------------------------------------------------
//   RooProdGenePdfBernsteinEffi * ProdPdfEffi = (RooProdGenePdfBernsteinEffi*) ws->function("ProdPdfEffi");
//   if ( !ProdPdfEffi || ProdPdfEffi->IsZombie() ) {
//     cout<<Form("PDF X Efficiency not found in file: %s",SaveFileName)<<endl;
//     exit(1);;
//   } 
  
   RooProdGenePdfBernsteinEffi * recoPdf = (RooProdGenePdfBernsteinEffi*) ws->function("ProdPdfEffi");
  if ( !recoPdf || recoPdf->IsZombie() ) {
    cout<<Form("PDF X Efficiency not found in file: %s",SaveFileName)<<endl;
    exit(1);;
  } 

  if (!TFile::Open(OutFileNameInputHisto,"READ"))
  {
    std::cout<<"File:"<<OutFileNameInputHisto<<" not found!!! create..."<<std::endl;
    CreateInputHistoFile(OutFileNameInputHisto);
  }else{
    std::cout<<"File:"<<OutFileNameInputHisto<<" Found!!!"<<std::endl;
   OutFileInputHisto = TFile::Open(OutFileNameInputHisto,"READ");
   TTree *RecoB0TreeOut    = (TTree*)OutFileInputHisto->Get(OutputRecoB0TreeName);
   if(!RecoB0TreeOut ){
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" not found!!! Suggestion: remove this file e try again..."<<endl;
     exit(1);
   }else{
     cout<<"TTree Reco Data: "<< OutputRecoB0TreeName <<" OK FOUND!!!"<<endl;
     double cos_theta_l    ;
     double cos_theta_k    ;
     double phi_kst_mumu   ;
     RecoB0TreeOut->SetBranchAddress("cos_theta_l"   ,&cos_theta_l );
     RecoB0TreeOut->SetBranchAddress("cos_theta_k"   ,&cos_theta_k );
     RecoB0TreeOut->SetBranchAddress("phi_kst_mumu"  ,&phi_kst_mumu);
//    int nentries = (int)RecoB0TreeOut->GetEntries();
     int nentries = 30000;
     cout<<"nentries: "<< nentries<<endl;
     for (Int_t i=0;i<nentries;i++) {
             RecoB0TreeOut->GetEntry(i);
    	     ctL->setVal(cos_theta_l);
    	     ctK->setVal(cos_theta_k);
    	     phi->setVal(phi_kst_mumu);
//	     std::cout<<ctL->getVal()<<" "<<ctK->getVal()<<" "<<phi->getVal()<<std::endl;
    	     recoData->add(RooArgSet(*ctL,*ctK,*phi));
     }
   }  
  }

 
  
//  RooArgList *coefLis = new RooArgList();
  std::vector<RooRealVar> parLis; 

 
  char varName[10];
  double parIni=0.0;
  std::cout<<"Try to open list of parameters :"<< ListParName <<std::endl;
  std::fstream *parListInput = new std::fstream(ListParName,std::ifstream::in);
  if(parListInput->is_open()){
     std::cout<<"List of initial parameters :"<< ListParName <<" FOUND!!!"<<std::endl;
     for (int i=0;i<numParameters;++i){
      	    sprintf(varName, "p%d", i);
 	    *parListInput >> parIni;
     			 parLis.emplace_back(varName,varName, parIni);
//     			 parLis.emplace_back(varName,varName, parIni,-100.0,100.0);
			 std::cout<<"Coeff ("<<varName<<") = "<<parLis[i].getVal()<<std::endl;
     }
     for (int i=numParameters;i<numParameters+NumNormIntegrals;++i){
	*parListInput >> parIni;
	std::cout<<"Norm["<<i-numParameters<<"] = "<<parIni<<std::endl;
      	                 sprintf(varName, "p%d", i);
     			 parLis.emplace_back(varName,varName, parIni);
     }
    parListInput->close();
  }else{
     std::cout<<"List of initial parameters "<< ListParName <<" not found!"<<std::endl;
     exit(1);
  }   
//   for (int i=0;i<numParameters+NumNormIntegrals;++i){
//    coefLis->add(parLis[i]);
//   }	      
  
//  RooArgList *coefList = new RooArgList();

  std::cout<<"Try to open list of initial parameters :"<< ListGenValues <<std::endl;
  std::fstream *parListInputPdf = new std::fstream(ListGenValues,std::ifstream::in);
  if(parListInputPdf->is_open()){
      char*argn[]={ListGenValues};
  
      std::map<std::string, std::string> mappa = ReadNamelist(1,argn );
      gSystem->Exec(Form("cp %s %s.tmp",ListGenValues,ListGenValues));
      Fl->setVal(atof  (mappa["F_L"].c_str() ) );
      P1->setVal(atof  (mappa["P_1"].c_str() ) );
      P2->setVal(atof  (mappa["P_2"].c_str() ) );
      P3->setVal(atof  (mappa["P_3"].c_str() ) );
      P4p->setVal(atof (mappa["P4p"].c_str() ) );
      P5p->setVal(atof (mappa["P5p"].c_str() ) );
      P6p->setVal(atof (mappa["P6p"].c_str() ) );
      P8p->setVal(atof (mappa["P8p"].c_str() ) );
    parListInputPdf->close();
   }else{
     std::cout<<"List of initial parameters "<< ListGenValues <<" not found, fit with "<<std::endl;
   }   
 
//    coefList->add(*FlS );
//    coefList->add(*P1S );
//    coefList->add(*P2S );
//    coefList->add(*P3S );
//    coefList->add(*P4pS);
//    coefList->add(*P5pS);
//    coefList->add(*P6pS);
//    coefList->add(*P8pS);
   
   
//    RooGenePdf GenePdf("GenePdf","GenePdf",*ctL,*ctK,*phi,*coefList);
//    RooEffProd recoPdf("recoPdf","recoPdf",GenePdf,*Effi);
   
   RooAbsReal::defaultIntegratorConfig()->setEpsAbs(SetIntPrecisionAbs) ;
   RooAbsReal::defaultIntegratorConfig()->setEpsRel(SetIntPrecisionRel) ;
   
   RooAbsReal* nll = recoPdf->createNLL(*recoData,RooFit::NumCPU(NCPU)); 
   RooMinuit Minuit(*nll) ;
//     P3S->setConstant(kTRUE);
//     P4pS->setConstant(kTRUE);
//     P5pS->setConstant(kTRUE);
//     P6pS->setConstant(kTRUE); 
//     P8pS->setConstant(kTRUE);
   Minuit.migrad() ;  

     
   RooPlot* Xframe=ctL->frame(RooFit::Bins(25));
   recoData->plotOn( Xframe, RooFit::MarkerStyle(kStar));
   recoPdf->plotOn(  Xframe, RooFit::LineColor(kRed),RooFit::ProjWData(RooArgSet(*ctL),*recoData));
   ClosurePlotsPdf->cd(1);
   Xframe->Draw();
   RooPlot* Yframe=ctK->frame(RooFit::Bins(25));
   recoData->plotOn( Yframe, RooFit::MarkerStyle(kStar));
   recoPdf->plotOn(  Yframe, RooFit::LineColor(kRed),RooFit::ProjWData(RooArgSet(*ctK),*recoData));
   ClosurePlotsPdf->cd(2);
   Yframe->Draw();
   RooPlot* Zframe=phi->frame(RooFit::Bins(25));
   recoData->plotOn( Zframe, RooFit::MarkerStyle(kStar));
   recoPdf->plotOn(  Zframe, RooFit::LineColor(kRed),RooFit::ProjWData(RooArgSet(*phi),*recoData));
   ClosurePlotsPdf->cd(3);
   Zframe->Draw();
   
   OutFile->cd();
   ClosurePlotsPdf->Write();
   OutFile->Close();
   Xframe->Delete();
   Yframe->Delete();
   Zframe->Delete();

   TimeWatch.Stop();
   TimeWatch.Print();
   
   std::cout<<"======= END OF TEST ======="<<std::endl;
   return 0;
}
//============================================================================
//============================================================================
//============================================================================
void CreateInputHistoFile(char* OutFileNameInputHisto){   
 
  TFile*OutFileNtupla = TFile::Open(OutFileNameInputHisto,"RECREATE");
  RecoB0TreeOut = new TTree(OutputRecoB0TreeName,OutputRecoB0TreeName) ;
  RecoB0TreeOut -> SetAutoSave(500000000);

  
//  TChain* RecoB0Tree = new TChain();

  int nfile = 0;
  
//   nfile = RecoB0Tree->Add(Form("%s/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/%s",RecoDir,InputRecoB0TreeName));  
//   if(nfile==0 ||  !RecoB0Tree->GetFile() ){
//     cout<<"Error:  no Reco files found!!!\n"<<endl;
//     exit(1);
//   }else{
//     printf("Try to open %s/2016MC_RECO_p1p2p3_newtag_LMNR_addW_add4BDT_addvars_bestBDTv4.root/%s \n",RecoDir,InputRecoB0TreeName);
//     cout<<"Opening "<<nfile <<" Reco files found!!!\n"<<endl;
//   }  
//   if(!RecoB0Tree ){
//     cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" not found!!!\n"<<endl;
//     exit(1);
//   }else{
//     cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" OK FOUND!!!\n"<<endl;
//   }  
  TChain* RecoB0Tree = new TChain();  
  nfile = RecoB0Tree->Add(Form("%s/%s/%s",RecoDir,InputFileNameMCReco,InputRecoB0TreeName));
  if( nfile==0 ||  !RecoB0Tree->GetFile() ){
    cout<<"Error:  no Reco files found!!!\n"<<endl;
    exit(1);
  }else{
    printf("Try to open %s/%s/%s\n",RecoDir,InputFileNameMCReco,InputRecoB0TreeName);
    cout<<"Opening "<<nfile <<" Reco files found!!!\n"<<endl;
  }  
  if(!RecoB0Tree ){
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" not found!!!\n"<<endl;
    exit(1);
  }else{
    cout<<"TTree Reco Data: "<< InputRecoB0TreeName <<" OK FOUND!!!\n"<<endl;
  }  

  TH1D* HxRecoMass= new TH1D( "HxRecoMass" , "B^{0} Mass",	xMassHBin2, XMinSign,  XMaxSign);
 
  TH3D* HxReco = new   TH3D( "HxReco"    , "B^{0} Reco correct tagged",  xCosLHBin, XMinCosThetaL, XMaxCosThetaL,
 									 xCosKHBin, XMinCosThetaK, XMaxCosThetaK,
									 xPhiHBin , XMinPhi, XMaxPhi );

//======================================================================
//======================================================================
//======================================================================
//
//			      RECONSTRUCTED EVENTS
//
//======================================================================
//======================================================================
//======================================================================
  double tagged_massF	;
  double mumuMassF	;
  double tagged_mass	;
  double cos_theta_l	;
  double cos_theta_k	;
  double phi_kst_mumu	;
  double mumuMass	;
  double recQ2          ;
  double truthMatchMum  ;
  double truthMatchMup  ;
  double truthMatchTrkm ;
  double truthMatchTrkp ;
  double genSignal	;
  double tagB0  	;
  double  mmk1		;
  double  mmk2		;
  RecoB0Tree->SetBranchAddress("tagged_mass"   ,&tagged_massF);
  RecoB0Tree->SetBranchAddress("cos_theta_l"   ,&cos_theta_l);
  RecoB0Tree->SetBranchAddress("cos_theta_k"   ,&cos_theta_k);
  RecoB0Tree->SetBranchAddress("phi_kst_mumu"  ,&phi_kst_mumu);
  RecoB0Tree->SetBranchAddress("mumuMass"      ,&mumuMassF);
  RecoB0Tree->SetBranchAddress("truthMatchMum" ,&truthMatchMum);
  RecoB0Tree->SetBranchAddress("truthMatchMup" ,&truthMatchMup);
  RecoB0Tree->SetBranchAddress("truthMatchTrkm",&truthMatchTrkm);
  RecoB0Tree->SetBranchAddress("truthMatchTrkp",&truthMatchTrkp);
  RecoB0Tree->SetBranchAddress("tagB0"         ,&tagB0);
  RecoB0Tree->SetBranchAddress("genSignal"     ,&genSignal);
  RecoB0Tree->SetBranchAddress("mmk1"          ,&mmk1);
  RecoB0Tree->SetBranchAddress("mmk2"          ,&mmk2);

//   RecoB0TreeOut->Branch("tagged_mass"   ,&tagged_massF   ,   "tagged_mass"   );
  RecoB0TreeOut->Branch("cos_theta_l"   ,&cos_theta_l    ,   "cos_theta_l/D"   );
  RecoB0TreeOut->Branch("cos_theta_k"   ,&cos_theta_k    ,   "cos_theta_k/D"   );
  RecoB0TreeOut->Branch("phi_kst_mumu"  ,&phi_kst_mumu   ,   "phi_kst_mumu/D"  );
//   RecoB0TreeOut->Branch("mumuMass"      ,&mumuMassF      ,   "mumuMass"	     );
//   RecoB0TreeOut->Branch("truthMatchMum" ,&truthMatchMum  ,   "truthMatchMum" );
//   RecoB0TreeOut->Branch("truthMatchMup" ,&truthMatchMup  ,   "truthMatchMup" );
//   RecoB0TreeOut->Branch("truthMatchTrkm",&truthMatchTrkm ,   "truthMatchTrkm");
//   RecoB0TreeOut->Branch("truthMatchTrkp",&truthMatchTrkp ,   "truthMatchTrkp");
//   RecoB0TreeOut->Branch("tagB0"	        ,&tagB0          ,   "tagB0"	     );
//   RecoB0TreeOut->Branch("genSignal"     ,&genSignal      ,   "genSignal"     );
  
  int nentries = (int)RecoB0Tree->GetEntries();
//  nentries=3000;
  for (Int_t i=0;i<nentries;i++) { 
   RecoB0Tree->GetEntry(i);
   tagged_mass   =  tagged_massF       ;
   mumuMass	 =  mumuMassF          ;  
   recQ2         =  mumuMass*mumuMass  ;
//
   if(truthMatchMum>0.&&truthMatchMup>0.&&truthMatchTrkm>0.&&truthMatchTrkp>0.){
//    HxMass  ->Fill(tagged_mass);
// folding theta_l e phi     
//     cos_theta_l=fabs(cos_theta_l);
//     cos_theta_k=fabs(cos_theta_k);
//       phi_kst_mumu=fabs(phi_kst_mumu);

    if(tagged_mass>=XMinSignW&&tagged_mass<=XMaxSignW){
// first Q^2 Bin
     if(recQ2>Q2Min&&recQ2<Q2Max){
      if(cos_theta_l>=XMinCosThetaL&&cos_theta_l<=XMaxCosThetaL&&cos_theta_k>=XMinCosThetaK&&cos_theta_k<=XMaxCosThetaK){
       if (! ( (mmk1 > 5.158 && mmk1 < 5.398) || (mmk2 > 5.158 && mmk2 < 5.398)) ){
 	HxRecoMass->Fill(tagged_mass);
        if( (tagB0 == 1 && genSignal == 1) || ( tagB0 == 0 && genSignal == 2)){  //(correctly tagged)
 	  HxReco->Fill(cos_theta_l,cos_theta_k,phi_kst_mumu);
	  RecoB0TreeOut->Fill();
	  ctL->setVal(cos_theta_l);
	  ctK->setVal(cos_theta_k);
	  phi->setVal(phi_kst_mumu);
	  recoData->add(RooArgSet(*ctL,*ctK,*phi));
	}    
       }    
      }   
     } 
    }    
   }
  }
  char TXT[200];
  sprintf(TXT,"Mass Reco   Entries = %7f",HxRecoMass->GetEntries());
  std::cout<<"***********************************"<<std::endl;
  std::cout<<"***** RECONSTRUCTED EVENTS ********\n"<<std::endl;
  std::cout<<"***********************************\n"<<std::endl;
  std::cout<<"RecoB0Tree   Entries      = "<<nentries<<std::endl;
  std::cout<<"HxReco       Entries      = "<<HxReco->GetEntries()<<std::endl;
  std::cout<<TXT<<endl;
  std::cout<<"\n***********************************"<<std::endl;
  std::cout<<"***********************************"<<std::endl;
  OutFileNtupla->cd();
  HxReco->Write();
  HxRecoMass->Write();
  RecoB0TreeOut->Write();
  OutFileNtupla->Close();
  return ;  
}
//=========================================================================================
//
// Namelist Routine
//
std::map<std::string, std::string> ReadNamelist(int argc, char** argv){
//    if ( argc>=1 && (strcmp(argv[0],"ListGenValues-2016-Q2Bin-")>=0) ){
//      std::cout<<"Defined namelist: "<<argv[0]<<std::endl;
//    }else{
//      std::cout<<"Namelist:"<<argv[0]<<"  should be named/renamed ListGenValues-2016-Q2Bin-*.txt "<<argc<<std::endl;
//      exit(1);
//    }
   vector<string> split( char *str, char c = ' ');
   ifstream indata;
   std::map<std::string, std::string> mappa;
   std::string line;
   vector<string>vstring ;
//
    indata.open(argv[0]);
   if(!indata) { // file couldn't be opened
   	std::cout <<__LINE__ <<argv[0]<< " Error: fileList can not be opened" << std::endl;
   	exit(1);
   }
   while(std::getline(indata, line)) {
	 line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
	 line.erase(std::remove(line.begin(), line.end(), ' ' ), line.end());

 	 char *cstr = new char [line.size()+1];


 	 strcpy (cstr, line.c_str());
//	 cout <<"stringa->"<< cstr << endl;
	 vstring = split(cstr,'=');
	 mappa.insert( std::pair<string,string>(vstring[0],vstring[1]) );
    }
    std::cout<<"#########################################################################################################"<<std::endl;
    for (map<string,string>::iterator imap = mappa.begin();
    			       imap != mappa.end();
    			       ++imap)
    {
   	std::cout <<"mappa->"<< (*imap).first<<" = "<<(*imap).second << std::endl;
    }
    std::cout<<"#########################################################################################################"<<std::endl;
    indata.close();	
  return mappa ;
}
vector<string> split( char *str, char c = ' ')
{
    vector<string> result;

    while(1)
    {
         char *begin = str;

        while(*str != c && *str)
                str++;

        result.push_back(string(begin, str));

        if(0 == *str++)
                break;
    }

    return result;
}
