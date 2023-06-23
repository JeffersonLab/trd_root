#define trdclass_cxx
#include "trdclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "PlotLib.C"

#define NPRT 1000
//#define NPRT 10
#define USE_TRK
#define  MAX_PRINT 1

// -- GEMTRD mapping --
int GetGEMChan(int ch, int slot) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  int invCardChannel = 23-cardChannel;
  if (slot<6 || (slot==6 && ch<24)) {
    return invCardChannel+cardNumber*24+(slot-3)*72.;
  }
  return -1;
}

// -- MMG-1 mapping --
int GetMMG1Chan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  int invCardChannel = 23-cardChannel;
  float dchan = invCardChannel+cardNumber*24+(slot-3)*72.;
	
  if (runNum<3148) { // -- Map #1
    if (slot==7 || slot==8 || (slot==6&&ch>23) || (slot==9&&ch<48)) {
      return dchan - 240.;
    }
  } else if (runNum>3147 && runNum<3262) { // -- Map #2
    if (slot==8 || (slot==9&&ch<48)) {
      return dchan - 360.;
    }
  } else if (runNum>3261 && runNum<3279) { // -- Map #3
    if (slot==7&&ch>23) {
      return dchan - 312.;
    }
    if (slot==8 || (slot==9&&ch<48)) {
      return dchan - 240.;
    }
  } else if (runNum>3278) { // -- Map #4
    if (slot==8 || (slot==7&&ch>23) || (slot==9&&ch<48)) {
      return dchan - 264.;
    }
  }
  return -1;
}

// -- MMG-2 mapping --
int GetMMG2Chan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  int invCardChannel = 23 - cardChannel;
  float dchan = invCardChannel+cardNumber*24+(slot-3)*72.;
	
  if (runNum>3261 && runNum<3279) { // -- Map #3
    if (slot==6&&ch>23&&ch<48) {
      return dchan - 144.;
    }
    if (slot==6&&ch>47) {
      return dchan - 192.;
    }
    if (slot==7&&ch<24) {
      return dchan - 240.;
    }
  } else if (runNum>3278) { // -- Map #4
    if (slot==6&&ch>23&&ch<48) {
      return dchan - 144.;
    }
    if (slot==6&&ch>47) {
      return dchan - 192.;
    }
    if (slot==7&&ch<24) {
      return dchan - 240.;
    }
  }
  return -1;
}

// -- uRWELLTRD mapping --
int GetRWELLChan(int ch, int slot, int runNum) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  int invCardChannel = 23 - cardChannel;
  float dchan = invCardChannel+cardNumber*24+(slot-3)*72.;
	
  if (runNum<3148) { // -- Map #1
    if (slot==10 || (slot==9&&ch>47)) {
      return dchan - 480.;
    }
    if (slot==13 || slot==14) {
      float specialChan = invCardChannel+cardNumber*24+(slot-5)*72.;
      return specialChan - 480.;
    }
  } else if (runNum>3147 && runNum<3262) { // -- Map #2
    if (slot==7 || (slot==6&&ch>23)) {
      return dchan - 240.;
    }
  }
  return -1;
}

void trdclass::Loop() {

  //   In a ROOT session, you can do:
  //      root> .L trdclass.C
  //      root> trdclass t(RunNum)
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
  //

  //    This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  	Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //    To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //	by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;

  //========= Book Histograms =============

   TList *HistList = new TList();

  //-----------------  canvas 0 Event Display ----------
  char c0Title[256]; sprintf(c0Title,"Event_Display_Run=%d",RunNum);
  TCanvas *c0 = new TCanvas("DISP",c0Title,200,200,1500,1300);
  c0->Divide(3,3); c0->cd(1);

  TF1 fx("fx","pol1",100,190);
  TF1 fx1("fx1","pol1",100,190);
  TF1 fx2("fx2","pol1",100,190);

  //-- GEMTRD - GEMTRK alignment --------
  double xx1=-40., yy1=-55.,  xx2=55., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1); 
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  //-------------------------------------

  hCal_occ = new TH1F("hCal_occ"," Calorimeter Occupancy",9,-0.5,8.5);         HistList->Add(hCal_occ);
  hCal_sum = new TH1F("hCal_sum"," Calorimeter Sum (GeV)",100.,0.,15.);        HistList->Add(hCal_sum);
  for (int cc=0; cc<NCAL; cc++) {
    char hName[128];  sprintf(hName,"hCal_adc%d",cc);
    char hTitle[128]; sprintf(hTitle,"Calorimeter ADC, cell%d",cc);
    hCal_adc[cc] = new TH1F(hName,hTitle,4096,-0.5,4095.5);                    HistList->Add(hCal_adc[cc]);
    //
    sprintf(hName,"hCal_cor%d",cc);  sprintf(hTitle,"Calorimeter ADC corr GemTrd X, cell%d",cc);   
    hCal_cor[cc] =  new TH2F(hName,hTitle,25,0.5,25.5,240,0.5,240.5);          HistList->Add(hCal_cor[cc]);
    //
    sprintf(hName,"hCal_trk%d",cc);  sprintf(hTitle,"Calorimeter corr GemTRK, cell%d",cc);   
    hCal_trk[cc] =  new TH2F(hName,hTitle,100,-55.,55.,100,-55.,55.);          HistList->Add(hCal_trk[cc]);
    //
    sprintf(hName,"hCal_cal%d",cc);  sprintf(hTitle,"Calorimeter ADC Calib, cell%d",cc);   
    hCal_cal[cc] = new TH2F(hName,hTitle,10,-0.5,4095.5,10,-5.,15.);           HistList->Add(hCal_cal[cc]);
  }

  hcount= new TH1D("hcount","Count",3,0,3);                                     HistList->Add(hcount);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
#if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
  hcount->SetCanExtend(TH1::kXaxis);
#else
  hcount->SetBit(TH1::kCanRebin);
#endif

  h250_size = new TH1F("h250_size"," fa250 Raw data size",4096,0.5,4095.5);                                        HistList->Add(h250_size); 
	
  hCal_sum_el = new TH1F("hCal_sum_el"," Calorimeter Sum for electrons",100,0.,25.);                               HistList->Add(hCal_sum_el);
  hCal_sum_pi = new TH1F("hCal_sum_pi"," Calorimeter Sum for pions",100,0.,25.);                                   HistList->Add(hCal_sum_pi);
  hCher_u_adc = new TH1F("hCher_u_adc"," Cherenkov Upstream ADC ; ADC Amplitude ",4096,-0.5,4095.5);               HistList->Add(hCher_u_adc);
  hCher_din_adc = new TH1F("hCher_din_adc"," Cherenkov Downstream (in) ADC ; ADC Amplitude ",4096,-0.5,4095.5);       HistList->Add(hCher_din_adc);
  hCher_dout_adc = new TH1F("hCher_dout_adc"," Cherenkov Downstream (out) ADC ; ADC Amplitude ",4096,-0.5,4095.5);    HistList->Add(hCher_dout_adc);
  
  hCher_u_time = new TH1F("hCher_u_time"," Cherenkov Upstream Time Response ",300,-0.5,299.5);                         HistList->Add(hCher_u_time);
  hCher_din_time = new TH1F("hCher_din_time"," Cherenkov Downstream (in) Time Response ",300,-0.5,299.5);              HistList->Add(hCher_din_time);
  hCher_dout_time = new TH1F("hCher_dout_time"," Cherenkov Downstream (out) Time Response ",300,-0.5,299.5);           HistList->Add(hCher_dout_time);
  
  hCCCor_u = new TH2F("hCCCor_u"," Cherenkov Calorimeter Corr ; Upstream ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);                 HistList->Add(hCCCor_u);
  hCCCor_dout = new TH2F("hCCCor_dout"," Cherenkov Calorimeter Corr ; Downstream (out) ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);   HistList->Add(hCCCor_dout);
  hCCor_ud = new TH2F("hCCor_ud"," Cherenkov Upstream/Downstream Corr ; Upstream ; Downstream (out) ",400,-0.5,4095.5,400,0.5,4095.5);   HistList->Add(hCCor_ud);
  
  f125_el = new TH1F("f125_el","GEM-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);                  HistList->Add(f125_el);
  f125_pi = new TH1F("f125_pi","GEM-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);                      HistList->Add(f125_pi);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg1_f125_el);
  mmg1_f125_pi = new TH1F("mmg1_f125_pi","MMG1-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_pi);
  urw_f125_el = new TH1F("urw_f125_el","uRW-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);          HistList->Add(urw_f125_el);
  urw_f125_pi = new TH1F("urw_f125_pi","uRW-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);              HistList->Add(urw_f125_pi);
  mmg2_f125_el = new TH1F("mmg2_f125_el","MMG2-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg2_f125_el);
  mmg2_f125_pi = new TH1F("mmg2_f125_pi","MMG2-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg2_f125_pi);

  f125_el_chi2 = new  TH1F("f125_el_chi2","TRK_el_chi2",100,0.,10000.);                            HistList->Add(f125_el_chi2);
  f125_pi_chi2 = new  TH1F("f125_pi_chi2","TRK_pi_chi2",100,0.,10000.);                            HistList->Add(f125_pi_chi2);
  
  //---- GEM-TRK --


  srs_ncl = new TH1F("srs_ncl"," Number SRS clusters per event",10,-0.5,9.5);                     HistList->Add(srs_ncl); 
  srs_trk_el = new TH2F("srs_trk_el","GEM-TRK , Electrons ; X ; Y ",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_trk_el);
  //srs_trk_el_->SetStats(0);  srs_trd_el->SetMinimum(THRESH);  srs_trk_el->SetMaximum(1000.);
  srs_trk_pi = new TH2F("srs_trk_pi","GEM-TRK , Pions ; X ; Y ",100,-55.,55.,100,-55.,55.);        HistList->Add(srs_trk_pi);
  //srs_trk_pi->SetStats(0);  srs_trk_pi->SetMinimum(THRESH);  srs_trk_pi->SetMaximum(1000.);

  srs_gem_x = new TH2F("srs_gem_x","Correlation TRD-TRK X ; X trk; X chan",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_gem_x);
  srs_gem_y = new TH2F("srs_gem_y","Correlation TRD-TRK Y ; Y trk; X chan",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_gem_y);

  srs_cal_corr = new TH2F("srs_cal_corr","Correlation TRK-CAL; X ; Y ",100,-55.,55.,100,-55.,55.);            HistList->Add(srs_cal_corr);
  srs_etrd_corr = new TH2F("srs_etrd_corr","Correlation TRK-energy TRD; X ; Y ",100,-55.,55.,100,-55.,55.);   HistList->Add(srs_etrd_corr);

  //---- GEM-TRD --
  int THRESH=120;
  f125_el_evt = new TH2F("f125_el_evt","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);  HistList->Add(f125_el_evt);
  f125_el_evt->SetStats(0);  f125_el_evt->SetMinimum(THRESH);  f125_el_evt->SetMaximum(1000.);
  f125_pi_evt = new TH2F("f125_pi_evt","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);      HistList->Add(f125_pi_evt);
  f125_pi_evt->SetStats(0); f125_pi_evt->SetMinimum(THRESH); f125_pi_evt->SetMaximum(1000.);
  f125_el_raw = new TH2F("f125_el_raw","GEM-TRD raw for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);    HistList->Add(f125_el_raw);
  f125_el_raw->SetStats(0);  f125_el_raw->SetMinimum(THRESH);   f125_el_raw->SetMaximum(1000.);
  f125_pi_raw = new TH2F("f125_pi_raw","GEM-TRD raw for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);        HistList->Add(f125_pi_raw);
  f125_pi_raw->SetStats(0); f125_pi_raw->SetMinimum(THRESH); f125_pi_raw->SetMaximum(1000.);
  
  f125_el_fit = new TH2F("f125_el_fit","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);     HistList->Add(f125_el_fit);
  f125_pi_fit = new TH2F("f125_pi_fit","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);         HistList->Add(f125_pi_fit);

  f125_el_amp2ds = new TH2F("f125_el_amp2ds","GEM-TRD Amp for Single Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);  HistList->Add(f125_el_amp2ds);
  f125_pi_amp2ds = new TH2F("f125_pi_amp2ds","GEM-TRD Amp for Single Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);      HistList->Add(f125_pi_amp2ds);
  
  f125_el_amp2d = new TH2F("f125_el_amp2d","GEM-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);           HistList->Add(f125_el_amp2d);
  f125_pi_amp2d = new TH2F("f125_pi_amp2d","GEM-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);               HistList->Add(f125_pi_amp2d);

  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);   HistList->Add(mmg1_f125_el_amp2d);
  mmg1_f125_pi_amp2d = new TH2F("mmg1_f125_pi_amp2d","MMG1-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);       HistList->Add(mmg1_f125_pi_amp2d);
  mmg2_f125_el_amp2d = new TH2F("mmg2_f125_el_amp2d","MMG2-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);   HistList->Add(mmg2_f125_el_amp2d);
  mmg2_f125_pi_amp2d = new TH2F("mmg2_f125_pi_amp2d","MMG2-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);       HistList->Add(mmg2_f125_pi_amp2d);
  urw_f125_el_amp2d = new TH2F("urw_f125_el_amp2d","uRW-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);      HistList->Add(urw_f125_el_amp2d);
  urw_f125_pi_amp2d = new TH2F("urw_f125_pi_amp2d","uRW-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);          HistList->Add(urw_f125_pi_amp2d);
  
  f125_el_clu2d = new TH2F("f125_el_clu2d","GEM-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);               HistList->Add(f125_el_clu2d);
  f125_pi_clu2d = new TH2F("f125_pi_clu2d","GEM-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);                   HistList->Add(f125_pi_clu2d);
  mmg1_f125_el_clu2d = new TH2F("mmg1_f125_el_clu2d","MMG1-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);    HistList->Add(mmg1_f125_el_clu2d);
  mmg1_f125_pi_clu2d = new TH2F("mmg1_f125_pi_clu2d","MMG1-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);        HistList->Add(mmg1_f125_pi_clu2d);
  mmg2_f125_el_clu2d = new TH2F("mmg2_f125_el_clu2d","MMG2-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);    HistList->Add(mmg2_f125_el_clu2d);
  mmg2_f125_pi_clu2d = new TH2F("mmg2_f125_pi_clu2d","MMG2-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);        HistList->Add(mmg2_f125_pi_clu2d);
  urw_f125_el_clu2d = new TH2F("urw_f125_el_clu2d","uRW-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);       HistList->Add(urw_f125_el_clu2d);
  urw_f125_pi_clu2d = new TH2F("urw_f125_pi_clu2d","uRW-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);           HistList->Add(urw_f125_pi_clu2d);
 
  //-------------------------------------------------------------------------

  //--- Calorimeters Calibration ----

  //-------Cal Cell        0      1      2      3      4      5      6 
  double CalCal10[NCAL]={ 4096., 3300., 3700., 2100., 2350., 2400., 1760. };
  double CalCal3[NCAL] ={ 1700., 1180., 1340.,  820.,  820.,  860.,  660. };
  double ac[NCAL],ab[NCAL];
  TF1 *fcal[NCAL]; 
  for (int i=0; i<NCAL; i++) {
    ac[i]=(10.-3.)/(CalCal10[i]-CalCal3[i]);
    ab[i]=3.-ac[i]*CalCal3[i];
    char fnam[64]; sprintf(fnam,"fcal%d",i);
    fcal[i] = new TF1(fnam,"[0]*x+[1]",-5.,4000.); 
    fcal[i]->SetParameter(0,ac[i]);
    fcal[i]->SetParameter(1,ab[i]);
  }
  double Ebeam=10.; // GeV 

  if  (3131 <= RunNum && RunNum <= 3171) {    Ebeam=10.;   }   // NEG 10 GeV
  if  (3172 <= RunNum && RunNum <= 3176) {    Ebeam=120.;   }   // protons
  if  (3177 <= RunNum && RunNum <= 3210) {    Ebeam=10.;   }   //  NEG 10 GeV
  if  (3211 <= RunNum && RunNum <= 3219) {    Ebeam=3.;    }   // NEG 3 GeV
  if  (3220 <= RunNum && RunNum <= 3226) {    Ebeam=120.;  }   // protons
  if  (3227 <= RunNum && RunNum <= 3227) {    Ebeam=10.;   }   // POS 10 GeV 
  if  (3228 <= RunNum && RunNum <= 3240) {    Ebeam=120.;  }   // protons 
  if  (3241 <= RunNum && RunNum <= 3253) {    Ebeam=3.;    }   // NEG 3GeV
  if  (3255 <= RunNum && RunNum <= 3261) {    Ebeam=120.;  }   // protons
  if  (3272 <= RunNum && RunNum <= 3288) {    Ebeam=10.;   }   // NEG 10 GeV
  if  (3289 <= RunNum && RunNum <= 3290) {    Ebeam=120.;  }   // protons

  double Ebeam_el=0.3*Ebeam;
  double Ebeam_pi=0.1*Ebeam; 

  //=========================================

  gem_zHist = new  TH1F("gem_zHist", "gem_zHist", 20, 80., 200.);
  mmg1_zHist = new  TH1F("mmg1_zHist", "mmg1_zHist", 20, 80., 200.);
  mmg2_zHist = new  TH1F("mmg2_zHist", "mmg2_zHist", 20, 80., 200.);
  urw_zHist = new  TH1F("urw_zHist", "urw_zHist", 20, 80., 200.);
  
  //=========================================
  TFile* fHits;
  int save_hits_root = 1;
  
  if (save_hits_root) {
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/trd_singleTrackHits_Run_%06d.root", RunNum);
    fHits = new TFile(hitsFileName, "RECREATE");
	  
    //--- Vector Branches -----
    EVENT_VECT_GEM = new TTree("gem_hits","GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&gem_nhit,"gem_nhit/I");
    //EVENT_VECT_GEM->Branch("trackID",&gem_trackID);
    EVENT_VECT_GEM->Branch("xpos",&gem_xpos);
    //EVENT_VECT_GEM->Branch("ypos",&gem_ypos);
    EVENT_VECT_GEM->Branch("zpos",&gem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&gem_dedx);
    EVENT_VECT_GEM->Branch("parID",&gem_parID);
    EVENT_VECT_GEM->Branch("zHist",&gem_zHist_vect);
	  
    EVENT_VECT_MMG1 = new TTree("mmg1_hits","MMG1 TTree with single track hit info");
    EVENT_VECT_MMG1->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG1->Branch("nhit",&mmg1_nhit,"mmg1_nhit/I");
    //EVENT_VECT_MMG1->Branch("trackID",&mmg1_trackID);
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    //EVENT_VECT_MMG1->Branch("ypos",&mmg1_ypos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
	  
    EVENT_VECT_MMG2 = new TTree("mmg2_hits","MMG2 TTree with single track hit info");
    EVENT_VECT_MMG2->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG2->Branch("nhit",&mmg2_nhit,"mmg2_nhit/I");
    //EVENT_VECT_MMG2->Branch("trackID",&mmg2_trackID);
    EVENT_VECT_MMG2->Branch("xpos",&mmg2_xpos);
    //EVENT_VECT_MMG2->Branch("ypos",&mmg2_ypos);
    EVENT_VECT_MMG2->Branch("zpos",&mmg2_zpos);
    EVENT_VECT_MMG2->Branch("dedx",&mmg2_dedx);
    EVENT_VECT_MMG2->Branch("parID",&mmg2_parID);
    EVENT_VECT_MMG2->Branch("zHist",&mmg2_zHist_vect);
	  
    EVENT_VECT_URW = new TTree("urw_hits","uRWELL TTree with single track hit info");
    EVENT_VECT_URW->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_URW->Branch("nhit",&urw_nhit,"urw_nhit/I");
    //EVENT_VECT_URW->Branch("trackID",&urw_trackID);
    EVENT_VECT_URW->Branch("xpos",&urw_xpos);
    //EVENT_VECT_URW->Branch("ypos",&urw_ypos);
    EVENT_VECT_URW->Branch("zpos",&urw_zpos);
    EVENT_VECT_URW->Branch("dedx",&urw_dedx);
    EVENT_VECT_URW->Branch("parID",&urw_parID);
    EVENT_VECT_URW->Branch("zHist",&urw_zHist_vect);
  }

  //=========================================

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  //int MAX_PRINT=1; //--- debug printing ---

  //================ Begin Event Loop ==============

  int N_trk_el=0;
  int N_trk_pi=0;
  Long64_t jentry=0;
  
  for (jentry=0; jentry<nentries; jentry++) { //-- Event Loop --
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
	  
    if (jentry<MAX_PRINT || !(jentry%NPRT))
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_raw_count=%llu %llu %llu %llu \n"
	     ,jentry,f125_wraw_count, f125_pulse_count, f250_wraw_count
	     ,srs_raw_count, gem_scluster_count, srs_prerecon_count, gem_peak_count);
    
    event_num = jentry;
	  
    gem_nhit=0;
    gem_xpos.clear();
    //gem_ypos.clear();
    gem_zpos.clear();
    gem_dedx.clear();
    gem_parID.clear();
    gem_zHist->Reset();
    gem_zHist_vect.clear();
	
    mmg1_nhit=0;
    mmg1_xpos.clear();
    //mmg1_ypos.clear();
    mmg1_zpos.clear();
    mmg1_dedx.clear();
    mmg1_parID.clear();
    mmg1_zHist->Reset();
    mmg1_zHist_vect.clear();
	 
    mmg2_nhit=0;
    mmg2_xpos.clear();
    //mmg2_ypos.clear();
    mmg2_zpos.clear();
    mmg2_dedx.clear();
    mmg2_parID.clear();
    mmg2_zHist->Reset();
    mmg2_zHist_vect.clear();
	
    urw_nhit=0;
    urw_xpos.clear();
    //urw_ypos.clear();
    urw_zpos.clear();
    urw_dedx.clear();
    urw_parID.clear();
    urw_zHist->Reset();
    urw_zHist_vect.clear();
	  
    //==================================================================================================
    //                    Show Event
    //==================================================================================================
	
    if (jentry<MAX_PRINT) {

      printf("-------------------- Pulse  125  ---------------------------\n");
      for (ULong64_t i=0;i<f125_pulse_count; i++) {
	printf("F125:pulse: i=%lld  sl=%d, ch=%d, npk=%d time=%d amp=%d ped=%d \n"
	       ,i,f125_pulse_slot->at(i),f125_pulse_channel->at(i),f125_pulse_npk->at(i)
	       ,f125_pulse_peak_time->at(i),f125_pulse_peak_amp->at(i),f125_pulse_pedestal->at(i));
      }
	   
      printf("-------------------- Pulse  250  n=%lld ---------------------------\n",f250_pulse_count);
      for (ULong64_t i=0;i<f250_pulse_count; i++) {
	printf("F250:: i=%lld  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n"
	       ,i,f250_pulse_slot->at(i),f250_pulse_channel->at(i),f250_pulse_pulse_number->at(i)
	       ,f250_pulse_course_time->at(i),f250_pulse_pulse_peak->at(i),f250_pulse_pedestal->at(i)/4.);
      }
      //if(f250_pulse_count>2) sleep(3);

      printf("-------------------- Raw  125  ---------------------------\n");
      for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
	printf("F125:raw: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
	       ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
	       ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
      }

      printf("-------------------- SRS     ---------------------------\n");
      for (ULong64_t i=0;i<gem_scluster_count; i++) { // --- SRS cluster loop
	printf("SRS:clusters:  i=%lld  X=%f Y=%f  \n",i,gem_scluster_x->at(i), gem_scluster_y->at(i));	      
      }
      //if (gem_scluster_count>0) sleep(3);
   }



    //==================================================================================================
    //                    Process Fa250  Pulse data
    //==================================================================================================

    //#define USE_250_PULSE
#ifdef USE_250_PULSE

    printf("-------------------- Pulse  250  n=%lld ---------------------------\n",f250_pulse_count);
    for (int i=0;i<f250_pulse_count; i++) {
      printf("F250:: i=%d  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n"
	     ,i,f250_pulse_slot->at(i),f250_pulse_channel->at(i),f250_pulse_pulse_number->at(i)
	     ,f250_pulse_course_time->at(i),f250_pulse_pulse_peak->at(i),f250_pulse_pedestal->at(i)/4.);
    } 
#endif
    //==================================================================================================
    //                    Process Fa250  RAW data
    //==================================================================================================

    if (jentry<MAX_PRINT) printf("------------------ Fadc250  wraw_count = %llu ---------\n", f250_wraw_count);
    h250_size->Fill(f250_wraw_count);
	
    double CalSum=0;
    double Ch_u=0;
    double Ch_in=0;
    double Ch_out=0;
    bool electron_ch=false;
    bool electron=false;
    bool pion=false;
    double Ecal[NCAL]; for (int i=0; i<NCAL; i++) Ecal[i]=0;
	
    for (ULong64_t i=0; i<f250_wraw_count; i++) { // --- fadc250 channels loop
      if (jentry<MAX_PRINT) printf("F250:: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
				   ,i,f250_wraw_slot->at(i),f250_wraw_channel->at(i)
				   ,f250_wraw_samples_index->at(i),f250_wraw_samples_count->at(i));
	    
      int fadc_chan = f250_wraw_channel->at(i);
      int fadc_window = f250_wraw_samples_count->at(i);
      hCal_occ->Fill(fadc_chan+0.);
		
      //printf("f250Loop:: fadc_window=%d\n",fadc_window);
      //if (fadc_window>5000) break;
      //printf("f250Loop:: fadc_window=%d  take it \n",fadc_window);
		
      int amax=0;
      int tmax=9;
      for (int si=0; si<fadc_window; si++) {
	//printf("f250Loop:: %d fadc_window=%d\n",si,fadc_window);
	int adc = f250_wraw_samples->at(f250_wraw_samples_index->at(i)+si); // printf(" sample=%d adc=%d \n",si,adc);
	if (adc>amax) {
	  amax=adc;
	  tmax=si;
	}
      } // --  end of samples loop
      if (fadc_chan<NCAL) {
	Ecal[fadc_chan]=fcal[fadc_chan]->Eval(amax);
	hCal_adc[fadc_chan]->Fill(amax);
	CalSum+=Ecal[fadc_chan];
      } else { // Cherenkov
	if (fadc_chan==13) { hCher_u_adc->Fill(amax);   hCher_u_time->Fill(tmax); Ch_u=amax; }
	if (fadc_chan==14) { hCher_din_adc->Fill(amax);  hCher_din_time->Fill(tmax); Ch_in=amax; }
	if (fadc_chan==15) { if(amax>300)electron_ch=true; hCher_dout_adc->Fill(amax);  hCher_dout_time->Fill(tmax);Ch_out=amax; Count("eCHR");}
      }
    } // -- end of fadc250 channels loop


    //=========================  set PID : Rich  + eCAL   =========

    if (electron_ch  && CalSum > Ebeam_el) { electron=true;  Count("elCC"); }
    if ( !electron_ch && CalSum < Ebeam_pi){ pion=true;  Count("piCC"); }

    //=======================  End Fa250 RAW  process Loop  =====================================================
	
    hCal_sum->Fill(CalSum);
    if(electron_ch){
      hCal_sum_el->Fill(CalSum);
    } else {
      hCal_sum_pi->Fill(CalSum);
    }
	
    //==================================================================================================
    //                    Process SRS data
    //==================================================================================================
    srs_ncl->Fill(gem_scluster_count);
    for (ULong64_t i=0;i<gem_scluster_count; i++) { // --- SRS cluster loop
      if (i==0) Count("nclSRS");
      double gemtrk_x=0., gemtrk_y=0., x, y;
      x=gem_scluster_x->at(i); if (x<=0) gemtrk_x=x+50.; else gemtrk_x=x-50.; gemtrk_x*=-1.;
      y=gem_scluster_y->at(i); if (y<=0) gemtrk_y=y+50.; else gemtrk_y=y-50.; gemtrk_y*=-1.;
      double gemtrk_E=gem_scluster_energy->at(i);
      //printf("SRS:clusters:  i=%lld  X=%f Y=%f E=%f  \n",i,gemtrk_x, gemtrk_y,gemtrk_E);
      if (electron_ch)      srs_trk_el->Fill(gemtrk_x,gemtrk_y,gemtrk_E);
      else                  srs_trk_pi->Fill(gemtrk_x,gemtrk_y,gemtrk_E);


      for (int cc=0; cc<NCAL; cc++) {
	  if (Ecal[cc]>0.8*Ebeam)  {
	    srs_cal_corr->Fill(gemtrk_x,gemtrk_y);
	    hCal_trk[cc]->Fill(gemtrk_x,gemtrk_y);
	  }
	}
    }
   
    //==================================================================================================
    //                    Process Fa125  Pulse  data
    //==================================================================================================
    if (!(jentry%NPRT)) {
      if(electron) {
	f125_el_evt->Reset();
	f125_el_raw->Reset();
      } else if (pion) {
	f125_pi_evt->Reset();
	f125_pi_raw->Reset();
      }
    }
    if (electron) {
      f125_el_fit->Reset();
      f125_el_amp2ds->Reset();
    } else {
      f125_pi_fit->Reset();
      f125_pi_amp2ds->Reset();
    }
	
    //==================================================================================================
    //                    Single Track Event Fitting
    //==================================================================================================
	
    int gem_chan_max=-1;
    //int mmg_chan_max=-1;
    //int urw_chan_max=-1;
    float f125_amp_max=0.;
    //float mmg_f125_amp_max=0.;
    //float urw_f125_amp_max=0.;
	
    for (ULong64_t i=0;i<f125_pulse_count; i++) {
      if (jentry<MAX_PRINT) printf("F125:: i=%lld  sl=%d, ch=%d, npk=%d time=%d amp=%d ped=%d \n"
				   ,i,f125_pulse_slot->at(i),f125_pulse_channel->at(i),f125_pulse_npk->at(i)
				   ,f125_pulse_peak_time->at(i),f125_pulse_peak_amp->at(i),f125_pulse_pedestal->at(i));
      //cout<<" ++++++++++++++++++++ f125_pulse_npk= "<<f125_pulse_npk->at(i)<<endl;
	
		
      //------ TRACK FITTING FILL ------
      float peak_amp = f125_pulse_peak_amp->at(i);
      float ped = f125_pulse_pedestal->at(i);
      if (0 > ped || ped > 200 ) ped = 100;
      float amp=peak_amp-ped;
      float time=f125_pulse_peak_time->at(i);   int itime=f125_pulse_peak_time->at(i);
      int fADCSlot = f125_pulse_slot->at(i);
      int fADCChan = f125_pulse_channel->at(i);
      int gemChan = GetGEMChan(fADCChan, fADCSlot);
      int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
      int mmg2Chan = GetMMG2Chan(fADCChan, fADCSlot, RunNum);
      int rwellChan = GetRWELLChan(fADCChan, fADCSlot, RunNum);
	
      if (amp<0) amp=0;
      int MM_THR=50;

      if(electron) {  //-----------------  electron ---------------
	if (gemChan>-1) {
#ifdef USE_TRK
	  f125_el_amp2ds->Fill(time,gemChan,amp);
	  f125_el_fit->Fill(time,gemChan,amp);
#else
	  f125_el_amp2d->Fill(time,gemChan,amp);
#endif
	}
      } else if (pion) {   //--------------- hadron / pion ---------------
	if (gemChan>-1) {
#ifdef USE_TRK
	  f125_pi_amp2ds->Fill(time,gemChan,amp);
	  f125_pi_fit->Fill(time,gemChan,amp);
#else
	  f125_pi_amp2d->Fill(time,gemChan,amp);
#endif
	}
      }
    } //-- End f125 pulse loop --

    double x0=0;
    double chi2_max=20000;
#ifdef USE_TRK
    //--- Calorimeter Correlation ---
    double chi2cc = TrkFit(f125_el_fit,fx,"fx");
    double a = fx.GetParameter(1);
    double b = fx.GetParameter(0);
    if (chi2cc>0. && chi2cc<chi2_max)  {
      x0=fx.Eval(100.)*0.4-50.;  // to [mm] ;
      for (ULong64_t i=0;i<gem_scluster_count; i++) {   // --- SRS cluster loop
	double gemtrk_x=0., gemtrk_y=0., x, y;
	x=gem_scluster_x->at(i); if (x<=0) gemtrk_x=x+50.; else gemtrk_x=x-50.; gemtrk_x*=-1.;
	y=gem_scluster_y->at(i); if (y<=0) gemtrk_y=y+50.; else gemtrk_y=y-50.; gemtrk_y*=-1.;
	double gemtrk_E=gem_scluster_energy->at(i);
	srs_gem_x->Fill(gemtrk_x,x0);
	srs_gem_y->Fill(gemtrk_y,x0);
      }
      //---------------------------------
      for (int cc=0; cc<NCAL; cc++) {   //--- Calorimeter cells ---
	if (Ecal[cc]>0.8*Ebeam) { 
	  if (electron_ch) { hCal_cor[cc]->Fill(5.,x0);} else { hCal_cor[cc]->Fill(7.,x0); }
	} else if (Ecal[cc]<0.1*Ebeam) {
	  if (electron_ch) { hCal_cor[cc]->Fill(10.,x0);}  else { hCal_cor[cc]->Fill(12.,x0); }
	}
      }
    } //-- chi2 --
#endif 
	
    //------  Chi^2 Track FIT -----

    bool isSingleTrack=false;
	
#ifdef USE_TRK
    //double chi2_max=20000;
    if (electron) {
      double chi2el = TrkFit(f125_el_fit, fx, "fx");
      f125_el_chi2->Fill(chi2el);
      //if ( 1 || chi2el>0. && chi2el<chi2_max) {
      if (chi2el>0. && chi2el<chi2_max) {
	f125_el_amp2d->Add(f125_el_amp2ds);
	Count("n_trk_el");
	N_trk_el++;
	isSingleTrack=true;
      }
    } else if (pion) {
      double chi2pi = TrkFit(f125_pi_fit, fx, "fx");
      f125_pi_chi2->Fill(chi2pi);
      //if ( 1 || chi2pi>0. && chi2pi<chi2_max) {
      if (chi2pi>0. && chi2pi<chi2_max) {
	f125_pi_amp2d->Add(f125_pi_amp2ds);
	Count("n_trk_pi");
	N_trk_pi++;
	isSingleTrack=true;
      }
    }
#endif

      //--- use GemTrk track (first or single ?? ----

    double gemtrk_x=-999.,  gemtrk_y=-999.,  gemtrk_E=0,  delta_x=1000.,  dx_thresh=5.;
    if (gem_scluster_count==1) {
      for (ULong64_t ic=0; ic<gem_scluster_count; ic++) {   // --- SRS cluster loop, actually only 0 ;
	double x=gem_scluster_x->at(ic); if (x<=0) gemtrk_x=x+50.; else gemtrk_x=x-50.;  gemtrk_x*=-1.;
	double y=gem_scluster_y->at(ic); if (y<=0) gemtrk_y=y+50.; else gemtrk_y=y-50.;  gemtrk_y*=-1.;
	gemtrk_E=gem_scluster_energy->at(ic);
      }
      delta_x=abs(ftrk.Eval(gemtrk_x)-x0);
    }
    
    //----- Single Track Event Histogram Filling -----
    if (isSingleTrack && gem_scluster_count==1 && delta_x<dx_thresh && 0.<gemtrk_x && gemtrk_x<20. && -20.<gemtrk_y && gemtrk_y <10.) { //--- single TRD track and single SRS hit --
      Count("1_TRK");

      for (ULong64_t i=0;i<f125_pulse_count; i++) {
		
	float peak_amp = f125_pulse_peak_amp->at(i);
	float ped = f125_pulse_pedestal->at(i);
	if (0 > ped || ped > 200 ) ped = 100;
	float amp=peak_amp-ped;
	float time=f125_pulse_peak_time->at(i);   int itime=f125_pulse_peak_time->at(i);
	int fADCSlot = f125_pulse_slot->at(i);
	int fADCChan = f125_pulse_channel->at(i);
			
	int gemChan = GetGEMChan(fADCChan, fADCSlot);
	int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
	int mmg2Chan = GetMMG2Chan(fADCChan, fADCSlot, RunNum);
	int rwellChan = GetRWELLChan(fADCChan, fADCSlot, RunNum);

	if (amp<0) amp=0;
	int MM_THR=50;
	if(electron) { //------------- electron -------
	  if (i==0) Count("1eTRK");
	  //printf("electron i=%llu gemChan=%d cal=%f ebeam=%f \n",i,gemChan, CalSum,Ebeam_el); sleep(1);
	  if (gemChan>-1) {
	    /*
#ifdef USE_TRK
	    f125_el_amp2ds->Fill(time,gemChan,amp);
	    f125_el_fit->Fill(time,gemChan,amp);
#else
	    f125_el_amp2d->Fill(time,gemChan,amp);
#endif
	    */
	    //------------ SRS - TRD clust correlation ---------------------
	    if (gem_scluster_count==1 && amp > 500 && time > 150. ) { //--- TR Radiator search
	      srs_etrd_corr->Fill(gemtrk_x,gemtrk_y,amp);
	    }
	    //--------------------------------------------------------------		

	    if (!(jentry%NPRT)) f125_el_evt->Fill(time,gemChan,amp);

	    f125_el->Fill(amp);
	    f125_el_clu2d->Fill(time,gemChan,1.);
	            	
	    gem_xpos.push_back(gemChan);
	    gem_dedx.push_back(amp);
	    gem_zpos.push_back(time);
	    gem_parID.push_back(electron);
	    gem_nhit++;
	    gem_zHist->Fill(time, amp);
	    //printf("Evt %d, Fill gem tree - Electron ... %d\n", event_num, gem_nhit);
	  }
	  if (amp>MM_THR && mmg1Chan>-1) {
	    mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
	    mmg1_f125_el->Fill(amp);
	    mmg1_f125_el_clu2d->Fill(time,mmg1Chan,1.);
	            	
	    mmg1_xpos.push_back(mmg1Chan);
	    mmg1_dedx.push_back(amp);
	    mmg1_zpos.push_back(time);
	    mmg1_parID.push_back(electron);
	    mmg1_nhit++;
	    mmg1_zHist->Fill(time, amp);
	  }
	  if (mmg2Chan>-1) {
	    mmg2_f125_el_amp2d->Fill(time,mmg2Chan,amp);
	    mmg2_f125_el->Fill(amp);
	    mmg2_f125_el_clu2d->Fill(time,mmg2Chan,1.);
		   		    
	    mmg2_xpos.push_back(mmg2Chan);
	    mmg2_dedx.push_back(amp);
	    mmg2_zpos.push_back(time);
	    mmg2_parID.push_back(electron);
	    mmg2_nhit++;
	    mmg2_zHist->Fill(time, amp);
	  }
	  if (rwellChan>-1) {
	    urw_f125_el_amp2d->Fill(time,rwellChan,amp);
	    urw_f125_el->Fill(amp);
	    urw_f125_el_clu2d->Fill(time,rwellChan,1.);
	            	
	    urw_xpos.push_back(rwellChan);
	    urw_dedx.push_back(amp);
	    urw_zpos.push_back(time);
	    urw_parID.push_back(electron);
	    urw_nhit++;
	    urw_zHist->Fill(time, amp);
	  }
	} else if (pion) {  //------ hadron/pion
	  if (i==0) Count("1pTRK");
	  if (gemChan>-1) {
	    /*
#ifdef USE_TRK
	    f125_pi_amp2ds->Fill(time,gemChan,amp);
	    //f125_pi_fit->Fill(time,gemChan,amp);
#else
	    f125_pi_amp2d->Fill(time,gemChan,amp);
#endif
	    */
	    if (!(jentry%NPRT)) f125_pi_evt->Fill(time,gemChan,amp);
	    f125_pi->Fill(amp);
	    f125_pi_clu2d->Fill(time,gemChan,1.);
	            	
	    gem_xpos.push_back(gemChan);
	    gem_dedx.push_back(amp);
	    gem_zpos.push_back(time);
	    gem_parID.push_back(electron);
	    gem_nhit++;
	    //printf("Evt %d, Fill gem tree - Pion ... %d\n", event_num, gem_nhit);
	  }
	  if (amp>MM_THR && mmg1Chan>-1) {
	    mmg1_f125_pi_amp2d->Fill(time,mmg1Chan,amp);
	    mmg1_f125_pi_clu2d->Fill(time,mmg1Chan,1.);
	    mmg1_f125_pi->Fill(amp);
	            	
	    mmg1_xpos.push_back(mmg1Chan);
	    mmg1_dedx.push_back(amp);
	    mmg1_zpos.push_back(time);
	    mmg1_parID.push_back(electron);
	    mmg1_nhit++;
	  }
	  if (mmg2Chan>-1) {
	    mmg2_f125_pi_amp2d->Fill(time,mmg2Chan,amp);
	    mmg2_f125_pi->Fill(amp);
	    mmg2_f125_pi_clu2d->Fill(time,mmg2Chan,1.);
		        	
	    mmg2_xpos.push_back(mmg2Chan);
	    mmg2_dedx.push_back(amp);
	    mmg2_zpos.push_back(time);
	    mmg2_parID.push_back(electron);
	    mmg2_nhit++;
	  }
	  if (rwellChan>-1) {
	    urw_f125_pi_amp2d->Fill(time,rwellChan,amp);
	    urw_f125_pi_clu2d->Fill(time,rwellChan,1.);
	    urw_f125_pi->Fill(amp);
	            	
	    urw_xpos.push_back(rwellChan);
	    urw_dedx.push_back(amp);
	    urw_zpos.push_back(time);
	    urw_parID.push_back(electron);
	    urw_nhit++;
	  }
	}
	 
	if (peak_amp-ped>f125_amp_max) {
	    f125_amp_max=peak_amp-ped;
	    gem_chan_max = fADCChan;
	}
			
	hCCor_ud->Fill(Ch_u,Ch_out);
	hCCCor_u->Fill(Ch_u,CalSum);
	hCCCor_dout->Fill(Ch_out,CalSum);
			
      } //--- end Fa125 Pulse Loop ---
		
      for (int i=1; i<21; i++) {
	gem_zHist_vect.push_back(gem_zHist->GetBinContent(i));
	mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
	mmg2_zHist_vect.push_back(mmg2_zHist->GetBinContent(i));
	urw_zHist_vect.push_back(urw_zHist->GetBinContent(i));
      }
		
    } //-- end single track condition --
	
    //======================= End Process Fa125 Pulse data ================================

    //==================================================================================================
    //                    Process Fa125  RAW data
    //==================================================================================================

#define USE_125_RAW
#ifdef USE_125_RAW
    if (jentry<MAX_PRINT) printf("------------------ Fadc125  wraw_count = %llu ---------\n", f125_wraw_count);
	
    for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
      if (jentry<MAX_PRINT) printf("F125:RAW: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
				   ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
				   ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
	
      int fadc_window = f125_wraw_samples_count->at(i);
		
      int fADCSlot = f125_wraw_slot->at(i);
      int fADCChan = f125_wraw_channel->at(i);
      int gemChan = GetGEMChan(fADCChan, fADCSlot);
      int amax=0;
      int tmax=0;
		
      for (int si=0; si<fadc_window; si++) {
	//printf("f125Loop:: %d fadc_window=%d\n",si,fadc_window);
	int time=si;
	int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si); // printf(" sample=%d adc=%d \n",si,adc);
	if (adc>amax) {
	  amax=adc;
	  tmax=si;
	}
	if (!(jentry%NPRT)) {
	  double adc_fill=adc;
	  if (electron)  {
	    f125_el_raw->Fill(time,gemChan,adc);
	  }else if (pion) {
	    f125_pi_raw->Fill(time,gemChan,adc);
	  }
	}
      } // --  end of samples loop
    } // -- end of fadc125 channels loop
  	  
    if (!(jentry%NPRT)) {
      int nx = f125_el_raw->GetNbinsX();
      int ny = f125_el_raw->GetNbinsY();
      double pedestal=100.;
      for (int ii=0; ii<nx; ii++) {
	for (int jj=0; jj<ny; jj++) {
	  if (electron) {
	    double cc = f125_el_raw->GetBinContent(ii, jj);
	    //printf("EL: %d %d cc=%f \n",ii,jj,cc);
	    if (cc == 0.) f125_el_raw->Fill(ii,jj,pedestal);
	  } else if (pion) {
	    double cc = f125_pi_raw->GetBinContent(ii, jj);
	    if (cc == 0.) f125_pi_raw->Fill(ii,jj,pedestal);
	  }
	}
      }
    }

    //=======================  End Fa125 RAW  process Loop  =====================================================
	
    if(electron){
      f125_el->Fill(f125_amp_max);
      //if((slot<6)||(slot==7&&chan<24))f125_el->Fill(f125_amp_max);
      //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_el->Fill(f125_amp_max);
    } else if (pion) {
      f125_pi->Fill(f125_amp_max);
      //if((slot<6)||(slot==7&&chan<24))f125_pi->Fill(f125_amp_max);
      //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_pi->Fill(f125_amp_max);
    }

#endif
	
    //=====================================================================================
    //===                     Event Display                                            ====
    //=====================================================================================
    
    if (jentry<MAX_PRINT || !(jentry%NPRT)) {
      c0->cd(1); f125_el_amp2d->Draw("colz");
      c0->cd(4); f125_pi_amp2d->Draw("colz");
      c0->cd(2); f125_el_evt->Draw("colz");         double chi2e = TrkFit(f125_el_evt,fx1,"fx1"); if (chi2e>0. && chi2e < chi2_max ) fx1.Draw("same"); 
      c0->cd(5); f125_pi_evt->Draw("colz");         double chi2p = TrkFit(f125_pi_evt,fx2,"fx2"); if (chi2p>0. && chi2p < chi2_max ) fx2.Draw("same"); 
      printf("========================>>>  Chi2 e=%f p=%f \n",chi2e,chi2p);
      //c0->cd(3); f125_el_chi2->Draw("colz");
      //c0->cd(6); f125_pi_chi2->Draw("colz");
      c0->cd(3); f125_el_raw->Draw("colz");  f125_el_evt->Draw("same"); 
      c0->cd(6); f125_pi_raw->Draw("colz");  f125_pi_evt->Draw("same"); 
      c0->cd(7); srs_trk_el->Draw("colz"); 
      c0->cd(8); srs_gem_x->Draw("colz");  ftrk.Draw("same"); 
      c0->cd(9); srs_etrd_corr->Draw("colz");
		
      c0->Modified();   c0->Update();
      if (NPRT<1000) sleep(1);
    }
    
	
    //--- Fill Track Hit Tree ---
    //printf("Filling hits tree ... ev=%d  Nhits=%d  Ntracks=%d \n", event_num, nhit, track_size);
    //printf("Filling hits tree ... ev=%d  Nhits=%d \n", event_num, nhit);
    if (gem_nhit>0) EVENT_VECT_GEM->Fill();
    if (mmg1_nhit>0)EVENT_VECT_MMG1->Fill();
    if (mmg2_nhit>0)EVENT_VECT_MMG2->Fill();
    if (urw_nhit>0)EVENT_VECT_URW->Fill();
    
  } // -- end of event loop
   
  cout<<" Total events= "<<jentry<< "  N_trk_el=" << N_trk_el++ << " N_trk_pi=" << N_trk_pi <<endl;
   
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================

  TFile* fOut;
  char rootFileName[256]; sprintf(rootFileName, "RootOutput/Run%06d_Output.root", RunNum);
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  HistList->Write("HistDQM", TObject::kSingleKey); 
  fOut->Close();
  delete fOut;
  
  //=====================================================================================
  //===                 S A V E   T R A C K   H I T   T T R E E                      ====
  //=====================================================================================

  if (save_hits_root) {
    printf(" Creating new TTree file... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_MMG1->Write();
    EVENT_VECT_MMG2->Write();
    EVENT_VECT_URW->Write();
    printf("OK, TTree File Write() ...  \n");
    fHits->Close();
    printf("OK, TTree File Close() ...  \n");
  }

  //=====================================================================================
  //===                 P L O T  H I S T O G R A M S                                  ===
  //=====================================================================================


  const char *OutputDir="RootOutput";
  int NN_MODE = 8;
  char ctit[120];
  //sprintf(G_DIR,"%s/SINGLERun_%06d",OutputDir,RunNum);
  sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
  sprintf(ctit,"File=%s",G_DIR);
  bool COMPACT=false;
  TCanvas *cc;
  int nxd=3;
  int nyd=5;

  //---  plot event display ---
  char pngname[120];  sprintf(pngname,"%s_evt.png",G_DIR);  c0->Print(pngname);
  char pdfname[120];  sprintf(pdfname,"%s_evt.pdf",G_DIR);  c0->Print(pdfname);

  //---------------------  page 1 --------------------
  htitle(" Calorimeter ");   // if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);                   hCal_occ->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[6]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_sum->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[3]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[4]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[5]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[0]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[1]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[2]->Draw();

  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_sum_el->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_sum_pi->Draw();

  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
 
 //---------------------  page 1a --------------------
  htitle(" Calorimeter Calib");    if (!COMPACT) cc=NextPlot(0,0);
  for (int i=0;i<NCAL; i++) {
    cc=NextPlot(nxd,nyd);  hCal_cal[i]->Draw();  fcal[i]->Draw("same"); 
  }
  //---------------------  page 1b --------------------
  htitle(" Calorimeter TRD_x correlation");    if (!COMPACT) cc=NextPlot(0,0);
  for (int i=0;i<NCAL; i++) {
    cc=NextPlot(nxd,nyd);  hCal_cor[i]->Draw("colz"); 
  }
  //---------------------  page 1c --------------------
  htitle(" Calorimeter TRK 2D correlation");    if (!COMPACT) cc=NextPlot(0,0);
  for (int i=0;i<NCAL; i++) {
    cc=NextPlot(nxd,nyd);  hCal_trk[i]->Draw("colz"); 
  }
   //---------------------  page 2 --------------------
  htitle(" Cherenkov (fa250)  ");   if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);  gPad->SetLogy();  hCher_u_adc->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy();  hCher_din_adc->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy();  hCher_dout_adc->Draw();
  cc=NextPlot(nxd,nyd);  hCher_u_time->Draw();
  cc=NextPlot(nxd,nyd);  hCher_din_time->Draw();
  cc=NextPlot(nxd,nyd);  hCher_dout_time->Draw();
  cc=NextPlot(nxd,nyd);  hCCor_ud->Draw("colz");
  cc=NextPlot(nxd,nyd);  hCCCor_u->Draw("colz");
  cc=NextPlot(nxd,nyd);  hCCCor_dout->Draw("colz");

  //---------------------  page 2a --------------------
  htitle(" SRS GEM-TRK  ");   if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);   gPad->SetLogy();  srs_ncl->Draw("");
  cc=NextPlot(nxd,nyd);  srs_trk_el->Draw("colz");
  cc=NextPlot(nxd,nyd);  srs_trk_pi->Draw("colz");
  cc=NextPlot(nxd,nyd);  srs_cal_corr->Draw("colz");
  cc=NextPlot(nxd,nyd);  srs_gem_x->Draw("colz");  ftrk.Draw("same"); 
  cc=NextPlot(nxd,nyd);  srs_gem_y->Draw("colz");
  cc=NextPlot(nxd,nyd);  srs_etrd_corr->Draw("colz");


 //---------------------  page 3 --------------------
  htitle("  GEMTRD (fa125) Amp ");    if (!COMPACT) cc=NextPlot(0,0);
  nxd=2; nyd=5;
  cc=NextPlot(nxd,nyd);   f125_el->Draw();
  cc=NextPlot(nxd,nyd);   f125_pi->Draw();
  cc=NextPlot(nxd,nyd);   mmg1_f125_el->Draw();
  cc=NextPlot(nxd,nyd);   mmg1_f125_pi->Draw();
  cc=NextPlot(nxd,nyd);   f125_el_amp2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   f125_pi_amp2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   mmg1_f125_pi_amp2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   urw_f125_el_amp2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   urw_f125_pi_amp2d->Draw("colz");

  //---------------------  page 4 --------------------
  htitle(" Clust   ");    if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);   f125_el_clu2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   f125_pi_clu2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   mmg1_f125_el_clu2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   mmg1_f125_pi_clu2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   urw_f125_el_clu2d->Draw("colz");
  cc=NextPlot(nxd,nyd);   urw_f125_pi_clu2d->Draw("colz");

  //---------------------  page 5 --------------------
  htitle(" Tracking   ");    if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);  f125_el_chi2->Draw("colz");
  cc=NextPlot(nxd,nyd);  f125_pi_chi2->Draw("colz");

  //--- close PDF file ----
  cc=NextPlot(-1,-1);
  //--- the end ---
  
}
//===============================================================
