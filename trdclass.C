#define trdclass_cxx
#include "trdclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSocket.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "PlotLib.C"
#include <TError.h>
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

#define NPRT 100000000
#define USE_TRK
//#define VERBOSE
#define SAVE_TRACK_HITS
#define MAX_PRINT 1
//#define SHOW_EVT_DISPLAY
//#define SAVE_PDF
//#define SHOW_EVTbyEVT
//#define USE_250_PULSE
#define USE_125_RAW
#define BUFSIZE 128000
#define MAX_CLUST 500
#define MAX_NODES 100
#define USE_GNN  1
#define USE_TCP  0
#define USE_FIT  0

//=================================================
//            TRD Prototype DAQ Mapping
//=================================================

// -- GEMTRD mapping --
int GetGEMChan(int ch, int slot) {
  int cardNumber = ch/24;
  int cardChannel = ch-cardNumber*24;
  int invCardChannel = 23-cardChannel;
  if (slot<6 || (slot==6 && ch<24)) {
    //--Remove noisy channels before returning the rest
    if (slot==3 && ch==8) {return -1;} else if (slot==3 && ch==10) {return -1;} else if (slot==6 && ch==16) {return -1;}
    else {return invCardChannel+cardNumber*24+(slot-3)*72.;}
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
      if (slot==9 && ch==40) {return -1;}
      else {return dchan - 240.;}
    }
  } else if (runNum>3261 && runNum<3279) { // -- Map #3
    if (slot==7&&ch>23) {
      if (slot==7 && ch==58) {return -1;}
      else {return dchan - 312.;}
    }
    if (slot==8 || (slot==9&&ch<48)) {
      if (slot==9 && ch==40) {return -1;}
      else {return dchan - 240.;}
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
    return -1;
  } else if (runNum>3147 && runNum<3262) { // -- Map #2
    if (slot==7 || (slot==6&&ch>23)) {
      if (slot==7 && ch==56) {return -1;} else {
      if (slot==6 && ch<48){return dchan - 144.;} else if (slot==6 && ch>47) {return dchan - 192.;} else if (slot==7 && ch<24) {return dchan - 240.;} else if (slot==7 && ch>47) {return dchan - 336.;} else {return dchan - 288.;} }
    }
  }
  return -1;
}
//==========  END TRD Prototype DAQ Mapping ===========

void trdclass::Loop() {
  
  if (fChain == 0) return;
  
//=================================================
//                TCP INITIALIZE
//=================================================
  
  #if  (USE_TCP==1)  //--- tcp send
  
  // Open connection to server
  TSocket *sock = new TSocket("localhost", 20250);
  // Wait till we get the start message
  char str[32];
  Float_t messlen  = 0;
  Float_t cmesslen = 0;
  static unsigned int BUFFER[BUFSIZE];
  float* FBUFFER = (float*) BUFFER;
  int LENEVENT;
  int nmod=1,hdr,itrg=0;
  unsigned int EVT_Type_BOR=(0x0&0x3)<<22;   //--  BOR=0x0
  unsigned int EVT_Type_EOR=(0x1&0x3)<<22;   //--  EOR=0x1
  unsigned int EVT_Type_DATA=(0x2&0x3)<<22;  //--  DATA=0x2
  static int HEADER[10];
  int modID=4;
#endif
//============== END TCP INITIALIZE ============

//==============================================================
//            B o o k    H i s t o g r a m s
//==============================================================

  TList *HistList = new TList();
  //--- Set TRD Prototype ADC thresholds
  int MM_THR=65;
  if (RunNum>3250) MM_THR=80;
  int URW_THR=105;
  int GEM_THR=120;

//============= Event Display (canvas 0) =============
  #ifdef SHOW_EVT_DISPLAY
    char c0Title[256];
    sprintf(c0Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c0 = new TCanvas("DISP",c0Title,1100,200,1500,1300);
    c0->Divide(4,3); c0->cd(1);
    TLine peak_line[100];
    //============== FPGA Display (canvas 1) ==========
    char c2Title[256]; sprintf(c2Title,"Event_Display_Run=%d",RunNum);
    TCanvas *c2 = new TCanvas("FPGA",c2Title,100,100,1000,1300);
    c2->Divide(1,3); c2->cd(1);
    
    //f125_el_evt_display = new TH2F("f125_el_evt_display","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  //HistList->Add(f125_el_evt_display);
    //f125_el_evt_display->SetStats(0);  f125_el_evt_display->SetMinimum(GEM_THR);  f125_el_evt_display->SetMaximum(1000.);
    //f125_pi_evt_display = new TH2F("f125_pi_evt_display","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  //HistList->Add(f125_pi_evt_display);
    //f125_pi_evt_display->SetStats(0); f125_pi_evt_display->SetMinimum(GEM_THR); f125_pi_evt_display->SetMaximum(1000.);
  #endif
  hcount= new TH1D("hcount","Count",3,0,3); HistList->Add(hcount);
  hcount->SetStats(0); hcount->SetFillColor(38);  hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
    hcount->SetCanExtend(TH1::kXaxis);
  #else
    hcount->SetBit(TH1::kCanRebin);
  #endif
  hNTracks = new TH1D("hNTracks","Number of Tracks in GEMTRD",12,-0.5,11.5);  HistList->Add(hNTracks);
  hNTracks_e = new TH1D("hNTracks_e","Number of Electron Tracks in GEMTRD",12,-0.5,11.5); HistList->Add(hNTracks_e);
  hNTracks_pi = new TH1D("hNTracks_pi","Number of Pion Tracks in GEMTRD",12,-0.5,11.5); HistList->Add(hNTracks_pi);
  
  f125_el_raw = new TH2F("f125_el_raw","GEM-TRD raw for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  //HistList->Add(f125_el_raw);
  f125_el_raw->SetStats(0);  f125_el_raw->SetMinimum(GEM_THR);   f125_el_raw->SetMaximum(1000.);
  f125_pi_raw = new TH2F("f125_pi_raw","GEM-TRD raw for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  //HistList->Add(f125_pi_raw);
  f125_pi_raw->SetStats(0); f125_pi_raw->SetMinimum(GEM_THR); f125_pi_raw->SetMaximum(1000.);
  
  //============ Calorimeter & Cherenkovs Plots ===============
  hCal_occ = new TH1F("hCal_occ"," Calorimeter Occupancy ; Cal. Cell Number ",8,-0.5,7.5);  HistList->Add(hCal_occ);
  hCal_sum = new TH1F("hCal_sum"," Calorimeter Sum (GeV)",100.,0.,25.); HistList->Add(hCal_sum);
  hCal_sum_el = new TH1F("hCal_sum_el"," Calorimeter Sum for electrons",100,0.,25.); HistList->Add(hCal_sum_el);
  hCal_sum_pi = new TH1F("hCal_sum_pi"," Calorimeter Sum for pions",100,0.,25.); HistList->Add(hCal_sum_pi);
  for (int cc=0; cc<NCAL; cc++) {
    char hName[128];  sprintf(hName,"hCal_adc%d",cc);
    char hTitle[128]; sprintf(hTitle,"Calorimeter ADC, cell%d",cc);
    hCal_adc[cc] = new TH1F(hName,hTitle,4096,-0.5,4095.5); HistList->Add(hCal_adc[cc]);
    //sprintf(hName,"hCal_cor%d",cc);  sprintf(hTitle,"Calorimeter ADC corr GemTrd X, cell%d",cc);
    //hCal_cor[cc] =  new TH2F(hName,hTitle,25,0.5,25.5,100,-55.,55.);  HistList->Add(hCal_cor[cc]);
    //sprintf(hName,"hCal_trk%d",cc);  sprintf(hTitle,"Calorimeter corr GemTRK, cell%d",cc);
    //hCal_trk[cc] =  new TH2F(hName,hTitle,100,-55.,55.,100,-55.,55.); HistList->Add(hCal_trk[cc]);
    sprintf(hName,"hCal_cal%d",cc);  sprintf(hTitle,"Calorimeter ADC Calib, cell%d",cc);
    hCal_cal[cc] = new TH2F(hName,hTitle,10,-0.5,4095.5,10,-5.,15.);  HistList->Add(hCal_cal[cc]);
  }
  h250_size = new TH1F("h250_size"," fa250 Raw data size",4096,0.5,4095.5); HistList->Add(h250_size);
  cal_el_evt = new TH2F("cal_el_evt"," Calorimeter Electron Event ; X ; Y ",3,-0.5,2.5,3,-0.5,2.5); HistList->Add(cal_el_evt);
  cal_el_evt->SetMinimum(-2.); cal_el_evt->SetMaximum(10.); cal_el_evt->SetStats(0);
  cal_pi_evt = new TH2F("cal_pi_evt"," Calorimeter Pion Event ; X ; Y ",3,-0.5,2.5,3,-0.5,2.5); HistList->Add(cal_pi_evt);
  cal_pi_evt->SetMinimum(-2.); cal_pi_evt->SetMaximum(10.); cal_pi_evt->SetStats(0);
/*
  //hCher_din_adc = new TH1F("hCher_din_adc"," Cherenkov Downstream (in) ADC ; ADC Amplitude ",4096,-0.5,4095.5); HistList->Add(hCher_din_adc);
  //hCher_din_time = new TH1F("hCher_din_time"," Cherenkov Downstream (in) Time ; Time Response (8ns)",300,-0.5,299.5); HistList->Add(hCher_din_time);
  //hCCCor_u = new TH2F("hCCCor_u"," Cherenkov Calorimeter Corr ; Upstream ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);  HistList->Add(hCCCor_u);
  //hCCCor_dout = new TH2F("hCCCor_dout"," Cherenkov Calorimeter Corr ; Downstream (out) ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);  HistList->Add(hCCCor_dout);
  */
  hCher_u_adc = new TH1F("hCher_u_adc"," Cherenkov Upstream ADC ; ADC Amplitude ",4096,-0.5,4095.5);  HistList->Add(hCher_u_adc);
  hCher_dout_adc = new TH1F("hCher_dout_adc"," Cherenkov Downstream (out) ADC ; ADC Amplitude ",4096,-0.5,4095.5);  HistList->Add(hCher_dout_adc);
  hCher_u_time = new TH1F("hCher_u_time"," Cherenkov Upstream Time ; Time Response (8ns)",300,-0.5,299.5);  HistList->Add(hCher_u_time);
  hCher_dout_time = new TH1F("hCher_dout_time"," Cherenkov Downstream (out) Time ; Time Response (8ns)",300,-0.5,299.5);  HistList->Add(hCher_dout_time);
  
  // ============= Detector Correlation Plots =============
  hCCor_ud = new TH2F("hCCor_ud"," Cherenkov Upstr./Downstr. Corr ; Upstream ; Downstream (out) ",400,-0.5,4095.5,400,0.5,4095.5);  HistList->Add(hCCor_ud);
  
  //-- GEM-TRKR & Prototype Correlations
  srs_gem_dx = new TH2F("srs_gem_dx","Correlation GEMTRD X & GEMTRKR X (Peaks) ; GEMTRD X [mm]; GEMTRKR Peak X [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(srs_gem_dx);
  srs_gem_dy = new TH2F("srs_gem_dy","Correlation GEMTRD X & GEMTRKR Y (Peaks) ; GEMTRD X [mm]; GEMTRKR Peak Y [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(srs_gem_dy);
  srs_mmg1_x = new TH2F("srs_mmg1_x","Correlation MMG1TRD Y & GEMTRKR X (Peaks); GEMTRKR Peak X [mm] ; MMG-1 Peak(SRS) Y [mm]",256,-64.,64.,256,-64.,64.);    HistList->Add(srs_mmg1_x);
  srs_mmg1_dx = new TH2F("srs_mmg1_dx","Correlation MMG1TRD X & GEMTRKR X (Peaks); MMG-1 X [mm]; GEMTRKR Peak X [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(srs_mmg1_dx);
  srs_mmg1_dy = new TH2F("srs_mmg1_dy","Correlation MMG1TRD X & GEMTRKR Y ; MMG-1 X [mm]; GEMTRKR Peak Y [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(srs_mmg1_dy);
  srs_mmg1_y = new TH2F("srs_mmg1_y","Correlation MMG1TRD Y & GEMTRKR Y (Peaks); GEMTRKR Peak Y [mm]; MMG-1 Peak(SRS) Y [mm]",256,-64.,64.,256,-64.,64.);    HistList->Add(srs_mmg1_y);
  if (RunNum>3147 && RunNum<3262) {
  srs_urw_x = new TH2F("srs_urw_x","Correlation uRWellTRD Y & GEMTRKR X (Peaks); GEMTRKR Peak X; uRWell Peak(SRS) Y [mm]",256,-64.,64.,128,-64.,64.);    HistList->Add(srs_urw_x);
  srs_urw_dx = new TH2F("srs_urw_dx","Correlation uRWellTRD X & GEMTRKR X (Peaks); uRWell X [mm]; GEMTRKR Peak X [mm]",120,-60.,60.,256,-64.,64.);    //HistList->Add(srs_urw_dx);
  srs_urw_dy = new TH2F("srs_urw_dy","Correlation uRWellTRD X & GEMTRKR Y (Peaks) ; uRWell X [mm] ; GEMTRKR Peak Y [mm]",120,-60.,60.,256,-64.,64.);    //HistList->Add(srs_urw_dy);
  srs_urw_y = new TH2F("srs_urw_y","Correlation uRWellTRD Y & GEMTRKR Y (Peaks); GEMTRKR Peak Y [mm]; uRWell Peak(SRS) Y [mm]",256,-64.,64.,128,-64.,64.);    HistList->Add(srs_urw_y);
  } else if (RunNum>3261) {
  srs_mmg2_x = new TH2F("srs_mmg2_x","Correlation MMG2TRD Y & GEMTRKR X (Peaks); GEMTRKR Peak X [mm]; MMG-2 Peak(SRS) Y [mm]",256,-64.,64.,128,-64.,64.);    HistList->Add(srs_mmg2_x);
  srs_mmg2_dx = new TH2F("srs_mmg2_dx","Correlation MMG2TRD X & GEMTRKR X (Peaks) ; MMG-2 X [mm]; GEMTRKR Peak X [mm]",72,-72.,72.,256,-64.,64.);    //HistList->Add(srs_mmg2_dx);
  srs_mmg2_dy = new TH2F("srs_mmg2_dy","Correlation MMG2TRD & GEMTRKR X&Y; MMG-2 X [mm] ; GEMTRKR Peak Y [mm]",72,-72.,72.,256,-64.,64.);    //HistList->Add(srs_mmg2_dy);
  srs_mmg2_y = new TH2F("srs_mmg2_y","Correlation MMG2TRD & GEMTRKR Y (Peaks); GEMTRKR Peak Y [mm]; MMG-2 Y [mm]",256,-64.,64.,128,-64.,64.);    HistList->Add(srs_mmg2_y);
  }
  //-- GEM-TRD & Other TRD Prototypes fADC Correlations
  if (RunNum>3147 && RunNum<3262) {
    gem_urw_x = new TH2F("gem_urw_x","Correlation GEMTRD X & uRWell X ; GEMTRD X [mm]; uRWell X [mm]",240,-60.,60.,120,-60.,60.);    //HistList->Add(gem_urw_x);
    gem_urw_y = new TH2F("gem_urw_y","Correlation GEMTRD X & uRWell Y ; GEMTRD X [mm] ; uRWell Y(SRS) [mm]",240,-60.,60.,128,-64.,64.);    //HistList->Add(gem_urw_y);
    urw_xy = new TH2F("urw_xy","Correlation uRWellTRD X&Y ; uRWell X [mm] ; uRWell Y(SRS) [mm]",120,-60.,60.,128,-64.,64.);    //HistList->Add(urw_xy);
    mmg1_urw_y = new TH2F("mmg1_urw_y","Correlation MMG1 Y & uRWell Y ; uRWell Y(SRS) [mm]; MMG-1 Y(SRS) [mm]",128,-64.,64.,256,-64.,64.);    //HistList->Add(mmg1_urw_y);
  } else if (RunNum>3261) {
    gem_mmg2_x = new TH2F("gem_mmg2_x","Correlation GEMTRD X & MMG2 X ; GEMTRD X [mm]; MMG-2 X [mm]",240,-60.,60.,72,-72.,72.);    //HistList->Add(gem_mmg2_x);
    gem_mmg2_y = new TH2F("gem_mmg2_y","Correlation GEMTRD X & MMG2 Y ; GEMTRD X [mm] ; MMG-2 Y(SRS) [mm]",240,-60.,60.,128,-64.,64.);    //HistList->Add(gem_mmg2_y);
    mmg2_xy = new TH2F("mmg2_xy","Correlation MMG2TRD X&Y ; MMG-2 X [mm] ; MMG-2 Y(SRS) [mm]",72,-72.,72.,128,-64.,64.);    //HistList->Add(mmg2_xy);
  }
  gem_mmg1_x = new TH2F("gem_mmg1_x","Correlation GEMTRD X & MMG1 X ; GEMTRD X [mm]; MMG-1 X [mm]",240,-60.,60.,240,-60.,60.);    //HistList->Add(gem_mmg1_x);
  gem_mmg1_y = new TH2F("gem_mmg1_y","Correlation GEMTRD X & MMG1 Y ; GEMTRD X [mm] ; MMG-1 Y(SRS) [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(gem_mmg1_y);
  mmg1_xy = new TH2F("mmg1_xy","Correlation MMG1TRD X&Y ; MMG-1 X [mm] ; MMG-1 Y(SRS) [mm]",240,-60.,60.,256,-64.,64.);    //HistList->Add(mmg1_xy);
  
  //-- GEM-TRKR & PID/Beam Correlations
  //srs_cal_corr = new TH2F("srs_cal_corr","Correlation GEMTRKR & CAL; X ; Y ",100,-55.,55.,100,-55.,55.);  HistList->Add(srs_cal_corr);
  //srs_gemtrd_el = new TH2F("srs_gemtrd_el","Correlation GEMTRKR & GEMTRD Electron Hit; GEMTRKR X [mm] ; GEMTRKR Y (Peaks) [mm] ",110,-55.,55.,110,-55.,55.);  HistList->Add(srs_gemtrd_el);
  //srs_etrd_beam = new TH2F("srs_etrd_beam","Correlation GEMTRKR & Beam; GEMTRKR Peak X [mm] ; GEMTRKR Peak Y [mm] ",100,-55.,55.,100,-55.,55.); HistList->Add(srs_etrd_beam);
  //srs_gemtrd_pion = new TH2F("srs_gemtrd_pion","Correlation GEMTRKR & GEMTRD Pion Hit; GEMTRKR X [mm] ; GEMTRKR Y [mm] ",110,-55.,55.,110,-55.,55.);  HistList->Add(srs_gemtrd_pion);
  //srs_etrd_ratio = new TH2F("srs_etrd_ratio","Correlation TRK ratio; X ; Y ",100,-55.,55.,100,-55.,55.);         HistList->Add(srs_etrd_ratio);
  
  //=============== Track Fitting & chi^2 ==================
  gErrorIgnoreLevel = kBreak; // Suppress warning messages from empty chi^2 fit data
  TF1 fx("fx","pol1",100,190); //-- Linear function fit (ax+b) over time response window
  TF1 fx1("fx1","pol1",100,190);
  TF1 fx2("fx2","pol1",100,190);
  TF1 fx_mmg1("fx_mmg1","pol1",80,190);
  TF1 fx_urw("fx_urw","pol1",80,190);
  TF1 fx_mmg2("fx_mmg2","pol1",80,190);
  
  //-- GEMTRD & GEMTRKR alignment
  double x_1=-48., y_1=-55., x_2=54., y_2=51.;
  double a_slope=(y_2-y_1)/(x_2-x_1);
  double b_intercept=y_1-a_slope*x_1;
  double x_cut1=-50., x_cut2=+50., y_cut1=-50., y_cut2=+50.;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,a_slope);
  ftrk.SetParameter(1,b_intercept);
  double gemtrkr_x2ch=-999.;
  
  //======== GEM-TRKR ========
  singleTrackIndex = new TH2F("singleTrackIndex","GEMTracker Correlated Hits with 1 Track; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);      HistList->Add(singleTrackIndex);
  singleTrackIndex_e = new TH2F("singleTrackIndex_e","GEMTracker Correlated Electrons with 1 Track; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);    HistList->Add(singleTrackIndex_e);
  singleTrackIndex_pi = new TH2F("singleTrackIndex_pi","GEMTracker Pions Hits with 1 Track; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);    HistList->Add(singleTrackIndex_pi);
  multiTrackIndex = new TH2F("multiTrackIndex","GEMTracker Correlated Hits with 2 Tracks; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);    HistList->Add(multiTrackIndex);
  multiTrackIndex_e = new TH2F("multiTrackIndex_e","GEMTracker Correlated Electrons with 2 Tracks; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);    HistList->Add(multiTrackIndex_e);
  multiTrackIndex_pi = new TH2F("multiTrackIndex_pi","GEMTracker Correlated Pions with 2 Tracks; GEMTRKR X Num Hits; GEMTRKR Y Num Hits",10,-0.5,9.5,10,-0.5,9.5);    HistList->Add(multiTrackIndex_pi);
  
  hgemtrkr_peak_xy = new TH2F("hgemtrkr_peak_xy","GEM-TRKR Peak X-Y Correlation (mm); Peak X [mm]; Peak Y [mm] ",256,-64.,64.,256,-64.,64.);    HistList->Add(hgemtrkr_peak_xy);
  hgemtrkr_ch_xy = new TH2F("hgemtrkr_ch_xy","GEM-TRKR Peak X-Y Correlation (CHANNELS); Peak X Channel Number; Peak Y Channel Number",256,-0.5,255.5,256,-0.5,255.5);    HistList->Add(hgemtrkr_ch_xy);
  hgemtrkr_peak_x = new TH1F("hgemtrkr_peak_x"," GEM-TRKR Peak X ; X [mm] ",256,-64.,64.);  HistList->Add(hgemtrkr_peak_x);
  hgemtrkr_peak_y = new TH1F("hgemtrkr_peak_y"," GEM-TRKR Peak Y ; Y [mm] ",256,-64.,64.);  HistList->Add(hgemtrkr_peak_y);
  
  //============= Prototype ADC Amplitude Distributions ============
  f125_el = new TH1F("f125_el","GEM-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);  HistList->Add(f125_el);
  f125_pi = new TH1F("f125_pi","GEM-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);  HistList->Add(f125_pi);
  f125_el_max = new TH1F("f125_el_max","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(f125_el_max);
  f125_pi_max = new TH1F("f125_pi_max","GEM-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(f125_pi_max);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg1_f125_el);
  mmg1_f125_pi = new TH1F("mmg1_f125_pi","MMG1-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg1_f125_pi);
  mmg1_f125_el_max = new TH1F("mmg1_f125_el_max","MMG1-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg1_f125_el_max);
  mmg1_f125_pi_max = new TH1F("mmg1_f125_pi_max","MMG1-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg1_f125_pi_max);
  if (RunNum>3147 && RunNum<3262) {
    urw_f125_el = new TH1F("urw_f125_el","uRWell-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(urw_f125_el);
    urw_f125_pi = new TH1F("urw_f125_pi","uRWell-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(urw_f125_pi);
    urw_f125_el_max = new TH1F("urw_f125_el_max","uRWell-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(urw_f125_el_max);
    urw_f125_pi_max = new TH1F("urw_f125_pi_max","uRWell-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(urw_f125_pi_max);
  } else if(RunNum>3261) {
    mmg2_f125_el = new TH1F("mmg2_f125_el","MMG2-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg2_f125_el);
    mmg2_f125_pi = new TH1F("mmg2_f125_pi","MMG2-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg2_f125_pi);
    mmg2_f125_el_max = new TH1F("mmg2_f125_el_max","MMG2-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg2_f125_el_max);
    mmg2_f125_pi_max = new TH1F("mmg2_f125_pi_max","MMG2-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096); HistList->Add(mmg2_f125_pi_max);
  }
  
  //-- Prototype ADC Amplitudes in Time (2D) for Single-Track Evts Only
  f125_el_amp2ds = new TH2F("f125_el_amp2ds","GEM-TRD Single-Trk Amp. in Time (Electrons); Time Response (8ns) ; GEMTRD Channel (X)",250,0.5,250.5,240,0.5,240.5);  HistList->Add(f125_el_amp2ds);
  f125_pi_amp2ds = new TH2F("f125_pi_amp2ds","GEM-TRD Single-Trk Amp. in Time (Pions); Time Response (8ns) ; GEMTRD Channel (X)",250,0.5,250.5,240,0.5,240.5);      HistList->Add(f125_pi_amp2ds);
  mmg1_f125_el_amp2ds = new TH2F("mmg1_f125_el_amp2ds","MMG1-TRD Single-Trk Amp. in Time (Electrons); Time Response (8ns) ; MMG1 Channel (X)",250,0.5,250.5,240,0.5,240.5);  HistList->Add(mmg1_f125_el_amp2ds);
  mmg1_f125_pi_amp2ds = new TH2F("mmg1_f125_pi_amp2ds","MMG1-TRD Single-Trk Amp. in Time (Pions); Time Response (8ns) ; MMG1 Channel (X)",250,0.5,250.5,240,0.5,240.5);      HistList->Add(mmg1_f125_pi_amp2ds);
  if (RunNum>3147 && RunNum<3262) {
    urw_f125_el_amp2ds = new TH2F("urw_f125_el_amp2ds","uRWell-TRD Single-Trk Amp. in Time (Electrons); Time Response (8ns) ; uRWell Channel (X)",250,0.5,250.5,120,0.5,120.5);  HistList->Add(urw_f125_el_amp2ds);
    urw_f125_pi_amp2ds = new TH2F("urw_f125_pi_amp2ds","uRWell-TRD Single-Trk Amp. in Time (Pions); Time Response (8ns) ; uRWell Channel (X)",250,0.5,250.5,120,0.5,120.5);      HistList->Add(urw_f125_pi_amp2ds);
  } else if (RunNum>3261) {
    mmg2_f125_el_amp2ds = new TH2F("mmg2_f125_el_amp2ds","MMG2-TRD Single-Trk Amp. in Time (Electrons); Time Response (8ns) ; MMG2 Channel (X)",250,0.5,250.5,64,48.5,112.5);  HistList->Add(mmg2_f125_el_amp2ds);
    mmg2_f125_pi_amp2ds = new TH2F("mmg2_f125_pi_amp2ds","MMG2-TRD Single-Trk Amp. in Time (Pions); Time Response (8ns) ; MMG2 Channel (X)",250,0.5,250.5,64,48.5,112.5);      HistList->Add(mmg2_f125_pi_amp2ds);
  }
  
  //-- Prototype Single-Track Amplitudes in Time (for NN)
  gem_zHist = new  TH1F("gem_zHist", "gem_zHist", 20, 80., 200.);
  mmg1_zHist = new  TH1F("mmg1_zHist", "mmg1_zHist", 20, 80., 200.);
  mmg2_zHist = new  TH1F("mmg2_zHist", "mmg2_zHist", 20, 80., 200.);
  urw_zHist = new  TH1F("urw_zHist", "urw_zHist", 20, 80., 200.);
  
  //---  Clustering Histograms (Pattern Recognition Track-Finding) ---
  int nx0=100;    int ny0=250;
  hevt  = new TH2F("hevt"," Event display; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.); hevt->SetStats( 0 ); hevt->SetMaximum(10.); //HistList->Add(hevt);
  hevtc = new TH2F("hevtc"," Clustering ; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);  //HistList->Add(hevtc);
  hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
  hevti = new TH2F("hevti"," ML-FPGA response; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.);  hevti->SetStats( 0 ); hevti->SetMaximum(10.);  //HistList->Add(hevti);
  hevtf = new TH2F("hevtf"," Clusters for FPGA ; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.);  hevtf->SetStats( 0 ); hevtf->SetMaximum(10.);  //HistList->Add(hevtf);
  
  //============================ End Histogram Booking =====================================
  
  //--- Calorimeter Calibration
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
  
  // ----- Declare Beam Energy -----
  double Ebeam=10.; // GeV
  if  (3131 <= RunNum && RunNum <= 3171) {    Ebeam=10.;   }   // NEG 10 GeV
  if  (3172 <= RunNum && RunNum <= 3176) {    Ebeam=120.;  }   // protons
  if  (3177 <= RunNum && RunNum <= 3210) {    Ebeam=10.;   }   //  NEG 10 GeV
  if  (3211 <= RunNum && RunNum <= 3219) {    Ebeam=3.;    }   // NEG 3 GeV
  if  (3220 <= RunNum && RunNum <= 3226) {    Ebeam=120.;  }   // protons
  if  (3227 <= RunNum && RunNum <= 3227) {    Ebeam=10.;   }   // POS 10 GeV
  if  (3228 <= RunNum && RunNum <= 3240) {    Ebeam=120.;  }   // protons
  if  (3241 <= RunNum && RunNum <= 3253) {    Ebeam=3.;    }   // NEG 3GeV
  if  (3255 <= RunNum && RunNum <= 3261) {    Ebeam=120.;  }   // protons
  if  (3272 <= RunNum && RunNum <= 3288) {    Ebeam=10.;   }   // NEG 10 GeV
  if  (3289 <= RunNum && RunNum <= 3290) {    Ebeam=120.;  }   // protons
  double Ebeam_el=0.2*Ebeam;
  double Ebeam_pi=0.1*Ebeam;
  
  //=======================================================
  //     Create TTrees of Single Track Evt Info for NN
  //=======================================================
  
  #ifdef SAVE_TRACK_HITS
    TFile* fHits;
    char hitsFileName[256]; sprintf(hitsFileName, "RootOutput/trd_singleTrackHits_Run_%06d.root", RunNum);
    fHits = new TFile(hitsFileName, "RECREATE");
    //-- GEM-TRD
    EVENT_VECT_GEM = new TTree("gem_hits","GEM TTree with single track hit info");
    EVENT_VECT_GEM->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_GEM->Branch("nhit",&gem_nhit,"gem_nhit/I");
    EVENT_VECT_GEM->Branch("xpos",&gem_xpos);
    EVENT_VECT_GEM->Branch("zpos",&gem_zpos);
    EVENT_VECT_GEM->Branch("dedx",&gem_dedx);
    EVENT_VECT_GEM->Branch("parID",&gem_parID);
    EVENT_VECT_GEM->Branch("zHist",&gem_zHist_vect);
    //-- MMG1-TRD
    EVENT_VECT_MMG1 = new TTree("mmg1_hits","MMG1 TTree with single track hit info");
    EVENT_VECT_MMG1->Branch("event_num",&event_num,"event_num/I");
    EVENT_VECT_MMG1->Branch("nhit",&mmg1_nhit,"mmg1_nhit/I");
    EVENT_VECT_MMG1->Branch("xpos",&mmg1_xpos);
    EVENT_VECT_MMG1->Branch("zpos",&mmg1_zpos);
    EVENT_VECT_MMG1->Branch("dedx",&mmg1_dedx);
    EVENT_VECT_MMG1->Branch("parID",&mmg1_parID);
    EVENT_VECT_MMG1->Branch("zHist",&mmg1_zHist_vect);
    if (RunNum>3261) {
    //-- MMG2-TRD
      EVENT_VECT_MMG2 = new TTree("mmg2_hits","MMG2 TTree with single track hit info");
      EVENT_VECT_MMG2->Branch("event_num",&event_num,"event_num/I");
      EVENT_VECT_MMG2->Branch("nhit",&mmg2_nhit,"mmg2_nhit/I");
      EVENT_VECT_MMG2->Branch("xpos",&mmg2_xpos);
      EVENT_VECT_MMG2->Branch("zpos",&mmg2_zpos);
      EVENT_VECT_MMG2->Branch("dedx",&mmg2_dedx);
      EVENT_VECT_MMG2->Branch("parID",&mmg2_parID);
      EVENT_VECT_MMG2->Branch("zHist",&mmg2_zHist_vect);
    } else if (RunNum>3147 && RunNum<3262) {
    //-- uRWell-TRD
      EVENT_VECT_URW = new TTree("urw_hits","uRWELL TTree with single track hit info");
      EVENT_VECT_URW->Branch("event_num",&event_num,"event_num/I");
      EVENT_VECT_URW->Branch("nhit",&urw_nhit,"urw_nhit/I");
      EVENT_VECT_URW->Branch("xpos",&urw_xpos);
      EVENT_VECT_URW->Branch("zpos",&urw_zpos);
      EVENT_VECT_URW->Branch("dedx",&urw_dedx);
      EVENT_VECT_URW->Branch("parID",&urw_parID);
      EVENT_VECT_URW->Branch("zHist",&urw_zHist_vect);
    }
  #endif
  
//===================================================================================
//                      E v e n t    L o o p
//===================================================================================
  
  int e_1TRK=0;
  int pi_1TRK=0;
  int _1TRK=0;
  int pi_CC=0;
  int el_CC=0;
  int e_CHR=0;
  int e_CHR_Up=0;
  int n_e_1=0, n_pi_1=0, n_e_2=0, n_pi_2=0;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  Long64_t jentry=0;
  printf("=============== Begin Event Loop: 1st Evt=%lld, MaxEvt=%lld ================ \n",FirstEvt,MaxEvt);
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) {
    
    event_num = jentry;
    Count("EVT");
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if (jentry<MAX_PRINT || !(jentry%1000)) {
      printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_raw_count=%llu, srs_prerecon_count=%llu, gem_peak_count=%llu \n", jentry, f125_wraw_count, f125_pulse_count, f250_wraw_count, srs_raw_count, srs_prerecon_count, gem_peak_count);
    }
    
//==================================================================
//                  E v e n t    D i s p l a y
//==================================================================
    #ifdef SHOW_EVT_DISPLAY
      if (jentry<MAX_PRINT) {
        printf("-------------------- Pulse  125  ---------------------------\n");
        for (ULong64_t i=0;i<f125_pulse_count; i++) {
          printf("F125:pulse: i=%lld  sl=%d, ch=%d, npk=%d time=%d amp=%d ped=%d \n"
           ,i,f125_pulse_slot->at(i),f125_pulse_channel->at(i),f125_pulse_npk->at(i)
           ,f125_pulse_peak_time->at(i),f125_pulse_peak_amp->at(i),f125_pulse_pedestal->at(i));
        }
        printf("-------------------- Pulse  250  n=%lld ---------------------------\n",f250_pulse_count);
        for (ULong64_t i=0;i<f250_pulse_count; i++) {
          printf("F250:pulse: i=%lld  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n"
           ,i,f250_pulse_slot->at(i),f250_pulse_channel->at(i),f250_pulse_pulse_number->at(i)
           ,f250_pulse_course_time->at(i),f250_pulse_pulse_peak->at(i),f250_pulse_pedestal->at(i)/4.);
        }
        printf("-------------------- Raw  125  ---------------------------\n");
        for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
          printf("F125:raw: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
           ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
           ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
        }
        printf("-------------------- SRS   ev=%lld   ---------------------------\n",jentry);
        printf("SRS:gem_scluster:  cnt=%lld \n",gem_scluster_count);
        for (ULong64_t i=0;i<gem_scluster_count; i++) { // --- SRS cluster loop
          printf("SRS:clusters:  i=%lld  X=%f Y=%f  \n",i,gem_scluster_x->at(i), gem_scluster_y->at(i));
        }
        printf("SRS:srs_prerecon:  cnt=%lld cnt_x=%ld cnt_y=%ld \n",srs_prerecon_count, srs_prerecon_x->size(), srs_prerecon_y->size() );
        for (ULong64_t i=0;i<srs_prerecon_count; i++) {
          printf("SRS:srs_prerecon: i=%lld %f %f \n",i,srs_prerecon_x->at(i),srs_prerecon_y->at(i));
        }
        printf("SRS:srs_peak:  cnt=%lld \n",gem_peak_count);
        for (ULong64_t i=0;i<gem_peak_count; i++) {
    	    printf("SRS:srs_peak: i=%lld id=%d name=%s idx=%d apv=%d Amp=%f wid=%f E=%f Pos=%f \n"
           ,i,gem_peak_plane_id->at(i),gem_peak_plane_name->at(i).c_str(),gem_peak_index->at(i), gem_peak_apv_id->at(i), gem_peak_height->at(i)
           ,gem_peak_width->at(i), gem_peak_area->at(i), gem_peak_real_pos->at(i));
        }
      }
    #endif
    
//====================================================================
//       Process Fa250 Pulse data (Calorimeter & Cherenkovs)
//====================================================================
    
    #ifdef USE_250_PULSE
      printf("-------------------- Pulse  250  n=%lld ---------------------------\n",f250_pulse_count);
      for (int i=0;i<f250_pulse_count; i++) {
        printf("F250:: i=%d  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n", i, f250_pulse_slot->at(i), f250_pulse_channel->at(i), f250_pulse_pulse_number->at(i), f250_pulse_course_time->at(i), f250_pulse_pulse_peak->at(i), f250_pulse_pedestal->at(i)/4.);
      }
      if (jentry<MAX_PRINT) printf("------------------ Fadc250  wraw_count = %llu ---------\n", f250_wraw_count);
      h250_size->Fill(f250_wraw_count);
    #endif
    
    double CalSum=0., _calsum=0., _calsum_ecut=0., _calsum_picut=0.;
    double Ch_u=0.;
    double Ch_out=0.;
    bool electron_chUp=false;
    bool electron_ch=false;
    bool electron=false;
    bool pion=false;
    double Ecal[NCAL]; for (int i=0; i<NCAL; i++) Ecal[i]=0;
    
    for (ULong64_t i=0; i<f250_wraw_count; i++) {
      int fadc_chan = f250_wraw_channel->at(i);
      int fadc_window = f250_wraw_samples_count->at(i);
      hCal_occ->Fill(fadc_chan+0.);
      int amax=0;
      int tmax=9;
      for (int si=0; si<fadc_window; si++) {
        int adc = f250_wraw_samples->at(f250_wraw_samples_index->at(i)+si);
        if (adc>amax) {
          amax=adc;
          tmax=si;
        }
      } //--  end of fADC 250 samples loop
      if (fadc_chan<NCAL) { //-- Cal Energy Sum
        Ecal[fadc_chan]=fcal[fadc_chan]->Eval(amax);
        hCal_adc[fadc_chan]->Fill(amax);
        CalSum+=Ecal[fadc_chan];
      } else { // Cherenkov
        if (fadc_chan==13) { if(amax>130)electron_chUp=true; hCher_u_adc->Fill(amax);  hCher_u_time->Fill(tmax); Ch_u=amax; Count("eCHR_Up"); e_CHR_Up++;}
        if (fadc_chan==15) { if(amax>300)electron_ch=true; hCher_dout_adc->Fill(amax);  hCher_dout_time->Fill(tmax); Ch_out=amax; Count("eCHR"); e_CHR++;}
      }
    } //--- End of Fadc250 Channel Loop
    
//=======================================================
//                   S e t    P I D
//=======================================================
    if (CalSum>0.) {Count("calSum"); _calsum++;}
    if (CalSum>Ebeam_el) {Count("calSumEl"); _calsum_ecut++;}
    if (CalSum<Ebeam_pi && CalSum>0.) {Count("calSumPi"); _calsum_picut++;}
    
    //--- CHERNEKOV ONLY PID
    //if (electron_ch  && electron_chUp) { electron=true;  Count("elCC"); el_CC++;}
    //if (!electron_ch && !electron_chUp) { pion=true;  Count("piCC"); pi_CC++;}

    //--- CHERENKOV & CALORIMETER PID
    if (electron_ch  && CalSum>Ebeam_el) { electron=true;  Count("elCC"); el_CC++;}
    if (!electron_ch && CalSum<Ebeam_pi) { pion=true;  Count("piCC"); pi_CC++;}
    
    hCal_sum->Fill(CalSum);
    if (electron) {
      hCal_sum_el->Fill(CalSum);
    } else if (pion) {
      hCal_sum_pi->Fill(CalSum);
    }

    
//============ END Process Fa250 Pulse data (Calorimeter & Cherenkovs) =============
    
//==========================================================================================
//                   Process  Fa125 Pulse Data (TRDs)
//==========================================================================================
    
    #ifdef SHOW_EVT_DISPLAY
      if (!(jentry%NPRT)) {
        if(electron) {
          //f125_el_evt_display->Reset();
          f125_el_raw->Reset();
          cal_el_evt->Reset();
          for (int cc=0; cc<NCAL; cc++) {
            int ix=cc%3; int iy=cc/3;
            if (cc<6) cal_el_evt->Fill(ix,iy,Ecal[cc]); else cal_el_evt->Fill(1,2,Ecal[cc]);
          }
        } else if (pion) {
          //f125_pi_evt_display->Reset();
          f125_pi_raw->Reset();
          cal_pi_evt->Reset();
          for (int cc=0; cc<NCAL; cc++) {
            if (cc<6) cal_pi_evt->Fill(cc%3,cc/3,Ecal[cc]); else cal_pi_evt->Fill(1,2,Ecal[cc]);
          }
        }
      }
    #endif
    
//=====================================================
//        Process Fa125 RAW data
//=====================================================
    
    #ifdef USE_125_RAW
      
//==============================================================
//           GEMTRD Pattern Recognition Tracking
//==============================================================
      
      hevt->Reset();
      hevtc->Reset();
      hevti->Reset();
      hevtf->Reset();
      
      for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
        
        int fadc_window = f125_wraw_samples_count->at(i);
        int fADCSlot = f125_wraw_slot->at(i);
        int fADCChan = f125_wraw_channel->at(i);
        int gemChan = GetGEMChan(fADCChan, fADCSlot);
        int amax=0;
        int tmax=0;
        if (gemChan<0) continue; ///TODO:IMPLEMENT IN OTHER LOOPS
        //DEDX_THR = 100;
        int TimeWindowStart = 95;
        
        for (int si=0; si<fadc_window; si++) {
          int time=si;
          int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si);
          if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
          if (adc>amax) {
            amax=adc;
            tmax=si;
          }
          if (adc>GEM_THR) {
            double adc_fill=adc;
            if (electron) {
              f125_el_raw->Fill(time,gemChan,adc);
            } else if (pion) {
              f125_pi_raw->Fill(time,gemChan,adc);
            }
            time-=TimeWindowStart;
        	  if ( 0 > time || time > 100 ) continue; // --- drop early and late hits ---
  	        hevtc->SetBinContent(100-time,gemChan,adc/100.);
  	        double xp = gemChan/250.*100-50.;
  	        double zp = (time)/(100.)*30;
  	        double ap = adc/100.;
  	        hevt->SetBinContent(100-time,gemChan,adc/100.);
          }
        } // --  end of samples loop
      } // -- end of fadc125 channels loop
      
      #if (USE_GNN==0)
        int nx = f125_el_raw->GetNbinsX();
        int ny = f125_el_raw->GetNbinsY();
        double pedestal=100.;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("+++                           RAW TRD DATA                                                         +++\n");
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        #endif
        for (int ii=0; ii<nx; ii++) {
          for (int jj=0; jj<ny; jj++) {
            if (electron) {
              double cc = f125_el_raw->GetBinContent(ii, jj);
            } else if (pion) {
              double cc = f125_pi_raw->GetBinContent(ii, jj);
            }
          }
        }
      #endif //================== End if (USE_GNN==0) ===========================
      
      #if (USE_GNN>0)
        // ------------------------   hist dist clustering         ------------------------
        float clust_Xpos[MAX_CLUST];
        float clust_Ypos[MAX_CLUST];
        float clust_Zpos[MAX_CLUST];
        float clust_dEdx[MAX_CLUST];
        float clust_Size[MAX_CLUST];
        float clust_Width[MAX_CLUST][3];  // y1, y2, dy ; strips
        float clust_Length[MAX_CLUST][3]; // x1, x2, dx ; time
        
        float hits_Xpos[500];
        float hits_Ypos[500];
        float hits_Zpos[500];
        float hits_dEdx[500];
        
        for (int k=0; k<MAX_CLUST; k++) {
          clust_Xpos[k]=0; clust_Ypos[k]=0; clust_Zpos[k]=0; clust_dEdx[k]=0;  clust_Size[k]=0;
          clust_Width[k][0]=999999;   	clust_Width[k][1]=-999999;   	clust_Width[k][2]=0;
          clust_Length[k][0]=999999;  	clust_Length[k][1]=-999999;  	clust_Length[k][2]=0;
        }
        float CL_DIST=2.7; // mm
        int nclust=0;
        
        TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
        TH2F* hpc = hevtc;
        
        int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY();
        double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx);
        double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny);
        double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
        #ifdef VERBOSE
          printf("nx=%d,ny=%d,xmi=%f,xma=%f,ymi=%f,yma=%f\n",nx,ny,xmi,xma,ymi,yma);
        #endif
        double THR2 = 1.2;
        for (int ix=0; ix<nx; ix++) {  //-------------------- clustering loop ------------------------------------
          for (int iy=0; iy<ny; iy++) {
            double c1 = hpc->GetBinContent(ix,iy);   // hpc->SetBinContent(ix,iy,5.);         // energy
            double x1=double(ix)/double(nx)*(xma-xmi)+xmi-binx/2.;    // drift time
            double y1=double(iy)/double(ny)*(yma-ymi)+ymi-biny/2.;    // X strip
            if (c1<THR2) continue;
            if (nclust==0) {
              clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++; continue;
            }
            int added=0;
            for (int k=0; k<nclust; k++) {
              double dist=sqrt(pow((y1-clust_Xpos[k]),2.)+pow((x1-clust_Zpos[k]),2.)); //--- dist hit to clusters
              if (dist<CL_DIST) {
                clust_Xpos[k]=(y1*c1+clust_Xpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new X pos
                clust_Zpos[k]=(x1*c1+clust_Zpos[k]*clust_dEdx[k])/(c1+clust_dEdx[k]);  //--  new Z pos
                clust_dEdx[k]=c1+clust_dEdx[k];  // new dEdx
                clust_Size[k]=1+clust_Size[k];  // clust size in pixels
                if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
                if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
                hpc->SetBinContent(ix,iy,k+1.);
                added=1; break;
              }
            }
            if (added==0) {
              if (nclust+1>=MAX_CLUST) continue;
              clust_Xpos[nclust]=y1; clust_Ypos[nclust]=0; clust_Zpos[nclust]=x1;  clust_dEdx[nclust]=c1;  clust_Size[nclust]=1;
              clust_Width[nclust][0]=y1;   	clust_Width[nclust][1]=y1;   	clust_Width[nclust][2]=0;
              clust_Length[nclust][0]=x1;  	clust_Length[nclust][1]=x1;  	clust_Length[nclust][2]=0;
              nclust++;
            }
          }
        } //----------------------- end  clustering loop ------------------------------
        
        int MinClustSize=10;
        double MinClustWidth=0.001;
        double MinClustLength=0.01;
        double zStart =  5.; // mm
        double zEnd   = 29.; // mm
        int ii=0;
        #ifdef VERBOSE
          printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
          printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
        #endif
        
        for (int k=0; k<nclust; k++) {
          #ifdef VERBOSE
            printf("%2d Clust(%2d): %6.1f %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Ypos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
            printf("                                  %6.1f %6.1f %6.1f %6.1f \n",clust_Width[k][0],clust_Width[k][1], clust_Length[k][0], clust_Length[k][1]);
          #endif
          
          //-------------  Cluster Filter -----------------
         if (clust_Size[k] >= MinClustSize && zStart < clust_Zpos[k] && clust_Zpos[k] < zEnd && clust_Width[k][2]>MinClustWidth ) {
            hits_Xpos[ii]=clust_Xpos[k];
            hits_Ypos[ii]=clust_Ypos[k];
            hits_Zpos[ii]=clust_Zpos[k];
            hits_dEdx[ii]=clust_dEdx[k];
            ii++;
            //printf("\n");
          } else {
            #ifdef VERBOSE
  	          printf(" <--- skip \n");
            #endif
  	      }
        }
        int nhits=ii;
        // -----  end hist dist clustering  -----
        
        //======================== Draw HITS and CLUST ==================================
        #ifdef SHOW_EVTbyEVT
        char hevtTitle[80]; sprintf(hevtTitle," Display Event: %lld   Run: %d; z pos,mm; y pos,mm ",jentry,RunNum);
        hevt->SetTitle(hevtTitle);
          #ifdef VERBOSE
            printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
          #endif
          if (jentry<1000) {
          //if (jentry<=1000 && pion) {
            c2->cd(1);   hevt->Draw("colz");
            c2->cd(2);   hevtf->Draw("text");
      	    c2->cd(3);   hevti->Draw("colz");
            c2->Modified(); c2->Update();
            c2->cd(1); gPad->Modified(); gPad->Update();
      	    int COLMAP[]={1,2,3,4,6,5};
      	    int pmt=22 ,pmt0 = 20; // PM type
      	    for(int i = 0; i < nclust; i++) {
      	      TMarker m = TMarker(clust_Zpos[i],clust_Xpos[i],pmt);
      	      int tcol=2; //min(tracks[i],6);
      	      if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
      	      int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);
      	      m.SetMarkerSize(0.7+clust_dEdx[i]/300);
      	      m.DrawClone();  gPad->Modified(); gPad->Update();
      	    }
      	    c2->Modified(); c2->Update();
          }
        #endif
        
        #if (USE_GNN==1)   // GNN MC
//==========================================================
//           Send to Model simulation
//==========================================================
          
          #ifdef VERBOSE
            printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
          #endif
          std::vector<int> tracks(nhits, 0);
          std::vector<float> Xcl;
          std::vector<float> Zcl;
          Xcl.clear();
          Zcl.clear();
          for (int n=0; n<nhits; n++) {
            Xcl.push_back(hits_Xpos[n]);
            Zcl.push_back(hits_Zpos[n]);
          }
          doPattern(Xcl, Zcl, tracks);  //---- call GNN ---
          #ifdef VERBOSE
            printf("**> End Model simulation \n"); //===================================================
          #endif
          
          #ifdef SHOW_EVTbyEVT
            c2->cd(2); gPad->Modified(); gPad->Update();
            int COLMAP[]={1,2,3,4,6,5};
            for(int i = 0; i < tracks.size(); i++) {
              #ifdef VERBOSE
        	      printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
              #endif
              TMarker m = TMarker(hits_Zpos[i],hits_Xpos[i],24);
              int tcol=min(tracks[i],6);
          	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5);
              m.DrawClone();  gPad->Modified(); gPad->Update();
            }
            #ifdef VERBOSE
              printf("\n\n");
              printf("**> End Cluster Plot \n");
            #endif
          #endif
//==================================================
//----           Track fitting                 -----
//==================================================
          
          #ifdef VERBOSE
            printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());
          #endif
          //----------------- tracks sorting -------------
          
          std::vector<std::vector<float>> TRACKS;  // -- [Nnodes][N_params]
          TRACKS.resize(nhits);
          std::vector<float> hit_coord(2,0);
          std::vector<int>  TRACKS_N(nhits, 0); // [Nnodes]
          for (int i = 0; i < nhits; i++)  {
            TRACKS_N[i] = 0;
          }
          std::vector<float> xz(2,0);
          
          for (int i2 = 0; i2 < nhits; i2++) {
            int num =  tracks[i2];
            int num2 = std::max(0, std::min(num, nhits - 1));
            #ifdef VERBOSE
              printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
            #endif
            xz[0]=Xcl[i2]; 	xz[1]=Zcl[i2];
            TRACKS[num2].push_back(Xcl[i2]); TRACKS[num2].push_back(Zcl[i2]);
            TRACKS_N[num2]++;
          }
          #if (DEBUG > 1)
            for (int i2 = 0; i2 < nhits; i2++) {
              printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
              for (int i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
                printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
              }
              if (TRACKS_N[i2]>0) printf("\n");
            }
          #endif //--- end tracks sorting
          
          //-----------------------------------
          //---       LSTM fitting          ---
          //-----------------------------------
          
          static TMultiGraph *mg;
          if (mg != NULL ) delete mg;
          mg = new TMultiGraph();
          int NTRACKS=0;
          
          for (int i2 = 1; i2 < nhits; i2++) {  // tracks loop; zero track -> noise
            if (TRACKS_N[i2]<2) continue;   //---- select 2 (x,z) and more hits on track ----
            #ifdef VERBOSE
              printf("==> fit: start trk: %d \r\n", i2);
            #endif
            std::vector<Double_t> x;
            std::vector<Double_t> y;
            
            for (int i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
              #ifdef VERBOSE
                printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
              #endif
              x.push_back(TRACKS[i2].at(i3+1));
              y.push_back(TRACKS[i2].at(i3));
            }
            #ifdef SHOW_EVTbyEVT
              c2->cd(3);   hevt->Draw("colz");
              TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
              TF1 *f = new TF1("f", "[1] * x + [0]");
              g->Fit(f,"q0");
              mg->Add(g,"p");
            #endif
            NTRACKS++;
          }  //---  end tracks loop
          
          #ifdef SHOW_EVTbyEVT
            cout<<"Event#="<<event_num<<" Electron="<<electron<<"  Pion="<<pion<<" CherenkovSignal="<<electron_ch<<" CalorimeterSum="<<CalSum<<" #ofTracks="<<NTRACKS<<endl;
            c2->cd(3); mg->Draw("APsame");  gPad->WaitPrimitive();
            mg->GetXaxis()->SetLimits(0.,30);
            mg->SetMinimum(-50.);
            mg->SetMaximum(+50.);
            gPad->Modified(); gPad->Update();
          #endif
        #endif // ============= End if (USE_GNN==1) (MC) ==============
        
        #if (USE_TCP==1)
//----------------------------------------------------------
//---                 Send to FPGA                      ----
//----------------------------------------------------------
          
//-----------------  send DATA  ----------------
          #ifdef VERBOSE
            printf(" send DATA  \n");
          #endif
          int DC_NROC=4;
          int LEN_HITS=nhits; //-- floats X,Y,Z,E
          if (LEN_HITS>50) LEN_HITS=50; //--- max hits to fpga
          int k=3; //-- start data 3 words header;
          #ifdef VERBOSE
            printf("-- BUFFER:: \n");
          #endif
          for (int n=0; n<LEN_HITS; n++) {
            FBUFFER[k++]= hits_Xpos[n];
            #ifdef VERBOSE
              printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
            #endif
            FBUFFER[k++]= hits_Zpos[n];
            #ifdef VERBOSE
              printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
              printf("=====>  Clust data: X,Y,Z,E=%f %f %f %f \n",hits_Xpos[n],hits_Ypos[n],hits_Zpos[n],hits_dEdx[n]);
            #endif
          }
          #ifdef VERBOSE
            printf("Filled buffer size=%d , data = %d hits =%d \n",k,k-3,(k-3)/2);
          #endif
          itrg++;
          
          for (int i=0;i<nmod;i++) {
            LENEVENT=k;
            hdr=(i & 0xF) << 24;
            
            BUFFER[0]=DC_NROC; // DC_NROC if send to EVB, or
            BUFFER[0] |= EVT_Type_DATA;
            BUFFER[1]=itrg;
            BUFFER[2]=LENEVENT;
            
            BUFFER[0]&= ~(0xff << 24);  // --  clear ModID
            unsigned int MODIDx=((modID+i)&0xff)<<24;   //-- ModID=modID+i
            BUFFER[0] |= MODIDx;
            
            unsigned int evtTrigID=BUFFER[1];
            int evtModID=(BUFFER[0]>>24)&0xff;
            int evtSize=BUFFER[2];
            #ifdef VERBOSE
              if ( itrg<200) {
                printf("==> SEND:: Trg=%d(%d,%d) Mod=%d(%d) siz=%d(%d)\n",evtTrigID,itrg,BUFFER[1],evtModID,i,evtSize,LENEVENT);
              }
            #endif
            
            unsigned int *DATA = BUFFER;
            int lenDATA=LENEVENT;
            int n = nmod;
            int k = i;
            unsigned int evtHDR = hdr;
            unsigned int TriggerID = itrg;
            
            HEADER[0]=0x5;  //---  buffered for evb
            HEADER[1]=0xAABBCCDD;
            HEADER[2]=lenDATA;
            HEADER[3]=evtHDR;
            HEADER[4]=TriggerID;
            HEADER[5]=n;
            HEADER[6]=k;
            HEADER[7]=k;
            
            sock->SendRaw((char*) HEADER, sizeof(HEADER), kDefault);
            sock->SendRaw((char*) DATA, lenDATA*4, kDefault);
            #ifdef VERBOSE
              printf("read GNN out, wait for FPGA data ... \n");  //=======================================================
            #endif
            unsigned int NDATA[MAX_NODES+10];
            int RHEADER[10];
            int COLMAP[]={1,2,3,4,6,5};
            sock->RecvRaw((char*) RHEADER, sizeof(RHEADER), kDefault);
            int lenNODES=RHEADER[2];
            #ifdef VERBOSE
              printf("RHEADER::"); for (int ih=0; ih<5; ih++) printf(" %d \n",HEADER[ih]);  printf(" LenDATA=%d \n",HEADER[2]);
            #endif
            sock->RecvRaw((char*) NDATA, lenNODES*4, kDefault);
            int nnodes=lenNODES-3;
            #ifdef VERBOSE
              printf("nodes return: %d (nclust=%d), TRKS: \n", nnodes,nclust);
            #endif
        	  unsigned int TDATA[2048];
            float *FTDATA = (float*) TDATA;
            
            #if (USE_FIT==1)
              #ifdef VERBOSE
                printf("read FIT out, wait for FPGA data ... \n");
              #endif
              sock->RecvRaw((char*) RHEADER, sizeof(RHEADER), kDefault);
              int lenFITS=RHEADER[2];
              #ifdef VERBOSE
                printf("RHEADER::"); for (int ih=0; ih<5; ih++) printf(" %d \n",RHEADER[ih]);  printf(" LenDATA=%d \n",RHEADER[2]);
              #endif
              sock->RecvRaw((char*) TDATA, lenFITS*4, kDefault);
              #ifdef VERBOSE
                for (int i=0; i<lenFITS; i++) {
                  printf("tracks fit return: i=%d  data=0x%x (%f)  \n", i,TDATA[i],FTDATA[i]);
                }
              #endif
              int nfits=lenFITS-3;
              #ifdef VERBOSE
                printf("tracks fit return: %d  \n", nfits);
              #endif
            #else //--- if USE_FIT !=1
              int nfits=nnodes;
              #ifdef VERBOSE
                printf("tracks fit return: %d  \n", nfits);
              #endif
            #endif  // --- End if USE_FIT==1
            
            if (nfits>0) {
              //=============== Draw FPGA Clust =============
              #ifdef SHOW_EVTbyEVT
                for (int nd=0; nd<min(nnodes,nhits); nd++) {
                  int trknum=NDATA[nd+3];
                  printf(" %u, ",trknum);
                  c2->cd(3);	gPad->Modified(); gPad->Update();
            	      if (trknum>0) {
                    TMarker m = TMarker(hits_Zpos[nd],hits_Xpos[nd],24);
                    int tcol=min(trknum,6);
                    int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(20);     m.SetMarkerSize(1.5);
                    m.DrawClone();  gPad->Modified(); gPad->Update();
                    #ifdef VERBOSE
                      printf("=====>  Draw:%d X,Y,Z,E=%f %f %f %f \n",nd,hits_Xpos[nd],hits_Ypos[nd],hits_Zpos[nd],hits_dEdx[nd]);
                    #endif
                  }
                }
                c2->Modified(); c2->Update();
                printf("\n");
              #endif
              
              #if (USE_FIT==1)
                //=============== Draw Tracks lines =============
                #ifdef SHOW_EVTbyEVT
                  int ntracks = nfits/3;
                  int cs = nfits % 3;
                  if (cs != 0)  { printf("==========>>>>   Error FIT results : %d %d %d \n",nfits, ntracks,cs);  break; };
                  
                  int cnt = 3; // word counter in data buffer
                  c2->cd(3);  gPad->Modified(); gPad->Update();
                  for (int i=0; i<ntracks; i++) {
                    int trknum = TDATA[cnt++];	  float aa = FTDATA[cnt++];	  float bb = FTDATA[cnt++];
                    #ifdef VERBOSE
                      printf(" Fit Track=%d aa=%f bb=%f \n",trknum,aa,bb);
                    #endif
                    TF1 ftrk("ftrk","[0]*x+[1]",zStart,zEnd);	  ftrk.SetParameter(0,aa);	  ftrk.SetParameter(1,bb);
                    ftrk.DrawClone("same");  gPad->Modified(); gPad->Update();
                  }
                  c2->Modified(); c2->Update();
                #endif
              #endif  // --- end if USE_FIT==1
              
            } else {
              printf(" No tracks to draw \n");
            } //---end if nfits>0
            
            //=============== Draw All Clust ================
            //----------------------  To FPGA ------------------
            #ifdef SHOW_EVTbyEVT
              c2->cd(2); gPad->Modified(); gPad->Update();
              //printf(" Draw clusters  \n");
              for (int k=0; k<nhits; k++) {
                TMarker mh = TMarker(hits_Zpos[k],hits_Xpos[k],24);
                int mhcol = 1;   mh.SetMarkerColor(mhcol);  mh.SetMarkerSize(1.5);
                mh.DrawClone();  gPad->Modified(); gPad->Update();
              }
              c2->Modified(); c2->Update();
              //----------------------  from FPGA ------------------
              c2->cd(3); gPad->Modified(); gPad->Update();
              //printf(" Draw clusters  \n");
              for (int k=0; k<nhits; k++) {
                TMarker mh = TMarker(hits_Zpos[k],hits_Xpos[k],24);
                int mhcol = 1;   mh.SetMarkerColor(mhcol);  mh.SetMarkerSize(1.5);
                mh.DrawClone();  gPad->Modified(); gPad->Update();
              }
              c2->Modified(); c2->Update();
            #endif
          } // -- End nmod loop ---
        #endif  // --- End if USE_TCP  ---
        
        #ifdef VERBOSE
          printf(" all done, click low right pad ...  \n");
        #endif
      #endif   // --- End USE_GNN>0 ---
      
//============ END GEMTRD Pattern Recognition Tracking ==================
//===============  END Process Fa125 RAW Data Loop  ============++++======
      
      #endif //====================== END if (USE_125_RAW) =============================
      
      
      
//===========================================================
//  GEMTracker (SRS) Correlation with TRD Prototypes
//===========================================================
      
      ULong64_t gt_idx_x = 0, gt_idx_y=0, mmg1_idx_y=0, urw_idx_y=0, mmg2_idx_y=0;
      //double x0_urw=-999;
      //double x0_mmg1=-999;
      //double x0_mmg2=-999;
      //double x0_gem=-999;
      
      double gemtrkr_peak_pos_y[gem_peak_count];
      double gemtrkr_peak_pos_x[gem_peak_count];
      double mmg1_peak_pos_y[gem_peak_count];
      double urw_peak_pos_y[gem_peak_count];
      double mmg2_peak_pos_y[gem_peak_count];
      int gemtrkr_peak_ch_y[gem_peak_count];
      int gemtrkr_peak_ch_x[gem_peak_count];
      
      //if (NTRACKS==2 && electron) n_e_2++; Count("n_e_2");
      //if (NTRACKS==2 && pion) n_pi_2++; Count("n_pi_2");
      
      for (ULong64_t i=0; i<gem_peak_count; i++) {
        gemtrkr_peak_pos_y[gem_peak_count] = 0.;
        gemtrkr_peak_pos_x[gem_peak_count] = 0;
        mmg1_peak_pos_y[gem_peak_count] = 0;
        urw_peak_pos_y[gem_peak_count] = 0;
        mmg2_peak_pos_y[gem_peak_count] = 0;
        gemtrkr_peak_ch_y[gem_peak_count] = 0;
        gemtrkr_peak_ch_x[gem_peak_count] = 0;
      }
      
      for (ULong64_t i=0; i<gem_peak_count; i++) { //-- SRS Peaks Loop
        
        if (gem_peak_plane_name->at(i) == "GEMTRKX") {
          gemtrkr_peak_pos_y[gt_idx_y] = gem_peak_real_pos->at(i);
          if (gemtrkr_peak_pos_y[gt_idx_y]<=0) gemtrkr_peak_pos_y[gt_idx_y]+=50.; else gemtrkr_peak_pos_y[gt_idx_y]-=50.;  gemtrkr_peak_pos_y[gt_idx_y]*=-1.;
          gemtrkr_peak_ch_y[gt_idx_y] = gem_peak_index->at(i);
          if (gemtrkr_peak_ch_y[gt_idx_y]<128) gemtrkr_peak_ch_y[gt_idx_y]+=128; else gemtrkr_peak_ch_y[gt_idx_y]-=128;
          gt_idx_y++;
        } if (gem_peak_plane_name->at(i) == "GEMTRKY") {
            gemtrkr_peak_pos_x[gt_idx_x] = gem_peak_real_pos->at(i);
            if (gemtrkr_peak_pos_x[gt_idx_x]<=0) gemtrkr_peak_pos_x[gt_idx_x]+=50.; else gemtrkr_peak_pos_x[gt_idx_x]-=50.;  gemtrkr_peak_pos_x[gt_idx_x]*=-1.;
            gemtrkr_peak_ch_x[gt_idx_x] = gem_peak_index->at(i);
            if (gemtrkr_peak_ch_x[gt_idx_x]<128) gemtrkr_peak_ch_x[gt_idx_x]+=128; else gemtrkr_peak_ch_x[gt_idx_x]-=128;
            gt_idx_x++;
        } if (gem_peak_plane_name->at(i) == "MMG1Y") {
            mmg1_peak_pos_y[mmg1_idx_y] = gem_peak_real_pos->at(i);
            if (mmg1_peak_pos_y[mmg1_idx_y]<=0) mmg1_peak_pos_y[mmg1_idx_y]+=50.; else mmg1_peak_pos_y[mmg1_idx_y]-=50.;  mmg1_peak_pos_y[mmg1_idx_y]*=-1.;
            mmg1_idx_y++;
        } if (gem_peak_plane_name->at(i) == "URWELLY") {
            urw_peak_pos_y[urw_idx_y] = gem_peak_real_pos->at(i);
            urw_peak_pos_y[urw_idx_y]*=-1.;
            urw_idx_y++;
        } if (gem_peak_plane_name->at(i) == "MMG2Y") {
            mmg2_peak_pos_y[mmg2_idx_y] = gem_peak_real_pos->at(i);
            if (mmg2_peak_pos_y[mmg2_idx_y]<=0) mmg2_peak_pos_y[mmg2_idx_y]+=50.; else mmg2_peak_pos_y[mmg2_idx_y]-=50.;  mmg2_peak_pos_y[mmg2_idx_y]*=-1.;
            mmg2_idx_y++;
        }
      } //--- End SRS Peaks Loop
      
      for (ULong64_t j=0; j<=gt_idx_y; j++) {
        hgemtrkr_peak_y->Fill(gemtrkr_peak_pos_y[j]);
        for (ULong64_t i=0; i<mmg1_idx_y; i++) {
          srs_mmg1_y->Fill(gemtrkr_peak_pos_y[j], mmg1_peak_pos_y[i]);
        }
        for (ULong64_t h=0; h<urw_idx_y; h++) {
          srs_urw_y->Fill(gemtrkr_peak_pos_y[j], urw_peak_pos_y[h]);
        }
        //srs_mmg2_y->Fill(gemtrkr_peak_pos_y[gt_idx_y], mmg2_peak_pos_y[i]);
        for (ULong64_t k=0; k<=gt_idx_x; k++) {
          if (NTRACKS==1) {
            singleTrackIndex->Fill(gt_idx_x, gt_idx_y);
            if (electron) singleTrackIndex_e->Fill(gt_idx_x, gt_idx_y);
            if (pion) singleTrackIndex_pi->Fill(gt_idx_x, gt_idx_y);
          }
          if (NTRACKS==2) {
            multiTrackIndex->Fill(gt_idx_x, gt_idx_y);
            if (electron) multiTrackIndex_e->Fill(gt_idx_x, gt_idx_y);
            if (pion) multiTrackIndex_pi->Fill(gt_idx_x, gt_idx_y);
          }
          hgemtrkr_peak_xy->Fill(gemtrkr_peak_pos_x[k], gemtrkr_peak_pos_y[j]);
          hgemtrkr_ch_xy->Fill(gemtrkr_peak_ch_x[k], gemtrkr_peak_ch_y[j]);
          //srs_gem_dy->Fill(x0_gem, gemtrkr_peak_pos_y[j]);
          //srs_gem_dx->Fill(x0_gem, gemtrkr_peak_pos_x[k]);
          //srs_mmg1_dx->Fill(x0_mmg1, gemtrkr_peak_pos_x[k]);
          //srs_mmg1_dy->Fill(x0_mmg1, gemtrkr_peak_pos_y[j]);
          //srs_urw_dx->Fill(x0_urw, gemtrkr_peak_pos_x[k]);
          //srs_urw_dy->Fill(x0_urw, gemtrkr_peak_pos_y[j]);
        }
      }
      for (ULong64_t k=0; k<gt_idx_x; k++) {
        hgemtrkr_peak_x->Fill(gemtrkr_peak_pos_x[k]);
        for (ULong64_t i=0; i<mmg1_idx_y; i++) {
          srs_mmg1_x->Fill(gemtrkr_peak_pos_x[k], mmg1_peak_pos_y[i]);
        }
        for (ULong64_t h=0; h<urw_idx_y; h++) {
          srs_urw_x->Fill(gemtrkr_peak_pos_x[k], urw_peak_pos_y[h]);
        }
        //srs_mmg2_x->Fill(gemtrkr_peak_pos_x[i], mmg2_peak_pos_y[i]);
      }
    //}
  //for (ULong64_t i=0; i<mmg1_idx_y; i++) {
  //  mmg1_xy->Fill(x0_mmg1, mmg1_peak_pos_y[i]);
  //  gem_mmg1_y->Fill(x0_gem, mmg1_peak_pos_y[i]);
  //}
  //for (ULong64_t h=0; h<urw_idx_y; h++) {
  //  urw_xy->Fill(x0_urw, urw_peak_pos_y[h]);
  //  gem_urw_y->Fill(x0_gem, urw_peak_pos_y[h]);
  //  for (ULong64_t i=0; i<mmg1_idx_y; i++) {
  //    mmg1_urw_y->Fill(urw_peak_pos_y[h], mmg1_peak_pos_y[i]);
  //  }
  //}
/*  srs_mmg2_dx->Fill(x0_mmg2, gemtrkr_peak_pos_x[i]);
    srs_mmg2_dy->Fill(x0_mmg2, gemtrkr_peak_pos_y[i]);
    mmg2_xy->Fill(x0_mmg2, mmg2_peak_pos_y[i]);
*/
      hNTracks->Fill(NTRACKS);
      if (electron) hNTracks_e->Fill(NTRACKS);
      if (pion) hNTracks_pi->Fill(NTRACKS);
      
//==========================================================
//    GEM TRKR Fuducial Area (Y-Direction) Selection
//==========================================================
    
      bool fiducialArea=false;
      
      if ((RunNum > 3200 && RunNum < 3203) || (RunNum >  3245 && RunNum < 3249)) { //-- GEMTRD Double Fleece
        x_cut1=-50., x_cut2=50., y_cut1=-20., y_cut2=15.;
        for (ULong64_t i=0; i<=gt_idx_y; i++) {
          if (y_cut2>gemtrkr_peak_pos_y[gt_idx_y] && gemtrkr_peak_pos_y[gt_idx_y]>y_cut1) {
            fiducialArea=true; 
            break;
          }
        }
      }
      if (RunNum > 3124 && RunNum < 3132) { //-- GEMTRD Single Fleece
        x_cut1=-50., x_cut2=50., y_cut1=-25., y_cut2=25.;
        for (ULong64_t i=0; i<=gt_idx_y; i++) {
          if (y_cut2>gemtrkr_peak_pos_y[gt_idx_y] && gemtrkr_peak_pos_y[gt_idx_y]>y_cut1) {
            fiducialArea=true; 
            break;
          }
        }
      }
      if ((RunNum > 3131 && RunNum < 3135) || (RunNum >  3249 && RunNum < 3254)) { //-- GEMTRD Single Foil (VU)
        x_cut1=-50., x_cut2=50., y_cut1=-50., y_cut2=50.;
        fiducialArea=true;
      }
      if (RunNum > 3274 && RunNum < 3289) { //-- GEMTRD Single Foil (TU)
        x_cut1=-50., x_cut2=50., y_cut1=-40., y_cut2=40.;
        for (ULong64_t i=0; i<=gt_idx_y; i++) {
          if (y_cut2>gemtrkr_peak_pos_y[gt_idx_y] && gemtrkr_peak_pos_y[gt_idx_y]>y_cut1) {
            fiducialArea=true; 
            break;
          }
        }
      }
      if ((RunNum >  3195 && RunNum < 3201) || (RunNum >  3214 && RunNum < 3219)) { //-- GEMTRD Double Foil
        x_cut1=-50., x_cut2=50., y_cut1=-40., y_cut2=40.;
        for (ULong64_t i=0; i<=gt_idx_y; i++) {
          if (y_cut2>gemtrkr_peak_pos_y[gt_idx_y] && gemtrkr_peak_pos_y[gt_idx_y]>y_cut1) {
            fiducialArea=true; 
            break;
          }
        }
      }
      if (RunNum > 3202 && RunNum < 3205) { //-- GEMTRD No Radiator
        x_cut1=-50., x_cut2=50., y_cut1=-50., y_cut2=50.;
        fiducialArea=true;
      }
      
//======== END GEM TRKR Fuducial Area Selection =========
    
//============================================================================
//                  Single Track Event fADC Processing
//============================================================================
    
      bool isSingleTrack=false;
      if (NTRACKS==1) isSingleTrack=true;
      if (isSingleTrack) { Count("_1TRK"); _1TRK++; }
      if (isSingleTrack && electron) { Count("e_1TRK");  e_1TRK++; }
      if (isSingleTrack && pion) { Count("pi_1TRK");  pi_1TRK++; }
      
      //--- Reset Single Track Hit Info for each event
      gem_nhit=0;
      gem_xpos.clear();
      gem_zpos.clear();
      gem_dedx.clear();
      gem_parID.clear();
      gem_zHist->Reset();
      gem_zHist_vect.clear();
      
      mmg1_nhit=0;
      mmg1_xpos.clear();
      mmg1_zpos.clear();
      mmg1_dedx.clear();
      mmg1_parID.clear();
      mmg1_zHist->Reset();
      mmg1_zHist_vect.clear();
     
      if (RunNum>3261) {
        mmg2_nhit=0;
        mmg2_xpos.clear();
        mmg2_zpos.clear();
        mmg2_dedx.clear();
        mmg2_parID.clear();
        mmg2_zHist->Reset();
        mmg2_zHist_vect.clear();
      }
      
      if (RunNum<3262 && RunNum>3147) {
        urw_nhit=0;
        urw_xpos.clear();
        urw_zpos.clear();
        urw_dedx.clear();
        urw_parID.clear();
        urw_zHist->Reset();
        urw_zHist_vect.clear();
      }
      
      if (isSingleTrack && fiducialArea) {
        
        double gem_amp_max = 0.;
        double mmg1_amp_max = 0.;
        double urw_amp_max = 0.;
        double mmg2_amp_max = 0.;
        n_1++; Count("n_1");
        if (electron) { n_e_1++; Count("n_e_1"); }
        if (pion) { n_pi_1++; Count("n_pi_1"); }
        
        for (ULong64_t i=0;i<f125_pulse_count; i++) { //--- Fadc125 Pulse Loop
          
          float peak_amp = f125_pulse_peak_amp->at(i);
          float ped = f125_pulse_pedestal->at(i);
          if (0 > ped || ped > 200 ) ped = 100;
          float amp=peak_amp-ped;
          if (amp<0) amp=0;
          float time=f125_pulse_peak_time->at(i);
          int fADCSlot = f125_pulse_slot->at(i);
          int fADCChan = f125_pulse_channel->at(i);
          
          int gemChan = GetGEMChan(fADCChan, fADCSlot);
          int mmg1Chan = GetMMG1Chan(fADCChan, fADCSlot, RunNum);
          int mmg2Chan = GetMMG2Chan(fADCChan, fADCSlot, RunNum);
          int rwellChan = GetRWELLChan(fADCChan, fADCSlot, RunNum);
          
          if (electron) {
            if (amp>GEM_THR && gemChan>-1) {
              if (gem_amp_max<amp) gem_amp_max=amp;
              f125_el->Fill(amp);
              f125_el_amp2ds->Fill(time, gemChan, amp);
              gem_xpos.push_back(gemChan);
              gem_dedx.push_back(amp);
              gem_zpos.push_back(time);
              gem_parID.push_back(electron);
              gem_nhit++;
              gem_zHist->Fill(time, amp);
            }
            if (amp>MM_THR && mmg1Chan>-1) {
              if (mmg1_amp_max<amp) mmg1_amp_max=amp;
              mmg1_f125_el->Fill(amp);
              mmg1_f125_el_amp2ds->Fill(time, mmg1Chan, amp);
              mmg1_xpos.push_back(mmg1Chan);
              mmg1_dedx.push_back(amp);
              mmg1_zpos.push_back(time);
              mmg1_parID.push_back(electron);
              mmg1_nhit++;
              mmg1_zHist->Fill(time, amp);
            }
            if (RunNum>3261) {
              if (amp>MM_THR && mmg2Chan>-1) {
                if (mmg2_amp_max<amp) mmg2_amp_max=amp;
                mmg2_f125_el->Fill(amp);
                mmg2_f125_el_amp2ds->Fill(time, mmg2Chan, amp);
                mmg2_xpos.push_back(mmg2Chan);
                mmg2_dedx.push_back(amp);
                mmg2_zpos.push_back(time);
                mmg2_parID.push_back(electron);
                mmg2_nhit++;
                mmg2_zHist->Fill(time, amp);
              }
            } else if (RunNum<3262 && RunNum>3147) {
              if (amp>URW_THR && rwellChan>-1) {
                if (urw_amp_max<amp) urw_amp_max=amp;
                urw_f125_el->Fill(amp);
                urw_f125_el_amp2ds->Fill(time, rwellChan, amp);
                urw_xpos.push_back(rwellChan);
                urw_dedx.push_back(amp);
                urw_zpos.push_back(time);
                urw_parID.push_back(electron);
                urw_nhit++;
                urw_zHist->Fill(time, amp);
              }
            }
            
          } else if (pion) {
            if (amp>GEM_THR && gemChan>-1) {
              if (gem_amp_max<amp) gem_amp_max=amp;
              f125_pi->Fill(amp);
              f125_pi_amp2ds->Fill(time, gemChan, amp);
              gem_xpos.push_back(gemChan);
              gem_dedx.push_back(amp);
              gem_zpos.push_back(time);
              gem_parID.push_back(electron);
              gem_nhit++;
            }
            if (amp>MM_THR && mmg1Chan>-1) {
              if (mmg1_amp_max<amp) mmg1_amp_max=amp;
              mmg1_f125_pi_amp2ds->Fill(time, mmg1Chan, amp);
              mmg1_f125_pi->Fill(amp);
              mmg1_xpos.push_back(mmg1Chan);
              mmg1_dedx.push_back(amp);
              mmg1_zpos.push_back(time);
              mmg1_parID.push_back(electron);
              mmg1_nhit++;
            }
            if (RunNum>3261) {
              if (amp>MM_THR && mmg2Chan>-1) {
                if (mmg2_amp_max<amp) mmg2_amp_max=amp;
                mmg2_f125_pi->Fill(amp);
                mmg2_f125_pi_amp2ds->Fill(time, mmg2Chan, amp);
                mmg2_xpos.push_back(mmg2Chan);
                mmg2_dedx.push_back(amp);
                mmg2_zpos.push_back(time);
                mmg2_parID.push_back(electron);
                mmg2_nhit++;
              }
            } else if (RunNum<3262 && RunNum>3147) {
              if (amp>URW_THR && rwellChan>-1) {
                if (urw_amp_max<amp) urw_amp_max=amp;
                urw_f125_pi->Fill(amp);
                urw_f125_pi_amp2ds->Fill(time, rwellChan, amp);
                urw_xpos.push_back(rwellChan);
                urw_dedx.push_back(amp);
                urw_zpos.push_back(time);
                urw_parID.push_back(electron);
                urw_nhit++;
              }
            }
          }
          hCCor_ud->Fill(Ch_u,Ch_out);
        } //--- end Fa125 Pulse Loop ---
        
        if (electron) {
          if (gem_amp_max!=0.) f125_el_max->Fill(gem_amp_max);
          if (mmg1_amp_max!=0.) mmg1_f125_el_max->Fill(mmg1_amp_max);
          if (urw_amp_max!=0.) urw_f125_el_max->Fill(urw_amp_max);
          if (mmg2_amp_max!=0.) mmg2_f125_el_max->Fill(mmg2_amp_max);
        } else if (pion) {
          if (gem_amp_max!=0.) f125_pi_max->Fill(gem_amp_max);
          if (mmg1_amp_max!=0.) mmg1_f125_pi_max->Fill(mmg1_amp_max);
          if (urw_amp_max!=0.) urw_f125_pi_max->Fill(urw_amp_max);
          if (mmg2_amp_max!=0.) mmg2_f125_pi_max->Fill(mmg2_amp_max);
        }
        for (int i=1; i<21; i++) {
          gem_zHist_vect.push_back(gem_zHist->GetBinContent(i));
          mmg1_zHist_vect.push_back(mmg1_zHist->GetBinContent(i));
          mmg2_zHist_vect.push_back(mmg2_zHist->GetBinContent(i));
          urw_zHist_vect.push_back(urw_zHist->GetBinContent(i));
        }
      } //=============== END Single Track Event fADC Processing =================
      
//=====================================================================================
//                        E v e n t    D i s p l a y
//=====================================================================================
  
    #ifdef SHOW_EVT_DISPLAY
      if (jentry<MAX_PRINT || !(jentry%NPRT)) {
        c0->cd(1); f125_el_amp2ds->Draw("colz");
        c0->cd(5); f125_pi_amp2ds->Draw("colz");
        //c0->cd(2); f125_el_evt_display->Draw("colz");
        //c0->cd(6); f125_pi_evt_display->Draw("colz");
        
        //c0->cd(3); f125_el_raw->Draw("colz");  f125_el_evt_display->Draw("same");
        //TLine lin1(110.,gemtrkr_x2ch,190.,gemtrkr_x2ch); lin1.Draw("same");   //--- draw  gemtrkr x
        ////printf("++++++++++++ Draw GEMTRK:: %f %f %f  \n",gemtrkr_x,ftrk.Eval(gemtrkr_x),gemtrkr_x2ch);
        //c0->cd(7); f125_pi_raw->Draw("colz");  f125_pi_evt_display->Draw("same");
        //TLine lin2(110.,gemtrkr_x2ch,190.,gemtrkr_x2ch); lin2.SetLineColor(kRed); lin2.Draw("same");   //--- draw  gemtrkr x
        /*
        if (electron || pion ) {
          int lc=0;
          for (int k=0; k<100; k++) {
            peak_line[k].SetX1 (100.);
            peak_line[lc].SetY1 (-10.);
            peak_line[lc].SetX2 (101.);
            peak_line[lc].SetY2 (-10.);
          }
          for (ULong64_t i=0; i<gem_peak_count; i++) {
            double pos = gem_peak_real_pos->at(i);
            if (pos<=0) pos=pos+50.; else pos=pos-50.;  pos*=-1.;
            double pos2ch=(ftrk.Eval(pos)+50.)/0.4;  // -- to gemtrd coordinate system
            peak_line[lc].SetX1 (110.);	    peak_line[lc].SetY1 (pos2ch);	    peak_line[lc].SetX2 (190.);	    peak_line[lc].SetY2 (pos2ch);
            if (gem_peak_plane_name->at(i) == "GEMTRKY" ) {
              if (lc<100) {
                peak_line[lc].SetLineColor(kGreen);
                peak_line[lc].Draw("same");
                lc++;
              }
            }
          }  //--- SRS Peak Loop
        }
        */
        c0->cd(9); hevt->Draw("colz");
        c0->cd(10); hevtc->Draw("colz");
        c0->Modified();   c0->Update();
        
        //---- Show Fiducial Area ----
        /*
        srs_gemtrd_el->Divide(srs_etrd_beam);
        //srs_etrd_ratio = (TH2F*)srs_gemtrd_el->Clone("srs_etrd_ratio");
        //srs_gemtrd_el->Copy((TH2F*)srs_etrd_ratio);
        //srs_etrd_ratio->GetXaxis()->SetTitle(" "); srs_etrd_ratio->GetYaxis()->SetTitle(" ");
        //srs_etrd_ratio->SetTitle("TR energy norm");  //srs_etrd_ratio->SetStats(1); srs_etrd_ratio->SetMinimum(0.8); srs_etrd_ratio->SetMaximum(1.35);
        #if 0
          srs_etrd_ratio->Divide(srs_etrd_beam);
        #else
          //srs_etrd_ratio->Add(srs_etrd_beam,-1.);
        #endif
        
        c0->cd(11); srs_etrd_ratio->DrawCopy("colz");
        c0->Modified();   c0->Update();
        //c0->cd(11); srs_gemtrd_el->Draw("colz");
        TBox fbox(x_cut1,y_cut1,x_cut2,y_cut2);  //---- draw box cut ---
        fbox.SetLineColor(kRed);
        fbox.SetFillStyle(0);
        fbox.SetLineWidth(2);
        fbox.DrawClone();
        */
        //c0->cd(12); srs_gem_dx->Draw("colz");
        c0->Modified();   c0->Update();
        c0->cd(12); gPad->WaitPrimitive();
      } //--- End jentry loop
    #endif //==================== END Event Display ========================
    
    #ifdef SAVE_TRACK_HITS
      //==== Fill Track Hit Info Trees ====
      if (gem_nhit>0) EVENT_VECT_GEM->Fill();
      if (mmg1_nhit>0) EVENT_VECT_MMG1->Fill();
      if (mmg2_nhit>0 && RunNum>3261) EVENT_VECT_MMG2->Fill();
      if (urw_nhit>0 && RunNum<3262 && RunNum>3147) EVENT_VECT_URW->Fill();
    #endif
  }
  //======================== End of Event Loop  ======================
  cout<<"========== End of Event Loop Reached - Last Evt:"<<event_num<<" =========="<<endl;
  
  cout<<"Total events= "<<jentry<<endl;
  cout<<"hcount values: 1_TRK="<<_1TRK<<" 1eTRK="<<e_1TRK<<" 1piTRK="<<pi_1TRK<<" elCC="<<el_CC<<" piCC="<<pi_CC<<" n_e_1="<<n_e_1<<" n_pi_1="<<n_pi_1<<" n_e_2="<<n_e_2<<" n_pi_2="<<n_pi_2<<endl;

//=====================================================================================
//===                 S A V E   H I S T O G R A M S                                ====
//=====================================================================================

  printf("Writing .root output file... \n");
  TFile* fOut;
  char rootFileName[256]; sprintf(rootFileName, "RootOutput/Run_%06d_Output.root", RunNum);
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  HistList->Write("HistDQM", TObject::kSingleKey);
  fOut->Close();
  printf("File written & closed OK \n");
  delete fOut;
  //HistList->Delete();
  
//=====================================================================================
//===        S A V E   S I N G L E   T R A C K   H I T   T T R E E S               ====
//=====================================================================================
  
  #ifdef SAVE_TRACK_HITS
    printf("Writing Hit Info TTree files... \n");
    fHits->cd();
    EVENT_VECT_GEM->Write();
    EVENT_VECT_MMG1->Write();
    if (RunNum>3261) {
      EVENT_VECT_MMG2->Write();
    } else if (RunNum<3262 && RunNum>3147) {
      EVENT_VECT_URW->Write();
    }
    fHits->Close();
    printf("TTree files written & closed OK \n");
    //fHits->Delete();
  #endif
  
//=====================================================================================
//===                 P L O T     H I S T O G R A M S                               ===
//=====================================================================================
  #ifdef SAVE_PDF
    const char *OutputDir="RootOutput";
    char ctit[120];
    sprintf(G_DIR,"%s/Run_%06d",OutputDir,RunNum);
    sprintf(ctit,"File=%s",G_DIR);
    bool COMPACT=false;
    TCanvas *cc;
    int nxd=3;
    int nyd=5;
    
    //======= Plot & Save Event Display ========
    //char pngname[120];  sprintf(pngname,"%s_evdisp.png",G_DIR);  //c0->Print(pngname);
    char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  //c0->Print(pdfname);
    
    //--------------------- new page --------------------
    htitle(" Cherenkov (Fadc250)  ");   //if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hcount->Draw();
    cc=NextPlot(nxd,nyd);  gPad->SetLogy();  hCher_u_adc->Draw();
    cc=NextPlot(nxd,nyd);  gPad->SetLogy();  hCher_dout_adc->Draw();
    cc=NextPlot(nxd,nyd);  hCCor_ud->Draw("colz");
    //--------------------- new page --------------------
    htitle(" GEM-Tracker (SRS) ");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  hgemtrkr_peak_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_ch_xy->Draw("colz");
    cc=NextPlot(nxd,nyd);  hgemtrkr_peak_x->Draw();
    cc=NextPlot(nxd,nyd);  hgemtrkr_peak_y->Draw();
    cc=NextPlot(nxd,nyd);  mmg1_xy->Draw("colz");
    if (RunNum<3262 && RunNum>3147) {urw_xy->Draw("colz");} else if (RunNum>3261) {mmg2_xy->Draw("colz");}
    /*  cc=NextPlot(nxd,nyd);  srs_gemtrd_el->Draw("colz");
    TBox fbox(x_cut1,y_cut1,x_cut2,y_cut2);  //---- draw box cut ---
    fbox.Draw("same");
    fbox.SetLineColor(kRed);
    fbox.SetFillStyle(0);
    fbox.SetLineWidth(1);
    cc=NextPlot(nxd,nyd);  srs_etrd_beam->Draw("colz");
    cc=NextPlot(nxd,nyd);  srs_etrd_ratio->Draw("colz");
    */
    //--------------------- new page --------------------
    htitle("  TRD Prototype (Fadc125) Amplitudes ");    if (!COMPACT) cc=NextPlot(0,0);
    cout<<"Amplitudes htitle..."<<endl;
    nxd=2; nyd=4;
    cout<<"Amplitudes nxd nyd OK"<<endl;
    cc=NextPlot(nxd,nyd); cout<<"NextPlot OK"<<endl;   gPad->SetLogy(); cout<<"SetLogY OK"<<endl;  f125_el->Draw(); cout<<"Draw _el OK"<<endl;
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();   f125_pi->Draw();
    cout<<"Amplitudes htitle _el OK"<<endl;
    cc=NextPlot(nxd,nyd);   f125_el_max->Draw();
    cc=NextPlot(nxd,nyd);   f125_pi_max->Draw();
    cout<<"Amplitudes htitle _max OK"<<endl;
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();   mmg1_f125_el->Draw();
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();   mmg1_f125_pi->Draw();
    cout<<"Amplitudes htitle _el mmg OK"<<endl;
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();   if (RunNum<3262 && RunNum>3147) {urw_f125_el->Draw();} else if (RunNum>3261) {mmg2_f125_el->Draw();}
    cc=NextPlot(nxd,nyd);   gPad->SetLogy();  if (RunNum<3262 && RunNum>3147) {urw_f125_pi->Draw();} else if (RunNum>3261) {mmg2_f125_pi->Draw();}
    cout<<"Amplitudes htitle Draw OK"<<endl;
    //--------------------- new page --------------------
    htitle("  TRD Prototype (Fadc125) Amplitudes - 2D");    if (!COMPACT) cc=NextPlot(0,0);
    nxd=2; nyd=3;
    cc=NextPlot(nxd,nyd);   f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   f125_pi_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   mmg1_f125_pi_amp2ds->Draw("colz");
    cc=NextPlot(nxd,nyd);   if (RunNum<3262 && RunNum>3147) {urw_f125_el_amp2ds->Draw("colz");} else if (RunNum>3261) {mmg2_f125_el_amp2ds->Draw("colz");}
    cc=NextPlot(nxd,nyd);   if (RunNum<3262 && RunNum>3147) {urw_f125_pi_amp2ds->Draw("colz");} else if (RunNum>3261) {mmg2_f125_pi_amp2ds->Draw("colz");}
    cout<<"2D Amplitudes htitle OK"<<endl;
    //--------------------- new page --------------------
    htitle(" SRS & TRD Prototypes - GEMTRD Correlations");    if (!COMPACT) cc=NextPlot(0,0);
    
    cc=NextPlot(nxd,nyd);  gem_mmg1_x->Draw("colz");  ftrk.Draw("same");
    if (RunNum<3262 && RunNum>3147) {cc=NextPlot(nxd,nyd); gem_urw_x->Draw("colz"); ftrk.Draw("same");} else if (RunNum>3261) {cc=NextPlot(nxd,nyd); gem_mmg2_x->Draw("colz"); ftrk.Draw("same");}
    cc=NextPlot(nxd,nyd);  gem_mmg1_y->Draw("colz");
    if (RunNum<3262 && RunNum>3147) { cc=NextPlot(nxd,nyd);  gem_urw_y->Draw("colz");} else if (RunNum>3261) { cc=NextPlot(nxd,nyd);  gem_mmg2_y->Draw("colz");}
    cc=NextPlot(nxd,nyd);  srs_gem_dy->Draw("colz");
    cout<<"SRS&TRD GEMTRD Correlations htitle OK"<<endl;
    //--------------------- new page --------------------
    htitle(" SRS & TRD Prototypes - Tracking");   if (!COMPACT) cc=NextPlot(0,0);
    cc=NextPlot(nxd,nyd);  singleTrackIndex->Draw("colz text");
    cc=NextPlot(nxd,nyd);  multiTrackIndex->Draw("colz text");
    cc=NextPlot(nxd,nyd);  srs_gem_dx->Draw("colz"); ftrk.Draw("same");
    cc=NextPlot(nxd,nyd);  srs_mmg1_dx->Draw("colz"); ftrk.Draw("same");
    cc=NextPlot(nxd,nyd);  srs_mmg1_dy->Draw("colz");
    if (RunNum<3262 && RunNum>3147) { cc=NextPlot(nxd,nyd); srs_urw_dx->Draw("colz"); ftrk.Draw("same"); cc=NextPlot(nxd,nyd); srs_urw_dy->Draw("colz");}
    else if (RunNum>3261) { cc=NextPlot(nxd,nyd); srs_mmg2_dx->Draw("colz"); ftrk.Draw("same"); cc=NextPlot(nxd,nyd); srs_mmg2_dy->Draw("colz");}
    cout<<"SRS&TRD Tracking htitle OK"<<endl;
    //--------------------- new page --------------------
    htitle(" SRS & TRD Prototypes - Y Correlations");   if (!COMPACT) cc=NextPlot(0,0);
    
    cc=NextPlot(nxd,nyd);  srs_mmg1_y->Draw("colz"); ftrk.Draw("same");
    cc=NextPlot(nxd,nyd);  srs_mmg1_x->Draw("colz");
    if (RunNum<3262 && RunNum>3147) {
      cc=NextPlot(nxd,nyd);  srs_urw_x->Draw("colz");
      cc=NextPlot(nxd,nyd);  srs_urw_y->Draw("colz"); ftrk.Draw("same");
      cc=NextPlot(nxd,nyd);  mmg1_urw_y->Draw("colz"); ftrk.Draw("same");
      
    } else if (RunNum>3261) {
      cc=NextPlot(nxd,nyd);  srs_mmg2_x->Draw("colz");
      cc=NextPlot(nxd,nyd);  srs_mmg2_y->Draw("colz"); ftrk.Draw("same");
    }
    //--- close PDF file ----
    cc=NextPlot(-1,-1);
    cout<<"Close PDF OK"<<endl;
  #endif
  cout<<"=========== END OF RUN "<<RunNum<<" ============"<<endl;
}
