#define trdclass_cxx
#include "trdclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TSocket.h"
#include "TMarker.h"
#include "TMultiGraph.h"
#include "PlotLib.C"
#include "GNN/gnn_model.h"
#include "GNN/gnn_model.cpp"
#include "GNN/toGraph.cpp"

#define NPRT 1
//#define NPRT 1
#define USE_TRK
#define  MAX_PRINT 1

#define USE_GNN  1
#define USE_TCP  0
#define USE_FIT  0

#define DEBUG    5

//---------  Large prototype (50x70) Mapping ----------

#define first_slot 3
#define gem_x_slot 0
#define gem_x_ch0 0

#define gem_y_slot 3
#define gem_y_ch0 24

#define last_slot 6
#define last_ch 47

int GetGEMXChan(int slot, int fch)
{
  int sl=slot-first_slot;
  if((sl>=gem_x_slot&&sl<gem_y_slot)||(sl==gem_y_slot&&fch<gem_y_ch0)){
    int ch=sl*72+fch-gem_x_ch0;
    return ch;
  } else {
    return -1;
  }

}
int GetGEMYChan(int slot, int fch)
{
  int sl=slot-first_slot;
  if((sl>=gem_y_slot&&sl<last_slot)||(sl==last_slot&&fch<=last_ch)){
    int ch=(sl-gem_y_slot)*72+fch-gem_y_ch0;
    return ch;
  } else {
    return -1;
  }

}
//--------------  FermiLab test mapping -----------------


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
   //+++++++++++++++++++++++++++++ TCP +++++++++++++++++++

#if  (USE_TCP==1)  //--  tcp send

	// Open connection to server
	TSocket *sock = new TSocket("localhost", 20250);

	// Wait till we get the start message
	char str[32];
	//sock->Recv(str, 32);

	// server tells us who we are
	//int idx = !strcmp(str, "go 0") ? 0 : 1;

	//printf("recv: %d %s \n",idx,str);
	Float_t messlen  = 0;
	Float_t cmesslen = 0;

#define BUFSIZE 128000
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

  //========= Book Histograms =============

   TList *HistList = new TList();

  //-----------------  canvas 0 Event Display ----------
  char c0Title[256]; sprintf(c0Title,"Event_Display_Run=%d",RunNum);
  TCanvas *c0 = new TCanvas("DISP",c0Title,1100,200,1500,1300);
  c0->Divide(4,3); c0->cd(1);
  //-----------------  canvas 1 FPGA Display ----------
  char c2Title[256]; sprintf(c2Title,"Event_Display_Run=%d",RunNum);
  TCanvas *c2 = new TCanvas("FPGA",c2Title,100,100,1000,1300);
  c2->Divide(1,3); c2->cd(1);

  TF1 fx("fx","pol1",100,190);
  TF1 fx1("fx1","pol1",100,190);
  TF1 fx2("fx2","pol1",100,190);

  //-- GEMTRD - GEMTRK alignment --------
  double xx1=-37., yy1=-55.,  xx2=53., yy2=44.;
  double aa=(yy2-yy1)/(xx2-xx1); 
  double bb=yy1-aa*xx1;
  TF1 ftrk("ftrk","[0]*x+[1]",-55.,55.);
  ftrk.SetParameter(0,aa);
  ftrk.SetParameter(1,bb);
  TF1 ftrkr("ftrk","(x-[1])/[0]",0.,255.);
  ftrkr.SetParameter(0,aa);
  ftrkr.SetParameter(1,bb);

  //----------------  GEM TRK fuducial area selection (box cut) ----------------------
  double xbc1=-50., xbc2=+50., ybc1=-50., ybc2=+50.; 
  double gemtrk_x2ch=-999.;
  TLine peak_line[100];
  //----------------------------------------------------------------------------------


  hCal_occ = new TH1F("hCal_occ"," Calorimeter Occupancy",9,-0.5,8.5);         HistList->Add(hCal_occ);
  hCal_sum = new TH1F("hCal_sum"," Calorimeter Sum (GeV)",100.,0.,15.);        HistList->Add(hCal_sum);
  for (int cc=0; cc<NCAL; cc++) {
    char hName[128];  sprintf(hName,"hCal_adc%d",cc);
    char hTitle[128]; sprintf(hTitle,"Calorimeter ADC, cell%d",cc);
    hCal_adc[cc] = new TH1F(hName,hTitle,4096,-0.5,4095.5);                    HistList->Add(hCal_adc[cc]);
    //
    sprintf(hName,"hCal_cor%d",cc);  sprintf(hTitle,"Calorimeter ADC corr GemTrd X, cell%d",cc);   
    hCal_cor[cc] =  new TH2F(hName,hTitle,25,0.5,25.5,100,-55.,55.);          HistList->Add(hCal_cor[cc]);
    //
    sprintf(hName,"hCal_trk%d",cc);  sprintf(hTitle,"Calorimeter corr GemTRK, cell%d",cc);   
    hCal_trk[cc] =  new TH2F(hName,hTitle,100,-55.,55.,100,-55.,55.);          HistList->Add(hCal_trk[cc]);
    //
    sprintf(hName,"hCal_cal%d",cc);  sprintf(hTitle,"Calorimeter ADC Calib, cell%d",cc);   
    hCal_cal[cc] = new TH2F(hName,hTitle,10,-0.5,4095.5,10,-5.,15.);           HistList->Add(hCal_cal[cc]);
  }
  cal_el_evt = new TH2F("cal_el_evt"," CAL el event ; X ; Y ",3,-0.5,2.5,3,-0.5,2.5);   HistList->Add(cal_el_evt);
  cal_el_evt->SetMinimum(-2.); cal_el_evt->SetMaximum(10.); cal_el_evt->SetStats(0);
  cal_pi_evt = new TH2F("cal_pi_evt"," CAL pi event ; X ; Y ",3,-0.5,2.5,3,-0.5,2.5);   HistList->Add(cal_pi_evt);
  cal_pi_evt->SetMinimum(-2.); cal_pi_evt->SetMaximum(10.); cal_pi_evt->SetStats(0);

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
  f125_el_max = new TH1F("f125_el_max","GEM-TRD f125 Max Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(f125_el_max);
  f125_pi_max = new TH1F("f125_pi_max","GEM-TRD f125 Max Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);               HistList->Add(f125_pi_max);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg1_f125_el);
  mmg1_f125_pi = new TH1F("mmg1_f125_pi","MMG1-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg1_f125_pi);
  urw_f125_el = new TH1F("urw_f125_el","uRW-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);          HistList->Add(urw_f125_el);
  urw_f125_pi = new TH1F("urw_f125_pi","uRW-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);              HistList->Add(urw_f125_pi);
  mmg2_f125_el = new TH1F("mmg2_f125_el","MMG2-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);       HistList->Add(mmg2_f125_el);
  mmg2_f125_pi = new TH1F("mmg2_f125_pi","MMG2-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);           HistList->Add(mmg2_f125_pi);

  f125_el_chi2 = new  TH1F("f125_el_chi2","TRK_el_chi2",100,0.,10000.);                            HistList->Add(f125_el_chi2);
  f125_pi_chi2 = new  TH1F("f125_pi_chi2","TRK_pi_chi2",100,0.,10000.);                            HistList->Add(f125_pi_chi2);
  f125_el_fita = new  TH1F("f125_el_fita","TRK_el_fita",100,-0.1,+0.1);                            HistList->Add(f125_el_fita);
  f125_pi_fita = new  TH1F("f125_pi_fita","TRK_pi_fita",100,-0.1,+0.1);                            HistList->Add(f125_pi_fita);
  
  //---- GEM-TRK --


  srs_ncl = new TH1F("srs_ncl"," Number SRS clusters per event",10,-0.5,9.5);                     HistList->Add(srs_ncl); 
  srs_trk_el = new TH2F("srs_trk_el","GEM-TRK , Electrons ; X ; Y ",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_trk_el);
  //srs_trk_el_->SetStats(0);  srs_trd_el->SetMinimum(THRESH);  srs_trk_el->SetMaximum(1000.);
  srs_trk_pi = new TH2F("srs_trk_pi","GEM-TRK , Pions ; X ; Y ",100,-55.,55.,100,-55.,55.);        HistList->Add(srs_trk_pi);
  //srs_trk_pi->SetStats(0);  srs_trk_pi->SetMinimum(THRESH);  srs_trk_pi->SetMaximum(1000.);

  srs_gem_dx = new TH2F("srs_gem_dx","Correlation TRD-TRK dX ; X trk; dX ",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_gem_dx);
  srs_gem_x = new TH2F("srs_gem_x","Correlation TRD-TRK X ; X trk; X chan",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_gem_x);
  srs_gem_y = new TH2F("srs_gem_y","Correlation TRD-TRK Y ; Y trk; X chan",100,-55.,55.,100,-55.,55.);    HistList->Add(srs_gem_y);

  srs_cal_corr = new TH2F("srs_cal_corr","Correlation TRK-CAL; X ; Y ",100,-55.,55.,100,-55.,55.);            HistList->Add(srs_cal_corr);
  srs_etrd_corr = new TH2F("srs_etrd_corr","Correlation TRK-energy TRD; X ; Y ",100,-55.,55.,100,-55.,55.);   HistList->Add(srs_etrd_corr);
  srs_etrd_beam = new TH2F("srs_etrd_beam","Correlation TRK beam; X ; Y ",100,-55.,55.,100,-55.,55.);         HistList->Add(srs_etrd_beam);
  srs_etrd_pion = new TH2F("srs_etrd_pion","Correlation TRK pion; X ; Y ",100,-55.,55.,100,-55.,55.);         HistList->Add(srs_etrd_pion);
  //srs_etrd_ratio = new TH2F("srs_etrd_ratio","Correlation TRK ratio; X ; Y ",100,-55.,55.,100,-55.,55.);         HistList->Add(srs_etrd_ratio);

  //---- GEM-TRD --
  int THRESH=0; 
  f125_el_evt = new TH2F("f125_el_evt","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);  HistList->Add(f125_el_evt);
  f125_el_evt->SetStats(0);  f125_el_evt->SetMinimum(THRESH);  f125_el_evt->SetMaximum(1000.);
  f125_pi_evt = new TH2F("f125_pi_evt","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);      HistList->Add(f125_pi_evt);
  f125_pi_evt->SetStats(0); f125_pi_evt->SetMinimum(THRESH); f125_pi_evt->SetMaximum(1000.);
  f125_el_raw = new TH2F("f125_el_raw","GEM-TRD raw for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);    HistList->Add(f125_el_raw);
  f125_el_raw->SetStats(0);  f125_el_raw->SetMinimum(THRESH);   f125_el_raw->SetMaximum(1000.);
  f125_pi_raw = new TH2F("f125_pi_raw","GEM-TRD raw for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,200,-0.5,249.5);        HistList->Add(f125_pi_raw);
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

  //---  clustering histograms  ---

    int nx0=100;    int ny0=250;
    hevt  = new TH2F("hevt"," Event display; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.); hevt->SetStats( 0 ); hevt->SetMaximum(10.);         HistList->Add(hevt);
    hevtc = new TH2F("hevtc"," Clustering ; FADC bins; GEM strips",nx0,-0.5,nx0-0.5,ny0,-0.5,ny0-0.5);                                 HistList->Add(hevtc);
    hevtc->SetStats(0);   hevtc->SetMinimum(0.07); hevtc->SetMaximum(40.);
    hevti = new TH2F("hevti"," ML-FPGA response; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.);  hevti->SetStats( 0 ); hevti->SetMaximum(10.);  HistList->Add(hevti);
    hevtf = new TH2F("hevtf"," Clusters for FPGA ; z pos,mm; y pos,mm ",nx0,0.,+30.,ny0,-50.,+50.);  hevtf->SetStats( 0 ); hevtf->SetMaximum(10.);  HistList->Add(hevtf);
 
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

  double Ebeam_el=0.4*Ebeam;
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
  printf("***>>>  Begin Event Loop 1st evt=%lld, MaxEvt=%lld \n",FirstEvt,MaxEvt);

  int N_trk_el=0;
  int N_trk_pi=0;
  Long64_t jentry=0;
  
  for (jentry=FirstEvt; jentry<nentries; jentry++) { //-- Event Loop --
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

      printf("-------------------- SRS   ev=%lld   ---------------------------\n",jentry);
      printf("SRS:: gem_scluster:  cnt=%lld \n",gem_scluster_count);
      for (ULong64_t i=0;i<gem_scluster_count; i++) { // --- SRS cluster loop
	printf("SRS:clusters:  i=%lld  X=%f Y=%f  \n",i,gem_scluster_x->at(i), gem_scluster_y->at(i));	      
      }
      printf("SRS:: srs_prerecon:  cnt=%lld cnt_x=%ld cnt_y=%ld \n",srs_prerecon_count, srs_prerecon_x->size(), srs_prerecon_y->size() );
      for (ULong64_t i=0;i<srs_prerecon_count; i++) { 
	printf("SRS:: srs_prerecon: i=%lld %f %f \n",i,srs_prerecon_x->at(i),srs_prerecon_y->at(i));
      }
      
      printf("SRS:0: srs_peak:  cnt=%lld \n",gem_peak_count);
      for (ULong64_t i=0;i<gem_peak_count; i++) { 
	printf("SRS:: srs_peak: i=%lld id=%d name=%s idx=%d apv=%d Amp=%f wid=%f E=%f Pos=%f \n"
	       ,i,gem_peak_plane_id->at(i),gem_peak_plane_name->at(i).c_str(),gem_peak_index->at(i), gem_peak_apv_id->at(i), gem_peak_height->at(i)
	       ,gem_peak_width->at(i), gem_peak_area->at(i), gem_peak_real_pos->at(i));
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
      if(electron_ch) {
	f125_el_evt->Reset();
	f125_el_raw->Reset();
	cal_el_evt->Reset();
	for (int cc=0; cc<NCAL; cc++) { int ix=cc%3; int iy=cc/3;
	  if (cc<6) cal_el_evt->Fill(ix,iy,Ecal[cc]); else cal_el_evt->Fill(1,2,Ecal[cc]); 
	  //double c = cal_el_evt->GetBinContent(ix+1, iy+1); printf("Ecal[%d]=%f ix=%d iy =%d e=%f \n",cc,Ecal[cc],ix,iy,c); 
	} 
	
      } else  {
	f125_pi_evt->Reset();
	f125_pi_raw->Reset();
	cal_pi_evt->Reset();	for (int cc=0; cc<NCAL; cc++) { if (cc<6) cal_pi_evt->Fill(cc%3,cc/3,Ecal[cc]); else cal_pi_evt->Fill(1,2,Ecal[cc]); } 
      }
    }
    if (electron_ch) {
      f125_el_fit->Reset();
      f125_el_amp2ds->Reset();
    } else {
      f125_pi_fit->Reset();
      f125_pi_amp2ds->Reset();
    }

    //==================================================================================================
    //                    SRS GemTRK 
    //==================================================================================================

    double gemtrk_x=-999.,  gemtrk_y=-999.,  gemtrk_E=0,  delta_x=1000.,  dx_thresh=5.;
    if (gem_scluster_count==1) {                         //--- use (first or single ?? ----
      for (ULong64_t ic=0; ic<gem_scluster_count; ic++) {   // --- SRS cluster loop, actually only 0 ;
	double x=gem_scluster_x->at(ic); if (x<=0) gemtrk_x=x+50.; else gemtrk_x=x-50.;  gemtrk_x*=-1.;
	double y=gem_scluster_y->at(ic); if (y<=0) gemtrk_y=y+50.; else gemtrk_y=y-50.;  gemtrk_y*=-1.;
	gemtrk_E=gem_scluster_energy->at(ic);
      }
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

      //if (jentry<MAX_PRINT) printf("F125:: i=%lld  sl=%d, ch=%d, npk=%d time=%d amp=%d ped=%d \n"
      //			   ,i,f125_pulse_slot->at(i),f125_pulse_channel->at(i),f125_pulse_npk->at(i)
      //			   ,f125_pulse_peak_time->at(i),f125_pulse_peak_amp->at(i),f125_pulse_pedestal->at(i));
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

      // ------- TR Radiator search -----
      if (gemChan>-1  && 100. < time && time < 185. ) { 
	srs_etrd_beam->Fill(gemtrk_x,gemtrk_y,1.);
	if (electron_ch) {
	  srs_etrd_corr->Fill(gemtrk_x,gemtrk_y,amp); 
	} else {  // --- pion ---
	  srs_etrd_pion->Fill(gemtrk_x,gemtrk_y,amp); 
	}
      }

      if(electron_ch) {  //-----------------  electron ---------------
	if (gemChan>-1) {
	  if (!(jentry%NPRT)) { 
	    f125_el_evt->Fill(time,gemChan,amp);     gemtrk_x2ch=(ftrk.Eval(gemtrk_x)+50.)/0.4; }
#ifdef USE_TRK
	  f125_el_amp2ds->Fill(time,gemChan,amp);
	  f125_el_fit->Fill(time,gemChan,amp);
#else
	  f125_el_amp2d->Fill(time,gemChan,amp);
#endif
	}
      } else  {   //--------------- hadron / pion ---------------
	if (gemChan>-1) {
	  if (!(jentry%NPRT)) {
	    f125_pi_evt->Fill(time,gemChan,amp);
	  }
#ifdef USE_TRK
	  f125_pi_amp2ds->Fill(time,gemChan,amp);
	  f125_pi_fit->Fill(time,gemChan,amp);
#else
	  f125_pi_amp2d->Fill(time,gemChan,amp);
#endif
	}
      }
    } //-------------------------------- End f125 pulse loop -----------------------------

    double x0=0;
    double chi2_max=20000;
#ifdef USE_TRK
    //--------------------- Calorimeter Correlation with GemTrk and GemTrd ----------------
    double chi2cc = TrkFit(f125_el_fit,fx,"fx",1);
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
      //--------------  SRS peaks  corr ---------
      for (ULong64_t i=0;i<gem_peak_count; i++) { 
	double pos = gem_peak_real_pos->at(i); 
	if (pos<=0) pos+=50.; else pos-=50.;  pos*=-1.; double pos2=ftrk.Eval(pos);
	if (gem_peak_plane_name->at(i) == "URWELLY" ) {  srs_gem_dx->Fill(x0,pos);  } 
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
      double chi2el = TrkFit(f125_el_fit, fx, "fx",0);
      Double_t p0x = fx.GetParameter(0);
      Double_t p1x = fx.GetParameter(1);
      f125_el_chi2->Fill(chi2el);
      //if ( 1 || chi2el>0. && chi2el<chi2_max) {
      if (chi2el>0. && chi2el<chi2_max) {
	f125_el_amp2d->Add(f125_el_amp2ds);   f125_el_fita->Fill(p1x); 
	if ( -0.04 < p1x && p1x < 0.02) {
	  Count("n_trk_el");
	  N_trk_el++;
	  isSingleTrack=true;
	}
      }
    } else if (pion) {
      double chi2pi = TrkFit(f125_pi_fit, fx, "fx",0);
      Double_t p0x = fx.GetParameter(0);
      Double_t p1x = fx.GetParameter(1);
      f125_pi_chi2->Fill(chi2pi);
      //if ( 1 || chi2pi>0. && chi2pi<chi2_max) {
      if (chi2pi>0. && chi2pi<chi2_max) {
	f125_pi_amp2d->Add(f125_pi_amp2ds); f125_pi_fita->Fill(p1x);
	if ( -0.04 < p1x && p1x < 0.02) {
	  Count("n_trk_pi");
	  N_trk_pi++;
	  isSingleTrack=true;
	}
      }
    }
#endif


    //----------------  GEM TRK fuducial area selection (box cut) ----------------------

    delta_x=abs(ftrk.Eval(gemtrk_x)-x0); 
    int BoxCut=1;

    if  (3200 <= RunNum && RunNum <= 3203) {  
      //xbc1=-25., xbc2=+22., ybc1=-40., ybc2=+10.;
      xbc1=0., xbc2=+22., ybc1=-20., ybc2=+10.;
      if ( gem_scluster_count!=1  ||  delta_x>dx_thresh || xbc1>gemtrk_x || gemtrk_x>xbc2 || ybc1 > gemtrk_y || gemtrk_y > ybc2)  BoxCut=0;
    }   
    if  (3248 <= RunNum && RunNum <= 3248) {  
      //xbc1=-25., xbc2=+22., ybc1=-40., ybc2=+10.;
      xbc1=10., xbc2=+26., ybc1=-5., ybc2=+15.;
      if ( gem_scluster_count!=1  ||  delta_x>dx_thresh || xbc1>gemtrk_x || gemtrk_x>xbc2 || ybc1 > gemtrk_y || gemtrk_y > ybc2)  BoxCut=0;
    }   

    //------------------- Single Track Event Histogram Filling -----

    if (isSingleTrack  && BoxCut) {                    //--- single TRD track and single SRS hit --

      Count("1_TRK");
      
      double amp_max=0;

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
	if(electron) {     //------------- electron: by Cherekov and emCal  -------
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

	    //--------------------------------------------------------------		

	    if (amp_max<amp) amp_max=amp;

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
	    if (amp_max<amp) amp_max=amp;

	    f125_pi->Fill(amp);
	    f125_pi_clu2d->Fill(time,gemChan,1.);
	            	
	    gem_xpos.push_back(gemChan);
	    gem_dedx.push_back(amp);
	    gem_zpos.push_back(time);
	    gem_parID.push_back(electron);
	    gem_nhit++;
	    gem_zHist->Fill(time, amp);
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

      if (electron) 
	f125_el_max->Fill(amp_max);
      else if (pion)
	f125_pi_max->Fill(amp_max);

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
    hevt->Reset();
    hevtc->Reset();
    hevti->Reset();
    hevtf->Reset();
	
    for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
      //  if (jentry<MAX_PRINT) printf("F125:RAW: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
      //			   ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
      //			   ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
	
      int fadc_window = f125_wraw_samples_count->at(i);
		
      int fADCSlot = f125_wraw_slot->at(i);
      int fADCChan = f125_wraw_channel->at(i);
      int gemChan = GetGEMChan(fADCChan, fADCSlot);
      int amax=0;
      int tmax=0;
      if (gemChan<0) continue;
      double DEDX_THR = 120;
      int TimeWindowStart = 95;
		
      for (int si=0; si<fadc_window; si++) {
	//printf("f125Loop:: %d fadc_window=%d\n",si,fadc_window);
	int time=si;
	int adc = f125_wraw_samples->at(f125_wraw_samples_index->at(i)+si); // printf(" sample=%d adc=%d \n",si,adc);
	if (adc>4090) printf("!!!!!!!!!!!!!!!!!!!!!! ADC 125 overflow: %d \n",adc);
	if (adc>amax) {
	  amax=adc;
	  tmax=si;
	}
	if (adc>DEDX_THR) {
	  double adc_fill=adc;
	  if (electron_ch)  {
	    f125_el_raw->Fill(time,gemChan,adc);
	  }else  {
	    f125_pi_raw->Fill(time,gemChan,adc);
	  }
	  time-=TimeWindowStart;
	  if ( 0 > time || time > 100 ) continue; // --- drop early and late hits ---

	  //hevtc->Fill(time-100,gemChan,adc/100.);
	  hevtc->SetBinContent(100-time,gemChan,adc/100.);

	  double xp = gemChan/250.*100-50.;
	  double zp = (time)/(100.)*30;
	  double ap = adc/100.;
	  //hevt->Fill(zp,xp,ap);
	  hevt->SetBinContent(100-time,gemChan,adc/100.);
	}
      } // --  end of samples loop
    } // -- end of fadc125 channels loop
  	  
      //    if (!(jentry%NPRT)) {
#if (USE_GNN==0)
      int nx = f125_el_raw->GetNbinsX();
      int ny = f125_el_raw->GetNbinsY();
      double pedestal=100.;
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
    printf("+++                           RAW TRD DATA                                                         +++\n");
    printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      for (int ii=0; ii<nx; ii++) {
	for (int jj=0; jj<ny; jj++) {
	if (electron_ch) {
	    double cc = f125_el_raw->GetBinContent(ii, jj);
	    //printf("EL: %d %d cc=%f \n",ii,jj,cc);
	  //if (cc == 0.) f125_el_raw->SetBinContent(ii,jj,pedestal);
	} else  {
	    double cc = f125_pi_raw->GetBinContent(ii, jj);
	  //printf("PI: %d %d cc=%f \n",ii,jj,cc);
	  //if (cc == 0.) f125_pi_raw->SetBinContent(ii,jj,pedestal);
	  }
	}
      }
#endif
    //=================================================================================================================
#if (USE_GNN>0) 
      // -------------------------------   hist dist clustering         ------------------------
#define MAX_CLUST 500
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

      //hevti->Reset();
      //for (int i=0; i<EVENT_SIZE; i++) {

      TH2F* hp = hevt; // -- hevt and hevtc should be same bin size
      TH2F* hpc = hevtc;

      int nx=hp->GetNbinsX();    int ny=hp->GetNbinsY(); 
      double xmi=hp->GetXaxis()->GetBinLowEdge(1);     double xma=hp->GetXaxis()->GetBinUpEdge(nx); 
      double ymi=hp->GetYaxis()->GetBinLowEdge(1);     double yma=hp->GetYaxis()->GetBinUpEdge(ny); 
      double binx = (xma-xmi)/nx;      double biny = (yma-ymi)/ny;
      printf("nx=%d,ny=%d,xmi=%f,xma=%f,ymi=%f,yma=%f\n",nx,ny,xmi,xma,ymi,yma);

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
	      //if (k==9) printf("L:1: k=%d y1=%f min=%f max=%f \n",k,y1, clust_Width[k][0], clust_Width[k][1] );
	      if (y1<clust_Width[k][0]) clust_Width[k][0]=y1; if (y1>clust_Width[k][1]) clust_Width[k][1]=y1; clust_Width[k][2]=clust_Width[k][1]-clust_Width[k][0];
	      if (x1<clust_Length[k][0]) clust_Length[k][0]=x1;if (x1>clust_Length[k][1]) clust_Length[k][1]=x1;clust_Length[k][2]=clust_Length[k][1]-clust_Length[k][0];
	      //if (k==9) printf("L:2: k=%d y1=%f min=%f max=%f \n",k,y1, clust_Width[k][0], clust_Width[k][1] );
	      // Var(X_n) = Var(X_{n-1}) + \frac{(X_n - \bar{X_{n-1}})^2}{n}
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
      } //----------------------------------- end  clustering loop -----------------------------------------------

      int MinClustSize=10;
      double MinClustWidth=0.001;
      double MinClustLength=0.01;
      double zStart =  5.; // mm
      double zEnd   = 29.; // mm
      //int ihit=0;
      int ii=0;
      printf("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
      printf("                Xpos   Ypos   Zpos       E    Width  Length   Size \n");
      //       0 Clust( 1):   43.8    0.0    3.9      5.3    0.0    0.0      4.0 
      for (int k=0; k<nclust; k++) {
	//hevt->Fill(clust_Zpos[k],clust_Xpos[k],clust_dEdx[k]);
	//hevti->Fill(hits_Zpos[k],hits_Xpos[k],dEdx[k]);
	//printf("%d Clust: X,Y,Z,E = %f %f %f %f \n",k,hits_Xpos[k],hits_Ypos[k],hits_Zpos[k],hits_dEdx[k]);
	printf("%2d Clust(%2d): %6.1f %6.1f %6.1f %8.1f %6.2f %6.2f %8.1f  ",k,k+1,clust_Xpos[k],clust_Ypos[k],clust_Zpos[k],clust_dEdx[k],clust_Width[k][2],clust_Length[k][2],clust_Size[k]);
	//printf("                                  %6.1f %6.1f %6.1f %6.1f \n",clust_Width[k][0],clust_Width[k][1], clust_Length[k][0], clust_Length[k][1]);

	//-------------  Cluster Filter -----------------

	if (clust_Size[k] >= MinClustSize && zStart < clust_Zpos[k] && clust_Zpos[k] < zEnd && clust_Width[k][2]>MinClustWidth ) {
	  hits_Xpos[ii]=clust_Xpos[k];
	  hits_Ypos[ii]=clust_Ypos[k];
	  hits_Zpos[ii]=clust_Zpos[k];
	  hits_dEdx[ii]=clust_dEdx[k];
	  ii++;
	  printf("\n");
	} else {
	  printf(" <--- skip \n");
	}
      }
      int nhits=ii;
      // -----------------------         end hist dist clustering          ----------------------------------------

      //=================================== Draw HITS and CLUST  ============================================

      char hevtTitle[80]; sprintf(hevtTitle," Display Event: %lld   Run: %d; z pos,mm; y pos,mm ",jentry,RunNum);
      hevt->SetTitle(hevtTitle);
      printf("hits_SIZE=%d  Clust size = %d \n",nhits,nclust);
      if (jentry<1000) {
	c2->cd(1);   hevt->Draw("colz");
	c2->cd(2);   hevtf->Draw("text");
	c2->cd(3);   hevti->Draw("colz");
	//c0->cd(10);   hevti->Draw("colz");
	c2->Modified(); c2->Update();
	c2->cd(1); gPad->Modified(); gPad->Update();
	int COLMAP[]={1,2,3,4,6,5};
	int pmt=22 ,pmt0 = 20; // PM type 
	for(int i = 0; i < nclust; i++){
	  //printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
	  TMarker m = TMarker(clust_Zpos[i],clust_Xpos[i],pmt); 
	  int tcol=2; //min(tracks[i],6);
	  if (clust_Size[i]<MinClustSize) pmt=22; else pmt=pmt0;
	  int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(pmt);    
	  //m.SetMarkerSize(0.5+clust_Size[i]/100); 
	  m.SetMarkerSize(0.7+clust_dEdx[i]/300); 
	  m.DrawClone();  gPad->Modified(); gPad->Update();  
	}
	c2->Modified(); c2->Update();  
      }

#if (USE_GNN==1)   // GNN MC
      //----------------------------------------------------------
      //--   Send to Model simulation 
      //----------------------------------------------------------
      
      /*
      if (jentry==0) {
	printf("++++++++++++++  test event from FPGA ++++++++++++");
	nclust=nodes_tst.size();
	for (int nt=0; nt<nclust; nt++) {
	  hits_Xpos[nt]=nodes_tst[nt][0];
	  hits_Zpos[nt]=nodes_tst[nt][1];
	  clust_Xpos[nt]=hits_Xpos[nt];
	  clust_Zpos[nt]=hits_Zpos[nt];
	}
      }
      */
      printf("**> Start Model simulation nclust=%d nhits=%d \n",nclust,nhits);
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
      
      // printVector("tracks", tracks);
      printf("**> End Model simulation \n"); //===================================================

      c2->cd(2); gPad->Modified(); gPad->Update();
      int COLMAP[]={1,2,3,4,6,5};
      for(int i = 0; i < tracks.size(); i++){
	printf("i=%d trk=%d |  %8.2f,%8.2f\n",i, tracks[i], Xcl[i], Zcl[i]);
	TMarker m = TMarker(hits_Zpos[i],hits_Xpos[i],24); 
	int tcol=min(tracks[i],6);
	int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(41);     m.SetMarkerSize(1.5); 
	m.DrawClone();  gPad->Modified(); gPad->Update();  
      }
      printf("\n\n");
      printf("**> End Cluster Plot \n");

      //--------------------------------------------------
      //----           Track fitting                 -----
      //--------------------------------------------------

      printf("==> GNN: tracks sort  : trk_siz=%ld \r\n", tracks.size());

      //-----------------   tracks sorting -------------
      // typedef std::vector<std::vector<float>> f2vec;
      //#define f2vec(X,Y) hls::vector<hls::vector<t_data,Y>,X>

      //float TRACKS[N_NODES_MAX][N_NODES_MAX*N_FEATURES];
      //std::vector<std::vector<std::vector<float>>> TRACKS;  // -- [Nnodes][N_params]
      std::vector<std::vector<float>> TRACKS;  // -- [Nnodes][N_params]
      TRACKS.resize(nhits);

      std::vector<float> hit_coord(2,0);

      //int TRACKS_N[N_NODES_MAX];
      std::vector<int>  TRACKS_N(nhits, 0); // [Nnodes]
      for (int i = 0; i < nhits; i++)  { TRACKS_N[i] = 0;  }
      std::vector<float> xz(2,0);

      for (int i2 = 0; i2 < nhits; i2++) {
	int num =  tracks[i2]; 
	int num2 = std::max(0, std::min(num, nhits - 1));
	printf("==> lstm3:track sort i=%d  : num=%d(%d) x=%f z=%f \n", i2, num, num2,  Xcl[i2],Zcl[i2]);
	xz[0]=Xcl[i2]; 	xz[1]=Zcl[i2]; 
	//TRACKS[num2].push_back(xz);
	TRACKS[num2].push_back(Xcl[i2]); TRACKS[num2].push_back(Zcl[i2]);
	TRACKS_N[num2]++;
      }
      
#if (DEBUG > 1)
      for (int i2 = 0; i2 < nhits; i2++) {
	printf(" trdID=%d n_hits=%d v_size=%d \n",i2,TRACKS_N[i2],TRACKS[i2].size());
	for (int i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
	  //printf(" i2=%d  i3=%d nhits-on-track=%d sz2=%d\n",i2,i3,TRACKS_N[i2], TRACKS[i2].size());
	  printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
	}
	if ( TRACKS_N[i2]>0) printf("\n");
      }
#endif
      //--- end tracks sorting ---
      

      //-----------------------------------
      //---       LSTM fitting          ---
      //-----------------------------------
      /*
TGraph *gr1 = new TGraph(...
TGraphErrors *gr2 = new TGraphErrors(...
TMultiGraph *mg = new TMultiGraph();
mg->Add(gr1,"lp");
mg->Add(gr2,"cp");
mg->Draw("a");
      */
      //float TRACKS_fit[N_NODES_MAX][LSTM_FEA+1];



      static TMultiGraph *mg;
      if (mg != NULL ) delete mg;
      //printf("!!!!!!!!!!!>>>>>>  mg pointer1 = %p ",mg);
      mg = new TMultiGraph();
      mg->SetTitle(" ML-FPGA response; z pos,mm; y pos,mm ");
      //printf(">>>>>>  pointer2 = %p \n ",mg);

      int NTRACKS=0;
      for (int i2 = 1; i2 < nhits; i2++) {  // tracks loop; zero track -> noise

	if (TRACKS_N[i2]<2) continue;   //---- select 2 (x,z) and more hits on track ----
	printf("==> fit: start trk: %d \r\n", i2);

	std::vector<Double_t> x;
	std::vector<Double_t> y;

	for (int i3 = 0; i3 < TRACKS[i2].size(); i3+=2) {
	  printf(" trkID=%d  hit=%d x=%f z=%f \n",i2,i3/2,TRACKS[i2].at(i3),TRACKS[i2].at(i3+1));
	  x.push_back(TRACKS[i2].at(i3+1));
	  y.push_back(TRACKS[i2].at(i3));
	}

	//c2->cd(3);   hevti->Draw("colz");

	TGraph *g = new TGraph(TRACKS_N[i2], &x[0], &y[0]);  g->SetMarkerStyle(21); g->SetMarkerColor(i2);
	/*
	  TMarker *m = new TMarker(nx[j], 0.5*ny[j], 22);
	  m->SetMarkerSize(2);
	  m->SetMarkerColor(31+j);
	  m->Draw();
	*/

	//c2->cd(3);  g->Draw("AC*");
	//TF1 *f = new TF1("f", "[2] * x * x + [1] * x + [0]"); 
	TF1 *f = new TF1("f", "[1] * x + [0]"); 
	g->Fit(f); 
	//c2->cd(3); g->Draw("AC*"); 

	//  --- get fit parameters ---
	TF1 *ffunc=g->GetFunction("f");
	Double_t p0=ffunc->GetParameter(0);
	Double_t p1=ffunc->GetParameter(1);
	printf("+++++>>  Track = %d fit: p0=%f p1=%f \n",i2,p0,p1);

	mg->Add(g,"p");

	NTRACKS++;

      }  //  end tracks loop
            
      c2->cd(3); mg->Draw("APsame"); 
      mg->GetXaxis()->SetLimits(0.,30);
      mg->SetMinimum(-50.);
      mg->SetMaximum(+50.);
      gPad->Modified(); gPad->Update(); 


#endif // USE_GNN MC


//******************************************************************************
//***
//******************************************************************************

#if (USE_TCP==1)
      //----------------------------------------------------------
      //---                 Send to FPGA                      ----
      //----------------------------------------------------------

      //-----------------  send DATA  ----------------
      //printf(" send DATA  \n"); 
      int DC_NROC=4;

      int LEN_HITS=nhits; //-- floats X,Y,Z,E 
      if (LEN_HITS>50) LEN_HITS=50; //--- max hits to fpga
      int k=3; //-- start data 3 words header;
      printf("-- BUFFER:: \n");
      for (int n=0; n<LEN_HITS; n++) {
	FBUFFER[k++]= hits_Xpos[n];
	printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
	//FBUFFER[k++]= hits_Ypos[n];
	//printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
	FBUFFER[k++]= hits_Zpos[n];
	printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
	//FBUFFER[k++]= hits_dEdx[n];
	//printf("%2d,0x%08x,f:%f ",(k-1),BUFFER[k-1],FBUFFER[k-1]);
	printf("=====>  Clust data: X,Y,Z,E=%f %f %f %f \n",hits_Xpos[n],hits_Ypos[n],hits_Zpos[n],hits_dEdx[n]);
      }

      printf("Filled bufer size=%d , data = %d hits =%d \n",k,k-3,(k-3)/2);

     itrg++;
      for (int i=0;i<nmod;i++) {
	LENEVENT=k;
	hdr=(i & 0xF) << 24;
              
	BUFFER[0]=DC_NROC; // DC_NROC if send to EVB, or  
	BUFFER[0] |= EVT_Type_DATA;
	BUFFER[1]=itrg;
	BUFFER[2]=LENEVENT;  //--
              
	BUFFER[0]&= ~(0xff << 24);             // --  clear ModID
	unsigned int MODIDx=((modID+i)&0xff)<<24;   //-- ModID=modID+i
	BUFFER[0] |= MODIDx;
              
              
	unsigned int evtTrigID=BUFFER[1];
	int evtModID=(BUFFER[0]>>24)&0xff;
	int evtSize=BUFFER[2];
              
	if ( itrg<200) {
	  printf("==> SEND:: Trg=%d(%d,%d) Mod=%d(%d) siz=%d(%d)\n"
		 ,evtTrigID,itrg,BUFFER[1],evtModID,i,evtSize,LENEVENT);
	}
              
	//rc=tcp_event_snd(BUFFER,LENEVENT,nmod,i,hdr,itrg);
	//if (rc<0) { printf(" ERROR send \n"); sleep(1); }
	//--------------------------------------------------
	//int tcp_event_snd( unsigned int *DATA, int lenDATA,int n,int k, unsigned int evtHDR, unsigned int TriggerID )

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

	printf("read GNN out, wait for FPGA data ... \n");  //=======================================================

#define MAX_NODES 100
	unsigned int NDATA[MAX_NODES+10];
	int RHEADER[10];
	int COLMAP[]={1,2,3,4,6,5};
	
	sock->RecvRaw((char*) RHEADER, sizeof(RHEADER), kDefault);
	int lenNODES=RHEADER[2];
	printf("RHEADER::"); for (int ih=0; ih<5; ih++) printf(" %d \n",HEADER[ih]);  printf(" LenDATA=%d \n",HEADER[2]);
	sock->RecvRaw((char*) NDATA, lenNODES*4, kDefault);
	
	int nnodes=lenNODES-3;
	printf("nodes return: %d (nclust=%d), TRKS: \n", nnodes,nclust);


	unsigned int TDATA[2048];
	float *FTDATA = (float*) TDATA;
#if (USE_FIT==1)
	printf("read FIT out, wait for FPGA data ... \n");  //=======================================================
	//RHEADER[10];
	
	sock->RecvRaw((char*) RHEADER, sizeof(RHEADER), kDefault);
	int lenFITS=RHEADER[2];
	printf("RHEADER::"); for (int ih=0; ih<5; ih++) printf(" %d \n",RHEADER[ih]);  printf(" LenDATA=%d \n",RHEADER[2]);
	sock->RecvRaw((char*) TDATA, lenFITS*4, kDefault);
	
	for (int i=0; i<lenFITS; i++) {

	  printf("tracks fit return: i=%d  data=0x%x (%f)  \n", i,TDATA[i],FTDATA[i]);

	}


	int nfits=lenFITS-3;
	printf("tracks fit return: %d  \n", nfits);

#else 
	int nfits=nnodes;
	printf("tracks fit return: %d  \n", nfits);
#endif  // --- end  if USE_FIT  ---

	if (nfits>0) { 

	  //=============== Draw FPGA Clust =============

	  for (int nd=0; nd<min(nnodes,nhits); nd++) {	  
	    int trknum=NDATA[nd+3];
	    printf(" %u, ",trknum); 
	    c2->cd(3);	gPad->Modified(); gPad->Update(); 
 	    if (trknum>0) {
	      //DrawPolyMarker (1, &hits_Xpos[nd], &hits_Zpos[n], Option_t *option="")
	      //TMarker* m = new TMarker(hits_Xpos[nd],hits_Zpos[n],24);  // memory leak !!!!!
	      TMarker m = TMarker(hits_Zpos[nd],hits_Xpos[nd],24); 
	      int tcol=min(trknum,6);
	      int mcol = COLMAP[tcol-1];   m.SetMarkerColor(mcol);   m.SetMarkerStyle(20);     m.SetMarkerSize(1.5); 
	      m.DrawClone();  gPad->Modified(); gPad->Update();  
	      printf("=====>  Draw:%d X,Y,Z,E=%f %f %f %f \n",nd,hits_Xpos[nd],hits_Ypos[nd],hits_Zpos[nd],hits_dEdx[nd]);
	    }
	  }
	  c2->Modified(); c2->Update(); 
	  printf("\n");

#if (USE_FIT==1)
	  //=============== Draw Tracks lines =============

	  int ntracks = nfits/3;
	  int cs = nfits % 3;
	  if (cs != 0)  { printf("==========>>>>   Error FIT results : %d %d %d \n",nfits, ntracks,cs);  break; };
	
	  int cnt = 3; // word counter in data buffer
	  c2->cd(3);  gPad->Modified(); gPad->Update(); 
	  for (int i=0; i<ntracks; i++) {
	  
	    int trknum = TDATA[cnt++];	  float aa = FTDATA[cnt++];	  float bb = FTDATA[cnt++];
	    printf(" Fit Track=%d aa=%f bb=%f \n",trknum,aa,bb);

	    TF1 ftrk("ftrk","[0]*x+[1]",zStart,zEnd);	  ftrk.SetParameter(0,aa);	  ftrk.SetParameter(1,bb);
	     ftrk.DrawClone("same");  gPad->Modified(); gPad->Update();  
	  }
	  c2->Modified(); c2->Update(); 
#endif  // --- end  if USE_FIT  ---

	} else {    //  if (nfits>0)
	  printf(" No tracks to draw \n");
	}

	//=============== Draw All Clust ================
	//----------------------  To FPGA ------------------
	c2->cd(2); gPad->Modified(); gPad->Update(); 
	printf(" Draw clusters  \n");  
	for (int k=0; k<nhits; k++) {
	  TMarker mh = TMarker(hits_Zpos[k],hits_Xpos[k],24); 
	  int mhcol = 1;   mh.SetMarkerColor(mhcol);  mh.SetMarkerSize(1.5); 
	  mh.DrawClone();  gPad->Modified(); gPad->Update(); 
	}
	c2->Modified(); c2->Update(); 
	//----------------------  from FPGA ------------------
	c2->cd(3); gPad->Modified(); gPad->Update(); 
	printf(" Draw clusters  \n");  
	for (int k=0; k<nhits; k++) {
	  TMarker mh = TMarker(hits_Zpos[k],hits_Xpos[k],24); 
	  int mhcol = 1;   mh.SetMarkerColor(mhcol);  mh.SetMarkerSize(1.5); 
	  //c2->cd(2); mh.DrawClone();  c2->Modified(); c2->Update(); 
	  mh.DrawClone();  gPad->Modified(); gPad->Update(); 
	}
	c2->Modified(); c2->Update(); 
      } // -- end nmod  ---
#endif  // --- end  if USE_TCP  -------------------------------------------------------------------------------------------

      printf(" all done, click low right pad ...  \n");
      //c2->cd(2); gPad->WaitPrimitive();
      //if (jentry<35) sleep(5); else sleep(1);
      
      
      //-------------------------------------------------------
 
#endif   // (USE_GNN>0) 

    //=======================  End Fa125 RAW  process Loop  =====================================================
	
    if(electron){
      //f125_el_max->Fill(f125_amp_max);
      //if((slot<6)||(slot==7&&chan<24))f125_el->Fill(f125_amp_max);
      //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_el->Fill(f125_amp_max);
    } else if (pion) {
      //f125_pi_max->Fill(f125_amp_max);
      //if((slot<6)||(slot==7&&chan<24))f125_pi->Fill(f125_amp_max);
      //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_pi->Fill(f125_amp_max);
    }

#endif   //  USE 125 RAW
	
    //=====================================================================================
    //===                     Event Display                                            ====
    //=====================================================================================
    
    if (jentry<MAX_PRINT || !(jentry%NPRT)) {
      c0->cd(1); f125_el_amp2d->Draw("colz");
      c0->cd(5); f125_pi_amp2d->Draw("colz");
      c0->cd(2); f125_el_evt->Draw("colz");         double chi2e = TrkFit(f125_el_evt,fx1,"fx1",0); if (chi2e>0. && chi2e < chi2_max ) fx1.Draw("same"); 
      c0->cd(6); f125_pi_evt->Draw("colz");         double chi2p = TrkFit(f125_pi_evt,fx2,"fx2",0); if (chi2p>0. && chi2p < chi2_max ) fx2.Draw("same"); 
      printf("========================>>>  Chi2 e=%f p=%f \n",chi2e,chi2p);
      //c0->cd(3); f125_el_chi2->Draw("colz");
      //c0->cd(6); f125_pi_chi2->Draw("colz");
      //
      
      //if (electron) { 
	c0->cd(3); f125_el_raw->Draw("colz");  f125_el_evt->Draw("same"); 
	TLine lin1(110.,gemtrk_x2ch,190.,gemtrk_x2ch); lin1.Draw("same");   //--- draw  gemtrk x 
	printf("++++++++++++ Draw GEMTRK:: %f %f %f  \n",gemtrk_x,ftrk.Eval(gemtrk_x),gemtrk_x2ch);
	//}
	//if (pion) { 
	c0->cd(7); f125_pi_raw->Draw("colz");  f125_pi_evt->Draw("same"); 
	TLine lin2(110.,gemtrk_x2ch,190.,gemtrk_x2ch); lin2.SetLineColor(kRed); lin2.Draw("same");   //--- draw  gemtrk x 
	//c0->cd(12); gPad->WaitPrimitive();
	//}
      if (electron || pion ) {
	//printf("--------------->  SRS:: srs_peak:  cnt=%lld evt=%llu \n",gem_peak_count,jentry);
	int lc=0; 
	for (int k=0; k<100; k++) {  peak_line[k].SetX1 (100.);	    peak_line[lc].SetY1 (-10.);	    peak_line[lc].SetX2 (101.);	    peak_line[lc].SetY2 (-10.); }
	for (ULong64_t i=0;i<gem_peak_count; i++) { 
	  //printf("SRS:: srs_peak: i=%lld id=%d name=%s idx=%d apv=%d Amp=%f wid=%f E=%f Pos=%f \n"
	  //	 ,i,gem_peak_plane_id->at(i),gem_peak_plane_name->at(i).c_str(),gem_peak_index->at(i), gem_peak_apv_id->at(i), gem_peak_height->at(i)
	  //	 ,gem_peak_width->at(i), gem_peak_area->at(i), gem_peak_real_pos->at(i));
	  double pos = gem_peak_real_pos->at(i); 
	  if (pos<=0) pos=pos+50.; else pos=pos-50.;  pos*=-1.;
	  double pos2ch=(ftrk.Eval(pos)+50.)/0.4;  // -- to gemtrd coordinate system
	  peak_line[lc].SetX1 (110.);	    peak_line[lc].SetY1 (pos2ch);	    peak_line[lc].SetX2 (190.);	    peak_line[lc].SetY2 (pos2ch);
	  //printf(" i=%llu pos=%f pos2ch=%f \n ",i,pos,pos2ch);
	  if (gem_peak_plane_name->at(i) == "URWELLX" ) { if (lc<100) {  peak_line[lc].SetLineColor(kGreen);  peak_line[lc].Draw("same");  lc++;   }   } 
	}  //--- peak Loop --
      }
      c0->cd(4); cal_el_evt->Draw("colz");
      c0->cd(8); cal_pi_evt->Draw("colz");

      //c0->cd(9); srs_trk_el->Draw("colz"); 
      //c0->cd(10); srs_gem_x->Draw("colz");  ftrk.Draw("same"); 
      //c0->cd(9); hevt->Draw("colz");
      //c0->cd(10); hevtc->Draw("colz");
      c0->Modified();   c0->Update();     

      //---------- fiducial area ---

      // srs_etrd_corr->Divide(srs_etrd_beam);  
      srs_etrd_ratio = (TH2F*)srs_etrd_corr->Clone("srs_etrd_ratio");  
      //srs_etrd_corr->Copy((TH2F*)srs_etrd_ratio);  
      srs_etrd_ratio->GetXaxis()->SetTitle(" "); srs_etrd_ratio->GetYaxis()->SetTitle(" "); 
      srs_etrd_ratio->SetTitle("TR energy norm");  //srs_etrd_ratio->SetStats(1); srs_etrd_ratio->SetMinimum(0.8); srs_etrd_ratio->SetMaximum(1.35); 
#if 0 
      srs_etrd_ratio->Divide(srs_etrd_beam);
#esle
      srs_etrd_ratio->Add(srs_etrd_beam,-1.);
#endif
      /*
      c0->cd(11); srs_etrd_ratio->DrawCopy("colz");
      c0->Modified();   c0->Update();     
      //c0->cd(11); srs_etrd_corr->Draw("colz");
      TBox fbox(xbc1,ybc1,xbc2,ybc2);  //---- draw box cut ---
      fbox.SetLineColor(kRed);
      fbox.SetFillStyle(0);
      fbox.SetLineWidth(2);
      fbox.DrawClone();
      */
     //---------
      c0->cd(12); srs_gem_dx->Draw("colz");
      c0->Modified();   c0->Update();
      //if (NPRT<1000) sleep(1);
      c0->cd(12); gPad->WaitPrimitive();

    }
    
	
    //--- Fill Track Hit Tree ---
    //printf("Filling hits tree ... ev=%d  Nhits=%d  Ntracks=%d \n", event_num, nhit, track_size);
    //printf("Filling hits tree ... ev=%d  Nhits=%d \n", event_num, nhit);
    if (gem_nhit>0) EVENT_VECT_GEM->Fill();
    if (mmg1_nhit>0)EVENT_VECT_MMG1->Fill();
    if (mmg2_nhit>0)EVENT_VECT_MMG2->Fill();
    if (urw_nhit>0)EVENT_VECT_URW->Fill();
    
    } // ------------------------ end of event loop  ------------------------------
   
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
  char pngname[120];  sprintf(pngname,"%s_evdisp.png",G_DIR);  c0->Print(pngname);
  char pdfname[120];  sprintf(pdfname,"%s_evdisp.pdf",G_DIR);  c0->Print(pdfname);

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
  cc=NextPlot(nxd,nyd);  srs_gem_dx->Draw("colz");  
  cc=NextPlot(nxd,nyd);  srs_gem_x->Draw("colz");  ftrk.Draw("same"); 
  cc=NextPlot(nxd,nyd);  srs_gem_y->Draw("colz");
  cc=NextPlot(nxd,nyd);  srs_etrd_corr->Draw("colz");  
  TBox fbox(xbc1,ybc1,xbc2,ybc2);  //---- draw box cut ---
  fbox.Draw("same");
  //fbox.SetFillColorAlpha(0,0);
  fbox.SetLineColor(kRed);
  fbox.SetFillStyle(0);
  fbox.SetLineWidth(1);
  cc=NextPlot(nxd,nyd);  srs_etrd_beam->Draw("colz");  
  //cc=NextPlot(nxd,nyd);  srs_etrd_ratio->Draw("colz");  

 //---------------------  page 3 --------------------
  htitle("  GEMTRD (fa125) Amp ");    if (!COMPACT) cc=NextPlot(0,0);
  nxd=2; nyd=4;
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   f125_el->Draw();
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   f125_pi->Draw();
  cc=NextPlot(nxd,nyd);     f125_el_max->Draw();
  cc=NextPlot(nxd,nyd);     f125_pi_max->Draw();
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   mmg1_f125_el->Draw();
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   mmg1_f125_pi->Draw();
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   urw_f125_el->Draw();
  cc=NextPlot(nxd,nyd);   gPad->SetLogy();   urw_f125_pi->Draw();

 //---------------------  page 3a --------------------
  htitle("  GEMTRD (fa125) Amp 2D");    if (!COMPACT) cc=NextPlot(0,0);
  nxd=2; nyd=3;
  cc=NextPlot(nxd,nyd);   f125_el_amp2d->Draw("box");
  cc=NextPlot(nxd,nyd);   f125_pi_amp2d->Draw("box");
  cc=NextPlot(nxd,nyd);   mmg1_f125_el_amp2d->Draw("box");
  cc=NextPlot(nxd,nyd);   mmg1_f125_pi_amp2d->Draw("box");
  cc=NextPlot(nxd,nyd);   urw_f125_el_amp2d->Draw("box");
  cc=NextPlot(nxd,nyd);   urw_f125_pi_amp2d->Draw("box");

  //---------------------  page 4 --------------------
  htitle(" Clust   ");    if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);   f125_el_clu2d->Draw("box");
  cc=NextPlot(nxd,nyd);   f125_pi_clu2d->Draw("box");
  cc=NextPlot(nxd,nyd);   mmg1_f125_el_clu2d->Draw("box");
  cc=NextPlot(nxd,nyd);   mmg1_f125_pi_clu2d->Draw("box");
  cc=NextPlot(nxd,nyd);   urw_f125_el_clu2d->Draw("box");
  cc=NextPlot(nxd,nyd);   urw_f125_pi_clu2d->Draw("box");

  //---------------------  page 5 --------------------
  htitle(" Tracking   ");    if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);  f125_el_chi2->Draw("colz");
  cc=NextPlot(nxd,nyd);  f125_pi_chi2->Draw("colz");
  cc=NextPlot(nxd,nyd);  f125_el_fita->Draw("colz");
  cc=NextPlot(nxd,nyd);  f125_pi_fita->Draw("colz");

  //--- close PDF file ----
  cc=NextPlot(-1,-1);
  //--- the end ---
  
  }
//===============================================================
