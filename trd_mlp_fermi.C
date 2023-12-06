//
// .x   trd_mlp.C++("hd_rawdata_000554_000.evio.root",1);   //   1-WC,  0-GEM
// .x   trd_mlp.C++(554,1);                                 // new
// .x   trd_mlp_fermi.C++(3200);    // fermi
#include <sstream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"
#include <iostream>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TString.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLinearFitter.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TPaveStats.h"
#include "TMultiGraph.h"
#include "TSystem.h"
#include "stdio.h"
#include "trd_mlp_fermi.h"
#include "PlotLib.C"

void Count(const char *tit);
void Count(const char *tit, double cut1);
void Count(const char *tit, double cut1, double cut2);
TH1D *hcount;

//--  0=dEdx 1=Params 2=dEdx+Params ; change MAXpar to 17 !
#define NN_MODE 3

const int NDEslices = 10;
const int NFixed = 7;

#if    NN_MODE == 0
const int MAXpar = NDEslices; //10; 17
#elif  NN_MODE == 1
const int MAXpar = NFixed; //10; 17
#elif  NN_MODE == 2
const int MAXpar = NFixed+NDEslices; //10; 17
#elif  NN_MODE == 3
const int MAXpar = NFixed+NDEslices; //10; 17
#else
  ne rabotaet !!!!
#endif

Float_t Par[MAXpar];
Int_t type, ievent, channel;
TCanvas *c2 = NULL;
TH2F  *disp, *dispe, *disppi;
int RunNum;

//----------------------------------------------------------------------
int hgive1(TH1 *hp, int *nx, double *xmi, double *xma) {
  //int hgive1(TH1F *hp) {
  //Int_t  nx=0;
  //double xmi=0., xma=0., sum=0.;
  
  *nx=hp->GetNbinsX();
  *xmi=hp->GetBinLowEdge(1);  //h1.GetXaxis().GetBinLowEdge(bin)
  *xma=hp->GetBinLowEdge(*nx)+hp->GetBinWidth(*nx); // h1.GetXaxis().GetBinWidth(bin)
  //sum=hp->Integral(1,nx);

  //cout << " hist: nx=" << nx << " xmi=" << xmi << " xma=" << xma << "  sum=" << sum  << endl;
  
  //hp->GetBinContent(int i);
  //  GetEntries()
  //  Integral(Int_t    binx1,  Int_t   binx2,  Option_t *      option = "" );
  //Double_t TH1::Interpolate ( Double_t  x )
  return hp->GetEntries();
}

//-------------------------------------------------------------------------------
int hscale(TH1 *he, TH1 *hpi, double scale, int NORM, int DRAW) {
  int ret = 0;
  int noent_e = he->GetEntries();
  int noent_pi = hpi->GetEntries();
  double escale = 1;
  if (scale>0) escale = scale;
  else if (noent_e > 0) escale = (double)noent_pi/(double)noent_e;
  else ret = 1; // -- 1 == error

  string name = hpi->GetName();
  printf("NORM :: hist=%s<<< noent e=%d pi=%d e-scale=%f \n", name.c_str(), noent_e, noent_pi, escale);
  if (NORM) he->Scale(escale);
  if (DRAW>0) {
    double maxe = he->GetMaximum();
    double maxpi = hpi->GetMaximum();
    if (maxe>maxpi) {
      he->Draw("hist");  hpi->Draw("histsames");
    } else {
      hpi->Draw("hist");  he->Draw("histsames");
    }
    he->SetLineColor(2);  hpi->SetLineColor(4);
  }
  if (DRAW>1) {
     //--- 2 stat boxes
     gPad->Update();
     TH1 *h1 = he;  TH1 *h0 = hpi;
     TPaveStats *ps1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
     ps1->SetY1NDC(0.57);  ps1->SetY2NDC(0.75); ps1->SetTextColor(kRed);
     TPaveStats *ps0 = (TPaveStats*)h0->GetListOfFunctions()->FindObject("stats");
     ps0->SetTextColor(kBlue);
     gPad->Modified(); gPad->Update();
     //------
  }
  return ret;
} //-- end hscale

//------- Rejection Factor Calculation -----------------------------------
double Reject(TH1 *hp, TH1 *he, double thr) {
  double s[202][2], r1, w, Rej, XMI, XMA, es, ps, e1=0. , e2=0., p1=0., p2=0. ;
  int NX, i0, noente=0, noentp=0;
  Rej=-1; r1=1-thr; i0=thr; w=thr-i0;

  //............ electron 10% of integral level ..........................
  noente = hgive1(he, &NX, &XMI, &XMA);
  TString etitle = he->GetTitle();
  cout << "===> " << etitle << " <===  NX e = " << NX <<  " w = " << w <<  " i0 = " << i0 << endl;
  if (NX>200) { cout << "error NX e = " << NX << endl; exit(1); }
  if (noente>0) {
    es=0;  for ( int i=0; i<=NX+1; i++)  {  es=es+he->GetBinContent(i); s[i][1]=es; };
    e2=0;
    cout << " es = " << es << " noente= " << noente << endl;
    if (es>0) {
      for (int i=0; i<=NX+1; i++) {
        e1=e2; e2=s[i][1]/es;
        //printf(" i=%2d  x=%5.2f cont=%5.1f sum=%6.2f e2=%4.2f\n",i, he->GetBinLowEdge(i), he->GetBinContent(i), s[i][1], e2);
        if(e2>r1) { i0=i-1; break; }
      }
      w=(r1-e1)/(e2-e1);  Rej = i0+w;
    }
  }

  //............ now pions level at this threshold .......................
  noentp = hgive1(hp, &NX, &XMI, &XMA);
  cout << "NX pi = " << NX <<  " i0=" << i0 << " thr=" << he->GetBinLowEdge(i0) << endl;
  if (noentp>0) { ps=0; for ( int i=0; i<=NX+1; i++) { ps=ps+hp->GetBinContent(i); s[i][0]=ps; };
    p1=s[i0][0]/ps;  p2=s[i0+1][0]/ps;  Rej = 1-(p1+(p2-p1)*w);
  }
  cout << "===> " << etitle << " <=== Rejection =" << 1./Rej << endl;
  return Rej;
}

//---------------------------------------------------------------------
int fill_trees( TTree *gem_hits, TTree *signal, TTree *background, TTree *sig_tst, TTree *bg_tst, int WC, int runnum, int *rtw1, int *rtw3) {

  // Set branch addresses.
  gem_hits->SetBranchAddress("event_num", &event_num, &b_event_num);
  gem_hits->SetBranchAddress("nhit", &gem_nhit, &b_gem_nhit);
  gem_hits->SetBranchAddress("xpos", &xpos, &b_xpos);
  gem_hits->SetBranchAddress("zpos", &zpos, &b_zpos);
  gem_hits->SetBranchAddress("dedx", &dedx, &b_dedx);
  gem_hits->SetBranchAddress("parID", &parID, &b_parID);
  gem_hits->SetBranchAddress("zHist", &zHist, &b_zHist);

  //========================================
  TH2F *hits2d_e = new TH2F("hits2d_e","Electron dEdx in Time ; Time (8ns) ; Channel",50,0.,300., 256,-0.5,255.5);
  TH2F *hits2d_p = new TH2F("hits2d_p","Pion dEdx in Time ; Time (8ns) ; Channel",50,0.,300., 256,-0.5,255.5);

  TH2F *aver2d_e = new TH2F("aver2d_e","aver-rms ",100,0.,240., 100,0.,100.);
  TH2F *aver2d_p = new TH2F("aver2d_p","aver-rms",100,0.,240., 100,0.,100.);

  //TH1F *hbeamX  = new TH1F("hbeamX","hbeamX",256,-0.5,255.5);
  TH1F *hNhits  = new TH1F("hNhits","Hits in GEMTRD",100,-0.5,99.5);

  //TH2F *rad2d   = new TH2F("rad2d","rad2d",50,0.,500., 50,0.,256.);
  //TH2F *norad2d = new TH2F("norad2d","norad2d",50,0.,500.,50,0.,256.);
  //TH2F *gem2d   = new TH2F("gem2d","gem2d",100,0.,256.,100,0.,256.);
  TH1F *e_amax    = new TH1F("e_amax","Max Amp per Hit ; ADC Amplitude",200,0.,4100.);
  TH1F *pi_amax  = new TH1F("pi_amax","Max Amp per Hit ; ADC Amplitude",200,0.,4100.);
  TH1F *e_etot   = new TH1F("e_etot","dEdx Total per Hit ; ADC Amplitude Sum",600,0.5,100000.);
  TH1F *pi_etot = new TH1F("pi_etot","dEdx Total per Hit ; ADC Amplitude Sum",600,0.5,100000.);
  
  //TH2F *h2xdiff = new TH2F("h2xdiff","h2xdiff; GemTRD; Ext Track",380,-124.5,255.5,380,-124.5,255.5);

  //TH1F *zrad2   = new TH1F("zrad2","zrad2",100,-0.5,99.5);
  //TH1F *znorad2 = new TH1F("znorad2","znorad2",100,-0.5,99.5);

  TH1F *time_e  = new TH1F("time_e","dEdx in Time ; Time (8ns)",330,0.5,330.5);
  TH1F *time_pi = new TH1F("time_pi","dEdx in Time ; Time (8ns)",330,0.5,330.5);
  //TH1F *ampl    = new TH1F("ampl","",500,0.5,500.5);

  TH1F *par_e[MAXpar];
  TH1F *par_pi[MAXpar];
 
  for (int ip=0; ip<MAXpar; ip++) {
    char hname[80];
    sprintf(hname, "par_e_%d ", ip);
    par_e[ip] = new TH1F(hname, hname, 100, -0.5, 99.5);
    sprintf(hname, "par_pi_%d ", ip);
    par_pi[ip] = new TH1F(hname, hname, 100, -0.5, 99.5);
  }

  //========================================
  TRandom3 *rndm = new TRandom3();
  int nx;
  double xmi,xma;
  int noent = hgive1(par_e[0], &nx, &xmi, &xma );
  cout << " hist: nx=" << nx << " xmi=" << xmi << " xma=" << xma << " noent= " << noent << endl;

  Long64_t nentries = gem_hits->GetEntries();
  Long64_t nbytes = 0;
  //int ntest=nentries*0.25;
  //cout << " ntest=" << ntest << endl;
  //==============================================================================
  //
  //==============================================================================
  const int NDE=NDEslices; // MAXpar;
  double dEdx[NDE];
  int NPF=MAXpar; // number of parameters filled
  int ntrk_e=0, ntrk_pi=0;
  int e_chan1=0;    //-- first TR channlel
  int e_chan2=0;    //-- last  TR channel
  int pi_chan1=0;   //-- first pion (no-rad) channel
  int pi_chan2=0;   //-- last  pion (no-rad) channel

  int ievOK=0;
  //--------------------------------------------------------------------------------
  //                   E V E N T    L O O P
  //--------------------------------------------------------------------------------
  //int Nfill=0;
  for (Long64_t iev = 0; iev < nentries; iev++) {
    nbytes += gem_hits->GetEntry(iev);
    ievent=iev;
    Count("EVT");
    if (!(iev % 1000)) cout << " event = "  << iev << " of " << nentries << endl;

    //double gemx=0, gemy=0;
    //gem2d->Fill(gemx, gemy);

    for (int ip=0; ip<MAXpar; ip++) Par[ip]=0;
    for (int i=0; i<NDE; i++) dEdx[i]=0;

    //------------------------------------------------------------------------------
    //                 W C
    //------------------------------------------------------------------------------
    int w2ch=22;
    int wch=22;
    int gemch=w2ch;
    int wcch=wch;
    //int trkch=dch;
    int trkch=wch;
    float zero_bin=0;
       
    //------------------------------------------------------------------------------
    //                 G E M
    //------------------------------------------------------------------------------
       
    channel = w2ch;
/*
    int USE_TRACK = 0;
    if ( USE_TRACK ) {
      if (abs(13-channel/27-trkch)<3.) { // tracking cut
        h2xdiff->Fill(channel,trkch);
        hbeamX->Fill(channel);
      }
    } else {
      h2xdiff->Fill(channel,trkch);
      hbeamX->Fill(channel);
    }
*/
    float amax2=0.;
    float emax2=0.;
    float tmax2=50.;
       
    //int nw2clust=gem_nhit; // w2nhit;

    int khit=0;
    int NTR=0;
    int tw1=110;
    int tw2=160;
    int tw3=185;

    float THR1=100.;
    float THR2=500.;
    float Escale=400.;
    float Ascale=40.;
       
    float atot=0.;
    float etot=0.;
    float etrzon=0.;

    pi_chan1=40;  //-- first pion (no-rad) channlel
    pi_chan2=220;  //-- last  pion (no-rad) channel
    //
    e_chan1=43;  //-- first TR channlel
    e_chan2=223;  //-- last  TR channel

    //---- run specific ---
    if ( 3000 < runnum && runnum <= 3999 ) { //------------------  FERMILAB   2023  ---------------------

      THR1=100.;
      THR2=350.;
      Escale=100.;
      Ascale=10.;

      e_chan1=40;          //-- first TR channlel
      e_chan2=220;         //-- last  TR channel
      pi_chan1=e_chan1+1;  //-- first pion (no-rad) channlel
      pi_chan2=e_chan2+1;  //-- last  pion (no-rad) channel

      //--- need hit from GEM TRK  and remove noisy channels!!!!
      //if (abs(channel-trkch)>85 || abs(channel-trkch)<55 || trkch == 223 || channel == 169 || channel == 170 ) continue;
      //if ( (trkch<(channel+50)) ||  channel == 49 || channel == 54 || channel == 56  ) continue;
      //if ( USE_TRACK ) { if ( abs(13-channel/27-trkch)<3 ||  channel == 49 || channel == 54 || channel == 56  ) continue; }

      tw1=110;
      tw2=160;
      tw3=185;

      switch (runnum) {

      case 3125:   tw1=110; tw2=160; tw3=185; e_chan1=80;   e_chan2=140;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Fleece (lrg)
      case 3126:   tw1=110; tw2=160; tw3=185; e_chan1=80;   e_chan2=140;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Fleece (lrg)
      case 3131:   tw1=110; tw2=160; tw3=185; e_chan1=80;   e_chan2=140;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Fleece (lrg)
      
      case 3132:   tw1=110; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Foil (VU)
      case 3133:   tw1=110; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Foil (VU)
      case 3134:   tw1=110; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Foil (VU)
      
      case 3196:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Foil
      case 3197:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Foil
      case 3198:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Foil
      case 3199:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Foil
      case 3200:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Foil
      
      case 3201:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=160;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Fleece
      case 3202:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=160;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Double Fleece
      
      case 3203:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- No Rad
      case 3204:   tw1=108; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- No Rad
      
      case 3287:   tw1=110; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Foil (TU)
      case 3288:   tw1=110; tw2=160; tw3=185; e_chan1=85;   e_chan2=165;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break; //-- Single Foil (TU)
      
      case 3216:   tw1=110; tw2=160; tw3=185; e_chan1=100;  e_chan2=185;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break;
      case 3248:   tw1=110; tw2=160; tw3=185; e_chan1=112;  e_chan2=134;  pi_chan1=e_chan1;   pi_chan2=e_chan2;   break;

      default:
      tw1=110;
      tw2=160;
      tw3=185;
    }

  } else {  //-- Runnum ---
    printf(" run %d not found in the GEM range \n",runnum);
    exit(1);
  }
  
  double dt=(tw3-tw1)/NDE;
  
  type=-1;
  if(parID->at(0)){
    type=1; ntrk_e++;
  } else {
    type=0; ntrk_pi++;
  }
  hNhits->Fill(gem_nhit);
  
  //=================== Loop to  calculate average and RMS  ========
  double xaver=0, xaver2=0;
  int naver=0;
  for (int i=0;i<gem_nhit;i++) {
    if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
    xaver+=xpos->at(i); naver++;
  }
  xaver=xaver/naver;
  
  //====================================================================
  for (int i=0;i<gem_nhit;i++) { // count fixed parameters
    
    Count("Hits");
    if (tw1 > zpos->at(i) || zpos->at(i) > tw3) continue;
    Count("zHits");
    xaver2+=((xpos->at(i)-xaver)*(xpos->at(i)-xaver));

    if (dedx->at(i)>THR1) {
      if (type==1) hits2d_e->Fill(zpos->at(i),xpos->at(i),dedx->at(i));   //--- ampl
      else if (type==0) hits2d_p->Fill(zpos->at(i),xpos->at(i),dedx->at(i));   //--- ampl
      //hits2d->Fill(zpos->at(i),gemch,w2ahit[i]); //--- energy
    }

    if(dedx->at(i)>amax2 && tw1<zpos->at(i) && zpos->at(i)<tw3) {
      amax2=dedx->at(i);
      tmax2=zpos->at(i);
    }
    if(dedx->at(i)>emax2) {  // currently same as amp
      emax2=dedx->at(i);
    }

    if (type>=0 ) {  //-- rad.pos. window OK
      if (tw1<zpos->at(i) && zpos->at(i)<tw3 && dedx->at(i)>THR1) {
        khit++;
        etot+=dedx->at(i);
        atot+=dedx->at(i);
        int ibin=(zpos->at(i)-tw1)/dt;  ibin=min(max(0,ibin),(NDE-1)); dEdx[ibin]+=dedx->at(i)/10; // w2ahit[i]/10.; w2mhit[i];
      }
      if (tw2 < zpos->at(i) && zpos->at(i) < tw3)  {  etrzon+=dedx->at(i); }
      if (dedx->at(i)>THR2) NTR++;
    }

    if (dedx->at(i)>THR1)  {
      if (type==1) time_e->Fill(zpos->at(i),dedx->at(i));
      if (type==0) time_pi->Fill(zpos->at(i),dedx->at(i));
    }
  } //--- Loop over gem hits

  //=================================================================

    xaver2=sqrt(xaver2/naver);
    if (type==1 ) {
      Count("el");
      aver2d_e->Fill(xaver,xaver2);
    } else {
      Count("pi");
      aver2d_p->Fill(xaver,xaver2);
    }
    //printf("++++>>>  e_chan1=%d xaver=%f  e_chan2=%d  \n", e_chan1,  xaver, e_chan2 );
    if ( e_chan1 > xaver || xaver > e_chan2 ) continue; //--- radiator  area  ; for Y - need a track
    Count("eChan");
    if (naver<3) continue; // --- too small number of hits ----
    Count("nHits");


    //--------------------------------------------------------------------------------
    //                    electron case
    //--------------------------------------------------------------------------------
    if(type==1) { //-- is electron
      e_amax->Fill(amax2);
      e_etot->Fill(etot);
    }
    
    //--------------------------------------------------------------------------------
    //                   pion case
    //--------------------------------------------------------------------------------
    if(type==0) { //-- is pion
      pi_amax->Fill(amax2);
      pi_etot->Fill(etot);
    }
    
    //-----------------------------------------------
    if ( type<0 ) continue;
    Count("type");
    ievOK++;

    if (NN_MODE==1) {  //-- Params only
      if (MAXpar<NFixed) { printf("ERROR :: MAXpar array too small =%d \n",MAXpar); exit(1); }
      Par[0]=amax2/Ascale;
      Par[1]=khit;
      Par[2]=emax2/Escale;
      Par[3]=etot/Escale/10.;
      Par[4]=etrzon/Escale/10.;
      Par[5]=NTR;
      Par[6]=atot/1000.;
    } else if (NN_MODE==2) {  //--  dEdx + Par
      if (MAXpar<(NFixed+NDE)) { printf("ERROR :: MAXpar array too small =%d \n",MAXpar); exit(1); }
      Par[0]=amax2/Ascale/5.;
      Par[1]=khit;
      Par[2]=emax2/Escale;
      Par[3]=etot/Escale/10.;
      Par[4]=etrzon/Escale/10.;
      Par[5]=NTR;
      Par[6]=atot/1000.;
      int np=NDE;
      double coef=Ascale*3.;  // Escale; Ascale;
      for (int ip=0; ip<np; ip++) {
        Par[ip+NFixed]=dEdx[ip]/coef;
        if ( dEdx[ip]<0.1 ) zero_bin++;
      }
    } else if (NN_MODE==3) {  //--  fermi dEdx(amp) + Par
      if (MAXpar<(NFixed+NDE)) { printf("ERROR :: MAXpar array too small =%d \n",MAXpar); exit(1); }
      Par[0]=amax2/Ascale/5.;
      Par[1]=khit;  // -- hits > THR1 = 100
      Par[2]=xaver2;
      Par[3]=etot/Escale/10.;
      Par[4]=etrzon/Escale/10.;
      Par[5]=NTR;
      Par[6]=atot/1000.;
      int np=NDE;
      double coef=Ascale*3.;  // Escale; Ascale;
      for (int ip=0; ip<np; ip++) {
        Par[ip+NFixed]=dEdx[ip]/coef;
        if ( dEdx[ip]<0.1 ) zero_bin++;
      }
    } else {   //-- dEdx only
      int np=min(MAXpar,NDE);
      double coef=Ascale*3.;  // Escale; Ascale;
      for (int ip=0; ip<np; ip++) {
        Par[ip]=dEdx[ip]/coef;
        if ( dEdx[ip]<0.1 ) zero_bin++;
      }
    }
    /*
      cout << "0: " << ievOK << " par0= " << amax2/Ascale << " par1= " << khit << " par2= " << emax2/Escale << " par3= " << etot/Escale/10. << endl;
      cout << "1: " << ievOK << " par0= " << Par[0]<< " par1= " << Par[1] << " par2= " << Par[2] << " par3= " << Par[3] << endl;
    */
    /*
      par1=rndm->Gaus(5.,2);
      par2=rndm->Gaus(10.,2);
      par3=rndm->Gaus(15.,2);
      signal->Fill();
      par1=rndm->Gaus(0.,2);
      par2=rndm->Gaus(20.,2);
      par3=rndm->Gaus(25.,2);
      background->Fill();
    */
       
    //return kTRUE;
    *rtw1=tw1;
    *rtw3=tw3;
    
    
    //cout << "2: " << ievOK << " par0= " << Par[0]<< " par1= " << Par[1] << " par2= " << Par[2] << " par3= " << Par[3] << endl;

    if (type==1) {    //-- electron
      //zrad2->Fill(zero_bin);
      //Nfill++;
      //if (Nfill>16000) break;  //-- Events LIMIT
      //if (iev<ntest)   sig_tst->Fill();
      if (rndm->Rndm()<0.1)   sig_tst->Fill();
      else  signal->Fill();
    
      for (int ip=0; ip<NPF; ip++) { par_e[ip]->Fill(Par[ip]); if (iev==6812 || iev==10724)  printf("t=1 i=%d dEdx = %f \n",ip,Par[ip]); }
    }
    if (type==0) {    //-- pion
      //znorad2->Fill((double)zero_bin);
      //if (iev<ntest)  bg_tst->Fill();
      if (rndm->Rndm()<0.1)   bg_tst->Fill();
      else                    background->Fill();
      for (int ip=0; ip<NPF; ip++) { par_pi[ip]->Fill(Par[ip]); if (iev==6812 || iev==10724)  printf("t=0 i=%d dEdx = %f \n",ip,Par[ip]); }
    }

     
  }//-- gemtrd hit loop
  
  //printf("----------------------------------------------------e_chan1=%d  e_chan2=%d\n",e_chan1,e_chan2);

  //------------ Plotter -------------------------------------

  double escale_trk = 1 ;   if (ntrk_e>0)  escale_trk = (double)ntrk_pi/(double)ntrk_e;
  printf("escale_trk=%f\n",escale_trk);
  int NORM=1; //-- scale hist using no entries
  int DRAW=2; // 1=draw 2= two stat boxes
  int nxd=3;
  int nyd=5;
  int COMPACT=0;
  TCanvas *c1,*c0;
  char ctit[120];
  sprintf(G_DIR,"hd_rawdata_%06d.root",runnum);
  sprintf(ctit,"File=%s",G_DIR);
  //sprintf(rootOut,"hd_rawdata_%06d.root",runnum);
  //sprintf(ctit,"File=%s",rootOut);
  htitle(ctit); //  if (!COMPACT) c1=NextPlot(0,0);

  c1=NextPlot(nxd,nyd);   hits2d_e->Draw("colz");
  TLine *lin1 = new TLine(0.,e_chan1,300.,e_chan1);    TLine *lin2 = new TLine(0.,e_chan2,300.,e_chan2);
  lin1->SetLineColor(kRed);   lin2->SetLineColor(kRed);   lin1->Draw();  lin2->Draw();
  gPad->Modified(); gPad->Update();

  c1=NextPlot(nxd,nyd);   hits2d_p->Draw("colz");
  TLine *lin1p = new TLine(0.,pi_chan1-1,300.,pi_chan1); TLine *lin2p = new TLine(0.,pi_chan2-1,300.,pi_chan2);
  lin1p->SetLineColor(kCyan); lin2p->SetLineColor(kCyan); lin1p->Draw(); lin2p->Draw();
  gPad->Modified(); gPad->Update();

  printf(" Draw Lines :: %d %d %d %d \n",e_chan1,e_chan2, pi_chan1, pi_chan2);
  //sleep(1);

  c1=NextPlot(nxd,nyd);  aver2d_e->Draw("colz");
  c1=NextPlot(nxd,nyd);  aver2d_p->Draw("colz");
  c1=NextPlot(nxd,nyd); hNhits->Draw("colz");

  if  (!WC) {
    //---------------------------------------------------------------------
    cout << " ++++  Amplitude Rejection ++++" << endl;
    double rej70 = Reject(pi_amax, e_amax, 0.7);
    double rej90 = Reject(pi_amax, e_amax, 0.9);
    cout << " Ampl: e=70% , Eff pi = " << rej70*100. << "% ,  Rejection =" << 1./rej70 << endl;
    cout << " Ampl: pi=90%, Eff pi = " << rej90*100. << "% ,  Rejection =" << 1./rej90 << endl;
  }

  //c1=NextPlot(nxd,nyd);    if  (WC) gem2d->Draw("colz"); else hscale(zrad2,znorad2,1.,NORM,2);
  //c1=NextPlot(nxd,nyd);   ampl->Draw("hist");
  //c1=NextPlot(nxd,nyd);   h2xdiff->Draw("colz"); gPad->SetLogz();
  //if  (WC) {
    //c1=NextPlot(nxd,nyd);
    //hscale(rad,norad,0.,NORM,1);
  //}
  if (!WC) {
    c1=NextPlot(nxd,nyd);
    hscale(e_amax,pi_amax,0.,NORM,2);
    c1=NextPlot(nxd,nyd);
    hscale(e_etot,pi_etot,0.,NORM,2);
  }
  c1=NextPlot(nxd,nyd);
  hscale(time_e,time_pi,escale_trk,NORM,2); //--- scale time hist here ---
  //c1=NextPlot(nxd,nyd);   hbeamX->Draw();

  for (int ip=0; ip<NPF; ip++) {
    c1=NextPlot(nxd,nyd);
    if (NN_MODE==0 || (NN_MODE > 1 && ip>=NFixed))  gPad->SetLogy();
    hscale(par_e[ip],par_pi[ip],escale_trk,1,2);  //-- scale no_trk
  }
  c1=NextPlot(-1,-1);
  //------------------------------------------------------------
  //char pngname[120];
  //sprintf(pngname,"mlpOutput/%s_dqm_m%d.png",G_DIR,NN_MODE);
  //c1->Print(pngname);
  char pdfname[120];
  sprintf(pdfname,"mlpOutput/%s_dqm_m%d.pdf",G_DIR,NN_MODE);
  c1->Print(pdfname);

  c0 = new TCanvas("Plot","Plot",400,200,1200,900);     c0->Divide(2,2);
  c0->cd(1);  hscale(time_e,time_pi,1.,0,2); // --- already scaled before
  
  c0->cd(3);  hits2d_e->Draw("colz"); // --- already scaled before
  lin1->Draw();  lin2->Draw();
  gPad->Modified(); gPad->Update();
  
  c0->cd(4);  hits2d_p->Draw("colz"); // --- already scaled before
  lin1p->Draw(); lin2p->Draw();
  gPad->Modified(); gPad->Update();
  
  printf("Draw Lines on c0: e1=%d e2=%d \n",e_chan1,e_chan2);
  c0->cd();
  c0->Modified(); c0->Update();
  char pdfname0[120];
  sprintf(pdfname0,"mlpOutput/%s_dqm_m%d_time.pdf",G_DIR,NN_MODE);
  c0->Print(pdfname0);

  //------------------------------------------------------------
  return NN_MODE;
} //--- fill _trees() ---

//==================================================================
//
//==================================================================

void trd_mlp_fermi(int runnum) {
  
  int WC = 0;
  
  hcount= new TH1D("hcount","Count",3,0,3);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
#if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
  hcount->SetCanExtend(TH1::kXaxis);
#else
  hcount->SetBit(TH1::kCanRebin);
#endif
  
  gStyle->SetTitleFontSize(.06);
  gStyle->SetLabelSize(.05, "XY");
  gStyle->SetTitleSize(0.05,"XY");
  
  char rootfile[256];
  //sprintf(rootfile,"DATA/trd_singleTrackHits_Run_%06d.root",runnum);
  sprintf(rootfile,"RootOutput/trd_singleTrackHits_Run_%06d.root",runnum);
  char basename[120];
  char *hd = strstr(rootfile,"/");
  strncpy(basename,&hd[1],120-1);   char *dot= strstr(basename,"."); *dot=0;
  printf("basename =%s \n",basename);
  
  Int_t ntrain=50;  //-- 100  // epoch
  int Nmod=3;

  // Prepare inputs
  // The 2 trees are merged into one, and a "type" branch,
  // equal to 1 for the signal and 0 for the background is added.
  //const char *fname = "hd_rawdata_000554_000.evio.root";
  TFile *input = 0;
  input = new TFile(rootfile);
  if (!input) return;

  TTree *gem_hits = (TTree *) input->Get("gem_hits");
  
  char mlpname[128];   sprintf(mlpname,"mlpOutput/mlp_run%06d.root",runnum);
  TFile* f = new TFile(mlpname,"RECREATE");
  
  TTree *sig_tst = new TTree("sig_tst", "Filtered Events");
  TTree *bg_tst = new TTree("bg_tst", "Filtered Events");
  TTree *signal = new TTree("signal", "Filtered Events");
  TTree *background = new TTree("background", "Filtered Events");
  TTree *simu = new TTree("simu", "Filtered Events");
  
  char ctit[120];   sprintf(ctit,"hd_rawdata_%06d_000",runnum);
  dispe =  new TH2F("dispe","disp e ",300,-0.5,299.5,200,40.5,240.5);
  disppi = new TH2F("disppi","disp #pi; time ; X strip ",100,0.5,100.5,350,-0.5,349.5);
  dispe->SetStats(0);
  disppi->SetStats(0);
  double fsiz=0.02;
  dispe->GetYaxis()->SetTitleSize(fsiz);   dispe->GetXaxis()->SetTitleSize(fsiz);  dispe->GetZaxis()->SetTitleSize(fsiz);
  dispe->GetYaxis()->SetLabelSize(fsiz);   dispe->GetXaxis()->SetLabelSize(fsiz);  dispe->GetZaxis()->SetLabelSize(fsiz);
  disppi->GetYaxis()->SetTitleSize(fsiz);   disppi->GetXaxis()->SetTitleSize(fsiz);  disppi->GetZaxis()->SetTitleSize(fsiz);
  disppi->GetYaxis()->SetLabelSize(fsiz);   disppi->GetXaxis()->SetLabelSize(fsiz);  disppi->GetZaxis()->SetLabelSize(fsiz);
  
   for (int ip=0; ip<MAXpar; ip++) {
     char bname[80], tname[80];
     sprintf(bname,"par%d",ip);     sprintf(tname,"par%d/F",ip);
     sig_tst->Branch(bname, &Par[ip], tname);
     bg_tst->Branch(bname, &Par[ip], tname);
     signal->Branch(bname, &Par[ip], tname);
     background->Branch(bname, &Par[ip], tname);
     simu->Branch(bname, &Par[ip], tname);
   }
   simu->Branch("type", &type, "type/I");
   simu->Branch("ievent", &ievent, "ievent/I");
   sig_tst->Branch("ievent", &ievent, "ievent/I");
   bg_tst->Branch("ievent", &ievent, "ievent/I");
   sig_tst->Branch("ievent", &channel, "channel/I");
   bg_tst->Branch("ievent", &channel, "channel/I");
   signal->Branch("ievent", &ievent, "ievent/I");
   background->Branch("ievent", &ievent, "ievent/I");
  
   cout << " Get trees signal=" << signal << endl;
   int rtw1,rtw3;
   int nn_mode = fill_trees( gem_hits, signal, background, sig_tst, bg_tst , WC, runnum, &rtw1, &rtw3);
   //signal->Print();
   //background->Print();
   //sig_tst->Print();
   //bg_tst->Print();

   //-----------------------------------------------
   type = 1; //-- electron
   Int_t i;
   for (i = 0; i < signal->GetEntries(); i++) {
      signal->GetEntry(i);
      simu->Fill();
   }
   //-----------------------------------------------
   type = 0; //-- pion
   for (i = 0; i < background->GetEntries(); i++) {
      background->GetEntry(i);
      simu->Fill();
   }
   
   //-----------------------------------------------
   // Build and train the NN par1 is used as a weight since we are primarly
   // interested  by high pt events.
   // The datasets used here are the same as the default ones.
   //----------------------------
   const int NPAR=MAXpar; // 5 , MAXpar=10
   string INL,NNcfg;
   for (int il=0; il<NPAR; il++) {
     stringstream ss;  ss << il;  string si = ss.str();
     INL=INL+"@par"+si;   if (il<(NPAR-1)) INL=INL+",";
   }
   //cout<<" ILN="<<INL<<endl;
   NNcfg=INL+":25:8:type";
   //cout << " NNcfg=" << NNcfg << endl;
   TMultiLayerPerceptron *mlp =
     new TMultiLayerPerceptron(NNcfg.data(),simu,"Entry$%2","(Entry$+1)%2");
  
   //===========================================================================
   // Use TMLPAnalyzer to see what it looks for
   TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis",1600,300,800,900);
   mlpa_canvas->Divide(3,4);
   int ipad=1;
   mlpa_canvas->cd(ipad++);
   TLatex latex;
   char text[120];
   latex.SetTextSize(0.05);
   latex.SetTextAlign(13);  //align at top
   double ystep=0.07, ypos=0.95;
   if   (WC>0) sprintf(text,"WireChamber TRD  Mode=%d",nn_mode);
   else        sprintf(text,"GEM TRD  Mode=%d",nn_mode);
   latex.DrawLatex(0.05,ypos-=ystep,text);
   latex.DrawLatex(0.05,ypos-=ystep,rootfile);
   latex.DrawLatex(0.05,ypos-=ystep,NNcfg.data());

   //===========================================================================
   //=====            Train                                             ========
   //===========================================================================
   mlpa_canvas->cd(ipad++);
   mlp->Train(ntrain, "text,graph,update=10,current");
   mlp->Export("trd_ann","C++");
   //-------------------------------
   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(ipad++);
   // shows how each variable influences the network
   ana.DrawDInputs();
   mlpa_canvas->cd(ipad++);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(ipad++);
   // draws the resulting network
   ana.DrawNetwork(0,"type==1","type==0");
   // Use the NN to plot the results for each sample
   // This will give approx. the same result as DrawNetwork.
   // All entries are used, while DrawNetwork focuses on
   // the test sample. Also the xaxis range is manually set.
   TH1F *bg = new TH1F("bgh", "NN output, single mod", 100, -0.05, 1.1);
   TH1F *sig = new TH1F("sigh", "NN output,single mod",100, -0.05, 1.1);

   char htit[120]; sprintf(htit,"NN output, Nmod=%d",Nmod);
   TH1F *bgm = new TH1F("bgm",htit, 100, -.05, 1.1);
   TH1F *sigm = new TH1F("sigm",htit,100, -.05, 1.1);
   TH1F *err = new TH1F("err","err",256, -.5, 255.5);
   //bg->SetDirectory(0);
   //sig->SetDirectory(0);

   //---------------------------------------------------------------------
   //---------           test net                              -----------
   //---------------------------------------------------------------------

   double tout1=0.44, tout2=0.46;

   int GAUSS=0; double g_mean1=0.2, g_sigma1=0.3,    g_mean2=0.8, g_sigma2=0.3  ;
   int DISP=0, DISP_THR=0;
   if (DISP>0) {
     c2 = new TCanvas("Event Dipsplay","Event Display",1200,10,2000,1000); c2->SetRightMargin(0.15);
     c2->Divide(2,1);
   }
   TH1F *trka = new TH1F("trka", " track angle ", 100, -0.5, 0.5);
   TH2F *trka2 = new TH2F("trka2", " track angle vs x ", 100, -0.5, 0.5, 100, 100., 200.);
   if (runnum==494) trka->GetXaxis()->SetTitle("angle, rad");
   else             trka->GetXaxis()->SetTitle("angle, popugay");
   TRandom3 *rndm = new TRandom3();
   Double_t params[NPAR];
   double sum=0;
   for (i = 0; i < bg_tst->GetEntries(); i++) {
     bg_tst->GetEntry(i);
     for (int ip=0; ip<NPAR; ip++ ) params[ip] = Par[ip];
     double out = mlp->Evaluate(0, params);
     if (GAUSS) out=rndm->Gaus(g_mean1,g_sigma1);
     sum+=out; if (!((i+1)%Nmod)) { bgm->Fill(sum/Nmod); sum=0; }
     bg->Fill(out);
     if (tout1<out&&out<tout2) {
       err->Fill(channel);
       //for (int ip=0; ip<NPAR; ip++) printf(" %f ",params[ip]); printf(" type = %d  out=%f iev=%d \n",type,out,ievent);
     }
     if (out>0.7 && (DISP==1 || DISP==3 )) {
       printf(" pi high : out=%f type=%d iev=%d par=%5.1f %5.0f  %5.1f %5.1f %5.1f  \n",out,type,ievent,Par[0],Par[1],Par[2],Par[3],Par[4]);
       gem_hits->GetEntry(ievent);
       c2->cd(); disppi->Reset();
       char htit[128]; sprintf(htit,"#pi event %d; time 8 ns/bin ; x-strip ",ievent);   disppi->SetTitle(htit);
       for (int i=0;i<gem_nhit;i++){
     if (dedx->at(i)>DISP_THR) disppi->Fill(xpos->at(i),zpos->at(i),dedx->at(i));
     //printf(" %d %f %f  \n",xpos->at(i),zpos->at(i),w2mhit[i]);
       }
       disppi->Draw("colz");
       TLine lin1(0.,rtw1,350.,rtw1);   TLine lin2(0.,rtw3,350.,rtw3);  lin1.SetLineColor(kBlue); lin2.SetLineColor(kBlue); lin1.Draw(); lin2.Draw();
       c2->Modified(); c2->Update(); // usleep(100);  sleep(1);
     }
   }
   sum=0;
   for (i = 0; i < sig_tst->GetEntries(); i++) {
     sig_tst->GetEntry(i);
     for (int ip=0; ip<NPAR; ip++ ) params[ip] = Par[ip];
     double out = mlp->Evaluate(0, params);
     if (GAUSS) out=rndm->Gaus(g_mean2,g_sigma2);
     sum+=out; if (!((i+1)%Nmod)) { sigm->Fill(sum/Nmod); sum=0; }
     sig->Fill(out);
     if (tout1<out&&out<tout2) {
       err->Fill(channel);
       for (int ip=0; ip<NPAR; ip++) printf(" %f ",params[ip]); printf(" type = %d out=%f iev=%d \n",type,out,ievent);
     }
     if (out<0.1 && DISP>1 ) {
       printf(" e low : out=%f type=%d iev=%d par=%5.1f %5.0f  %5.1f %5.1f %5.1f  \n",out,type,ievent,Par[0],Par[1],Par[2],Par[3],Par[4]);
       gem_hits->GetEntry(ievent);
       c2->cd(1);  dispe->Reset();
       char htit[128]; sprintf(htit,"electrons event %d ; time 8 ns/bin ; x-strip",ievent);  dispe->SetTitle(htit);
       int ii=0; double x[1000], y[1000],ey[1000];
       double yaver=0;
       for (int i=0;i<gem_nhit;i++){
     if (dedx->at(i)>DISP_THR) {
       dispe->Fill(zpos->at(i),xpos->at(i),dedx->at(i));
       if ( 110 < zpos->at(i) && zpos->at(i) < 185 ) {
         x[ii]=zpos->at(i);  y[ii]=xpos->at(i); if (dedx->at(i)==0) ey[ii]=100; else ey[ii]=1000./dedx->at(i); ii++; yaver+=y[ii];
       }
     }
     //printf(" %d %f %f  \n",xpos->at(i),zpos->at(i),dedx->at(i));
       }
       yaver/=ii;
       //TLinearFitter *lf=new TLinearFitter(ii);
       //lf->SetFormula("p1");
       TF1 ffit1("ffit1", "pol1", rtw1, rtw3);  ffit1.SetLineColor(kBlue);
       TF1 ffit2("ffit2", "pol1", rtw1, rtw3);  ffit2.SetLineColor(kRed);
           
       TGraphErrors grr(ii, x, y, 0, ey);   grr.SetMinimum(150);   grr.SetMaximum(250);
       Double_t p1=0,p0=0;
       grr.Draw("ap");
       if (ii>2) {
     grr.Fit(&ffit1);       grr.Fit(&ffit2, "+rob=0.75");
     p1=ffit2.GetParameter(1);
     p0=ffit2.GetParameter(0);
     trka2->Fill(p1*2.,p0);
     printf("p0=%f yaver=%f \n",p0,yaver);
       }
       if (runnum==494) trka->Fill(p1*2.); else trka->Fill(p1);
       c2->cd(2);
       dispe->Draw("colz");
       TLine lin1(rtw1,40.,rtw1,240.);   TLine lin2(rtw3,40.,rtw3,240.);  lin1.SetLineColor(kRed); lin2.SetLineColor(kRed); lin1.Draw(); lin2.Draw();
       ffit1.Draw("same");
       ffit2.Draw("same");
       c2->Modified(); c2->Update();  //usleep(100);  sleep(1); getchar();
       c2->WaitPrimitive();
     }
   }
   //---------------------------------------------------------------------
   //---------               plot NN out                      ------------
   //---------------------------------------------------------------------
   mlpa_canvas->cd(ipad++);
   bg->SetLineColor(kBlue);
   bg->SetFillStyle(3345);   bg->SetFillColor(kBlue);
   sig->SetLineColor(kRed);
   sig->SetFillStyle(3354); sig->SetFillColor(kRed);
   bg->SetStats(0);
   sig->SetStats(0);
   bg->Draw();
   sig->Draw("same");
   TLegend *legend = new TLegend(.75, .80, .95, .95);
   legend->AddEntry(bg, "Background (pions)");
   legend->AddEntry(sig, "Signal (electrons)");
   legend->Draw();
   //---------  plot NN mod sum  ------------
   mlpa_canvas->cd(ipad++);
   bgm->SetLineColor(kBlue);
   bgm->SetFillStyle(3345);   bg->SetFillColor(kBlue);
   sigm->SetLineColor(kRed);
   sigm->SetFillStyle(3354); sig->SetFillColor(kRed);
   bgm->SetStats(0);
   sigm->SetStats(0);
   bgm->Draw();
   sigm->Draw("same");
   TLegend *legend2 = new TLegend(.75, .80, .95, .95);
   legend2->AddEntry(bgm, "Background (no rad)");
   legend2->AddEntry(sigm, "Signal (rad)");
   legend2->Draw();
   //---------------------------------------------------------------------
   
   //---------- GRAPH 1 mod ------------------------
   const Int_t ngr = 9;
   const Int_t kNMAX = 9;
   Double_t *Xgr = new Double_t[kNMAX];
   Double_t *Ygr = new Double_t[kNMAX];
   for ( int i=0; i<ngr; i++) {
     double eeff=(i+1)*0.1;
     Xgr[i] = eeff;
     Ygr[i] = Reject(bg, sig, eeff);
     cout << " i=" << i << " x=" << Xgr[i] << " y=" << Ygr[i] << endl;
   }
   TGraph *gr = new TGraph(ngr,Xgr,Ygr); gr->SetName("e #pi efficiency");gr->SetTitle("Efficiency single module");
   /*
   gr->GetXaxis()->SetTitle("electon eff");
   gr->GetYaxis()->SetTitle("#pi efficiency");
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlue);
   gr->SetMarkerSize(0.5);
   */
   //
   TMultiGraph *mg = new TMultiGraph();
   mg->Add(gr,"lp");
   mlpa_canvas->cd(ipad++);
   //gr->Draw("ALP");
   mg->Draw("a");

   //---------------------------------------------------------------------
   /*
   //---------- GRAPH n-mod ------------------------
   const Int_t ngr2 = 9;
   const Int_t kNMAX2 = 9;
   Double_t *Xgr2 = new Double_t[kNMAX2];
   Double_t *Ygr2 = new Double_t[kNMAX2];
   for ( int i=0; i<ngr2; i++) {
     double eeff=0.5+(i+1)*0.05;
     Xgr2[i] = eeff;
     Ygr2[i] = Reject(bgm, sigm, eeff);
     cout << " i=" << i << " x=" << Xgr2[i] << " y=" << Ygr2[i] << endl;
   }
   TGraph *gr2 = new TGraph(ngr2,Xgr2,Ygr2); gr2->SetName("e #pi efficiency");
   char grtit[120]; sprintf(grtit," Efficiency %d modules ",Nmod);   gr2->SetTitle(grtit);
   gr2->GetXaxis()->SetTitle("electon eff");
   gr2->GetYaxis()->SetTitle("#pi efficiency");
   gr2->SetMarkerStyle(20);
   gr2->SetMarkerColor(kBlue);
   gr2->SetMarkerSize(0.5);
   mlpa_canvas->cd(ipad++);
   gr2->Draw("ALP");
   //---------- end GRAPH --------------------
   */
   mlpa_canvas->cd(ipad++);   trka->Draw();
   mlpa_canvas->cd(ipad++);   trka2->Draw("colz");
   mlpa_canvas->cd(ipad++);   err->Draw("colz");
   mlpa_canvas->cd(ipad++);   gPad->SetLogy(); hcount->Draw();
   //---------------------------------------------------------------------
   double rej70 = Reject(bg, sig, 0.7);
   double rej80 = Reject(bg, sig, 0.8);
   double rej85 = Reject(bg, sig, 0.85);
   double rej90 = Reject(bg, sig, 0.9);
   cout << " 1-Mod e=70% , Eff pi = " << rej70*100. << "% ,  Rejection =" << 1./rej70 << endl;
   cout << " 1-Mod pi=90%, Eff pi = " << rej90*100. << "% ,  Rejection =" << 1./rej90 << endl;
   //---------------------------------------------
   double rej70m = Reject(bgm, sigm, 0.7);
   double rej90m = Reject(bgm, sigm, 0.9);
   cout << "Nmod=" << Nmod << " Nmod=" << Nmod << " e=70% , Eff pi = " << rej70m*100. << "% ,  Rejection =" << 1./rej70m << endl;
   cout << " Nmod=" << Nmod << " pi=90%, Eff pi = " << rej90m*100. << "% ,  Rejection =" << 1./rej90m << endl;
   //---------------------------------------------
   mlpa_canvas->cd(1);
   stringstream ss;   ss << " Nmod=" << 1 << " e=70% , Eff #pi = " << rej70*100. << "% ,  Rej =" << 1./rej70 ;  string str = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str.data());
   //--
   ss.str("");  ss.clear();
   ss << " Nmod=" << 1 << " e=80% , Eff #pi = " << rej80*100. << "% ,  Rej =" << 1./rej80 ;  string str0 = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str0.data());
   //--
   ss.str("");  ss.clear();
   ss << " Nmod=" << 1 << " e=85% , Eff #pi = " << rej85*100. << "% ,  Rej =" << 1./rej85 ;  string str1 = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str1.data());
   //--
   ss.str("");  ss.clear();
   ss << " Nmod=" << 1 << " e=90% , Eff #pi = " << rej90*100. << "% ,  Rej =" << 1./rej90 ;  string str2 = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str2.data());
   latex.DrawLatex(0.05,ypos-=ystep,"--------------");
   //--
   ss.str("");  ss.clear();
   ss << " Nmod=" << Nmod << " e=70% , Eff #pi = " << rej70m*100. << "% ,  Rej =" << 1./rej70m ;  str2 = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str2.data());
   //--
   ss.str("");  ss.clear();
   ss << " Nmod=" << Nmod << " e=90% , Eff #pi = " << rej90m*100. << "% ,  Rej =" << 1./rej90m ;  str2 = ss.str();
   latex.DrawLatex(0.05,ypos-=ystep,str2.data());

   //---------------------------------------------
   sprintf(text,"mlpOutput/%s_m%d.pdf",basename,nn_mode);
   mlpa_canvas->Print(text);
   mlpa_canvas->cd(0);
 
   //f->Write(0, TObject::kWriteDelete);
  f->Write();
 
   delete input;
}


//==================================================================
void Count(const char *tit) {
  hcount->Fill(tit,1);
}
void Count(const char *tit, double cut1) {
  char clab[20];
  sprintf(clab,"%s_%.1f",tit,cut1);
  hcount->Fill(clab,1);
}
void Count(const char *tit, double cut1, double cut2) {
  char clab[20];
  sprintf(clab,"%s_%.1f_%.1f",tit,cut1,cut2);
  hcount->Fill(clab,1);
}
//------------------------------------------------------------------

