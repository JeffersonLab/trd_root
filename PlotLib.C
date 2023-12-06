//************************************************************************
//                   P L O T   L I B R A R Y  
//************************************************************************
#include "TPostScript.h"
#include "TPDF.h"
#include "TPaveLabel.h"

TCanvas *NextPlot(int nx, int ny);
void htitle(const char *tit);
//---------------------------------------------------
//        GLOBAL VARIABLES
//---------------------------------------------------
int G_RUN1=0, G_RUN2=0, G_VERS;
char G_DIR[128];
#define USE_PDF
#ifdef USE_PDF
TPDF *ps;
#else
TPostScript *ps;
#endif
int  PS=1;      //--  ON/OFF postscript file
int  HTIT=0;    //--  User Title
char HTitle[128];
int NORM=1;

//------------------------------------------------------------------------

void Process(TFile * hfile,TH1 *hist) {
  int noent= hist->GetEntries();
  if (noent>0 && NORM)  hist->Scale(1./noent);
}

//----------  ver 4.0 -------
TCanvas *NextPlot(int nx, int ny) {
  
  static int iplot, ipad, ican, LX, LY, FIRST, EXIT;
  static TCanvas *c1;
  static TPad *pad[50*50];
  static char Cname[128], PSname[128];
  memset(Cname,0,sizeof(Cname));
  
  int run1=G_RUN1;
  int run2=G_RUN2;
  
  if (nx>0 && ny>0) { LX=nx; LY=ny; }
  if (ican>0 && nx==0 && ny==0) { ipad=LX*LY; return c1; }
  if (ican>0 && nx<0 && ny<0) { EXIT=1; }
  
  if (ipad==0 || ipad>=(LX*LY)  || EXIT ) { //--- next canvas --
    if (ican>0 && PS>0 ) c1->Update();
    if (EXIT && PS>0 ) {  //-- print and close previous ps  --
      ps->Close();
      char CMD[580]; sprintf(CMD,"evince  %s &",PSname);
      //printf("Exec PSname=%s  CMD=%s \n",PSname, CMD);
      gSystem->Exec(CMD);
    }
    if (EXIT) {  EXIT=0; iplot=0; ipad=0; ican=0; LX=0; LY=0;  FIRST=0;  return c1; }
    ican++; ipad=0;
    //----------------------------------------------------------
    char sltime[256];
    struct tm *date;
    time_t  timer;
    time(&timer);
    //printf("timer=%ld \n",timer);
    strcpy(sltime,ctime(&timer));
    //printf("Time = %s",sltime);
    date=localtime(&timer);
    strftime(sltime,80," %d-%b-%Y  %H:%M:%S ",date);
    sprintf(Cname,"%s    page %d    %s",G_DIR,ican,sltime);
    //printf("open next canvas:%s, plot=%d, pad=%d  ican=%d \n",Cname,iplot,ipad,ican);
    //----------------------------------------------------------
    if (iplot==0 || PS==0) { //------  Create Canvas --- Open PS file ---
      c1=new TCanvas(Cname,Cname,50+(ican-1)*10,10+(ican-1)*10,695*1.3,900*1.3); //- Letter 695x900  +30%
      c1->SetFillColor(0); 
      if (iplot==0 && PS>0) {
#ifdef USE_PDF
	      sprintf(PSname,"%s-v%d.pdf",G_DIR,G_VERS);
	      ps = new TPDF(PSname,111);
#elif USE_PDF_MLP
        sprintf(PSname,"mlpOutput/%s-v%d.pdf",G_DIR,G_VERS);
        ps = new TPDF(PSname,111);
#else
        sprintf(PSname,"%s-v%d.ps",G_DIR,G_VERS);
        ps = new TPostScript(PSname,-100111);
#endif
	      printf("Opening new file PSname=%s for histograms...\n",PSname);
      }
    } 
#ifndef USE_PDF
    if (iplot>0 && PS>0)  ps->NewPage();
#endif
    if (PS>0)  c1->Clear();
    
    //----------- Draw a global picture title -----------
    TPaveLabel *title = new TPaveLabel(0.03,0.98,0.97,0.99,Cname);
    title->SetFillColor(0);
    title->SetTextFont(52);
    title->SetTextSize(0.9);
    title->SetBorderSize(0);
    title->Draw();
    //----------- user title -----------------
    if (HTIT) {
      TPaveLabel *title = new TPaveLabel(0.03,0.95,0.97,0.97,HTitle);
      title->SetFillColor(0);
      title->SetTextFont(52);
      title->SetTextSize(0.9);
      title->SetBorderSize(1);
      title->Draw();
    }
    //--------------  Create PADs ::  zon 3 6  ----------------
    //printf("Create PADs:: ipad=%d  LY=%d LX=%d \n",ipad,LY,LX);
    int ip=0;
    for (int ypad=(LY-1);ypad>=0;ypad--) {
      for (int xpad=0;xpad<LX;xpad++) {
	      char Cpad[80];
	      sprintf(Cpad,"Pad%d-%d-%d",xpad,ypad,ican);
	      double xp1=xpad*1./LX;	     double xp2=(xpad+1)*1./LX;
	      double yp1=ypad*0.94/LY;     double yp2=(ypad+1)*0.94/LY;
	      pad[ip] = new TPad(Cpad,"The pad",xp1,yp1,xp2,yp2,0);
	      pad[ip]->Draw();
	      ip++;
      }
    }
  } //-- if not next canvas -->> if (ipad==0 || ipad>=(LX*LY)  || EXIT )
  pad[ipad]->cd();
  ipad++;
  iplot++;
  return c1;
}
//---------------------------------------------------
void htitle(const char *tit) {
  if (tit) { HTIT=1; sprintf(HTitle,"%s",tit);
  } else { HTIT=0; }
}
//---------------------------------------------------
