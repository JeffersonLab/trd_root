//************************************************************************
//                   P L O T   L I B R A R Y  
//************************************************************************
#include "TPostScript.h"
#include "TPDF.h"
#include "TPaveLabel.h"

TCanvas *NextPlot(int nx, int ny);
void htitle(const char *tit);
//---------------------------------------------------
//        GLOBAL VARIABLEs
//---------------------------------------------------
int G_RUN1=0,  G_RUN2=0, G_VERS;
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
  //  printf(" Enter Process .... \n");
  //TH1 *hist_norm;
  //hist_norm = (TH1*)hfile->Get("h2_dedx");
  //int noent= hist_norm->GetEntries(); 
  int noent= hist->GetEntries(); 
  //printf(" noent=%d \n",noent);
  if (noent>0 && NORM)  hist->Scale(1./noent);

}

//-------------------------------------------------------  ver 4.0 -------
//------------------------------------------------------------------------
TCanvas *NextPlot(int nx, int ny) {

  static int iplot,ipad, ican, LX,LY,FIRST,EXIT;
  static TCanvas *c1;
  static TPad *pad[50*50];
  static char Cname[128], PSname[128];   
  memset(Cname,0,sizeof(Cname)); 
  
  printf("next0: Plot=%d, pad=%d  ican=%d\n",iplot,ipad,ican);

  int run1=G_RUN1;
  int run2=G_RUN2;
  
  if (nx>0 && ny>0) { LX=nx; LY=ny; }
  if (ican>0 && nx==0 && ny==0) { ipad=LX*LY; /*c1->Update();*/ return c1; }
  if (ican>0 && nx<0 && ny<0) {  EXIT=1;  }
  
  if (ipad==0 || ipad>=(LX*LY)  || EXIT ) {  //--- next canvas --
    if (ican>0 && PS>0 ) c1->Update();
    if (EXIT && PS>0 ) {  //-- print and close previous ps  --
      //c1->Update();
      //ps->NewPage();
      ps->Close();
      char CMD[580]; sprintf(CMD,"evince  %s &",PSname);
      printf("Exec PSname=%s  CMD=%s \n",PSname, CMD);
      gSystem->Exec(CMD);
    } 
    if (EXIT) {  EXIT=0; iplot=0; ipad=0; ican=0; LX=0; LY=0;  FIRST=0;  return c1; }
    //--
    ican++; ipad=0;
    //----------------------------------------------------------
    char sltime[256];
    struct tm *date; time_t  timer;   time(&timer);    printf(" timer=%ld \n",timer);
    strcpy(sltime,ctime(&timer));    printf(" Date=%s",sltime);      
    date=localtime(&timer);    strftime(sltime,80," %d-%b-%Y  %H:%M:%S ",date);
    //printf(" DATE=%s \n",sltime);

    //sprintf(Cname,"%s   Run-%d-%d-p.%d    %s",G_DIR,run1,run2,ican,sltime);
    sprintf(Cname,"%s  p.%d    %s",G_DIR,ican,sltime);
    printf("open next canvas:%s, plot=%d, pad=%d  ican=%d \n",Cname,iplot,ipad,ican);
    //----------------------------------------------------------
    if (iplot==0 || PS==0) { //------  Create Canvas --- Open PS file ---
      //if (iplot>=0 ) { //------  Create Canvas --- Open PS file ---
      c1=new TCanvas(Cname,Cname,50+(ican-1)*10,10+(ican-1)*10,695*1.3,900*1.3); //- Letter 695x900  +30%
      //c1=new TCanvas(Cname,Cname,50+(ican-1)*10,10+(ican-1)*10,800,800);            
      c1->SetFillColor(0); 
      if (iplot==0 && PS>0) {
#ifdef USE_PDF
	sprintf(PSname,"%s-v%d.pdf",G_DIR,G_VERS);
	ps = new TPDF(PSname,111);
	//ps = new TPDF(PSname,-100111); // Letter
#else
	sprintf(PSname,"%s-v%d.ps",G_DIR,G_VERS);
	ps = new TPostScript(PSname,-100111);
#endif
	//ps->Range(30.,30.);
	printf("Open new file PSname=%s \n",PSname);
      }
    } 
#ifndef USE_PDF
    if (iplot>0 && PS>0)  ps->NewPage();
    //if (PS>0)  c1->Update();
    //if (PS>0)  ps->NewPage();
#endif
    if (PS>0)  c1->Clear();
    
    //----------- Draw a global picture title -----------
    printf("Draw a global picture title ipad=%d  LY=%d LX=%d \n",ipad,LY,LX);
    TPaveLabel *title = new TPaveLabel(0.03,0.98,0.97,0.99,Cname);
    title->SetFillColor(0);    title->SetTextFont(52);    title->SetTextSize(0.9);
    title->SetBorderSize(0);
    title->Draw();
    //----------- user title -----------------
    if (HTIT) {
      TPaveLabel *title = new TPaveLabel(0.03,0.95,0.97,0.97,HTitle);
      title->SetFillColor(0);    title->SetTextFont(52);    title->SetTextSize(0.9);
      title->SetBorderSize(1);
      title->Draw();
    }
    //--------------  Create PADs ::  zon 3 6  ----------------
    printf("Create PADs:: ipad=%d  LY=%d LX=%d \n",ipad,LY,LX);
    int ip=0;
    for (int ypad=(LY-1);ypad>=0;ypad--) {
      for (int xpad=0;xpad<LX;xpad++) {
	char Cpad[80];
	sprintf(Cpad,"Pad%d-%d-%d",xpad,ypad,ican);
	double xp1=xpad*1./LX;	     double xp2=(xpad+1)*1./LX;
	double yp1=ypad*0.94/LY;     double yp2=(ypad+1)*0.94/LY;
	//printf("new pad=%d %f %f %f %f\n",ip,xp1,xp2,yp1,yp2);
	pad[ip] = new TPad(Cpad,"The pad",xp1,yp1,xp2,yp2,0);
	pad[ip]->Draw();
	//pad[ip]->Update();
	ip++;
      }
    }
  } //-- if not next canvas -->> if (ipad==0 || ipad>=(LX*LY)  || EXIT )
  pad[ipad]->cd(); 
  //gPad->SetGrid();
  //gPad->SetFillColor(0);
  //gStyle->SetFillColor(0);

  printf("next: Plot=%d, pad=%d  ican=%d\n",iplot,ipad,ican);
  ipad++; iplot++;   
  return c1;
}
//---------------------------------------------------
void htitle(const char *tit) {
  if (tit) { HTIT=1; sprintf(HTitle,"%s",tit);
  } else   { HTIT=0; }
}
//---------------------------------------------------
