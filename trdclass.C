#define trdclass_cxx
#include "trdclass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "PlotLib.C"

#define NPRT 1000
//#define NPRT 10
#define USE_TRK

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

//-----------------  canvas 0 Event Display ----------
  //char c0Title[256]; sprintf(c0Title,"Event_Display_Run=%d",RunNum);
  //TCanvas *c0 = new TCanvas("DISP",c0Title,200,200,1500,1300);
  //c0->Divide(3,2); c0->cd(1);

  TF1 fx("fx","pol1",100,190);
  TF1 fx1("fx1","pol1",100,190);
  TF1 fx2("fx2","pol1",100,190);

  hCal_occ = new TH1F("hCal_occ"," Calorimeter Occupancy",9,-0.5,8.5);
  hCal_sum = new TH1F("hCal_sum"," Calorimeter Sum",4096,0.5,4095.5);
  for (int cc=0; cc<NCAL; cc++) {
  	  char hName[128];  sprintf(hName,"hCal_adc%d",cc);
	  char hTitle[128]; sprintf(hTitle,"Calorimeter ADC, cell%d",cc);
	  hCal_adc[cc] = new TH1F(hName,hTitle,4096,-0.5,4095.5);
  }

  hcount= new TH1D("hcount","Count",3,0,3);
  hcount->SetStats(0);   hcount->SetFillColor(38);   hcount->SetMinimum(1.);
  #if ROOT_VERSION_CODE > ROOT_VERSION(6,0,0)
  hcount->SetCanExtend(TH1::kXaxis);
  #else
  hcount->SetBit(TH1::kCanRebin);
  #endif

  h250_size = new TH1F("h250_size"," fa250 Raw data size",4096,0.5,4095.5);
	
  hCal_sum_el = new TH1F("hCal_sum_el"," Calorimeter Sum for electrons",4096,0.5,4095.5);
  hCal_sum_pi = new TH1F("hCal_sum_pi"," Calorimeter Sum for pions",4096,0.5,4095.5);
  hCher_u_adc = new TH1F("hCher_u_adc"," Cherenkov Upstream ADC ; ADC Amplitude ",4096,-0.5,4095.5);
  hCher_din_adc = new TH1F("hCher_din_adc"," Cherenkov Downstream (in) ADC ; ADC Amplitude ",4096,-0.5,4095.5);
  hCher_dout_adc = new TH1F("hCher_dout_adc"," Cherenkov Downstream (out) ADC ; ADC Amplitude ",4096,-0.5,4095.5);
  
  hCher_u_time = new TH1F("hCher_u_time"," Cherenkov Upstream Time Response ",300,-0.5,299.5);
  hCher_din_time = new TH1F("hCher_din_time"," Cherenkov Downstream (in) Time Response ",300,-0.5,299.5);
  hCher_dout_time = new TH1F("hCher_dout_time"," Cherenkov Downstream (out) Time Response ",300,-0.5,299.5);
  
  hCCCor_u = new TH2F("hCCCor_u"," Cherenkov Calorimeter Corr ; Upstream ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);
  hCCCor_dout = new TH2F("hCCCor_dout"," Cherenkov Calorimeter Corr ; Downstream (out) ; Calorimeter ",400,0.5,4095.5,400,0.5,4095.5);
  hCCor_ud = new TH2F("hCCor_ud"," Cherenkov Upstream/Downstream Corr ; Upstream ; Downstream (out) ",400,-0.5,4095.5,400,0.5,4095.5);
  
  f125_el = new TH1F("f125_el","GEM-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);
  f125_pi = new TH1F("f125_pi","GEM-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);
  mmg1_f125_el = new TH1F("mmg1_f125_el","MMG1-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);
  mmg1_f125_pi = new TH1F("mmg1_f125_pi","MMG1-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);
  urw_f125_el = new TH1F("urw_f125_el","uRW-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);
  urw_f125_pi = new TH1F("urw_f125_pi","uRW-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);
  mmg2_f125_el = new TH1F("mmg2_f125_el","MMG2-TRD f125 Peak Amp for Electrons ; ADC Amplitude ; Counts ",100,0.,4096);
  mmg2_f125_pi = new TH1F("mmg2_f125_pi","MMG2-TRD f125 Peak Amp for Pions ; ADC Amplitude ; Counts ",100,0.,4096);

  f125_el_chi2 = new  TH1F("f125_el_chi2","TRK_el_chi2",100,0.,10000.);
  f125_pi_chi2 = new  TH1F("f125_pi_chi2","TRK_pi_chi2",100,0.,10000.);
  
  int THRESH=200;
  f125_el_evt = new TH2F("f125_el_evt","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);
  f125_el_evt->SetStats(0);  f125_el_evt->SetMinimum(THRESH);  f125_el_evt->SetMaximum(1000.);
  f125_pi_evt = new TH2F("f125_pi_evt","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);
  f125_pi_evt->SetStats(0); f125_pi_evt->SetMinimum(THRESH); f125_pi_evt->SetMaximum(1000.);
  f125_el_raw = new TH2F("f125_el_raw","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);
  f125_el_raw->SetStats(0);  f125_el_raw->SetMinimum(THRESH);   f125_el_raw->SetMaximum(1000.);
  f125_pi_raw = new TH2F("f125_pi_raw","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",100,100.5,200.5,160,40.5,200.5);
  f125_pi_raw->SetStats(0); f125_pi_raw->SetMinimum(THRESH); f125_pi_raw->SetMaximum(1000.);
  
  f125_el_fit = new TH2F("f125_el_fit","GEM-TRD track for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  f125_pi_fit = new TH2F("f125_pi_fit","GEM-TRD track for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);

  f125_el_amp2ds = new TH2F("f125_el_amp2ds","GEM-TRD Amp for Single Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  f125_pi_amp2ds = new TH2F("f125_pi_amp2ds","GEM-TRD Amp for Single Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  f125_el_amp2d = new TH2F("f125_el_amp2d","GEM-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  f125_pi_amp2d = new TH2F("f125_pi_amp2d","GEM-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);

  mmg1_f125_el_amp2d = new TH2F("mmg1_f125_el_amp2d","MMG1-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  mmg1_f125_pi_amp2d = new TH2F("mmg1_f125_pi_amp2d","MMG1-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  mmg2_f125_el_amp2d = new TH2F("mmg2_f125_el_amp2d","MMG2-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  mmg2_f125_pi_amp2d = new TH2F("mmg2_f125_pi_amp2d","MMG2-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  urw_f125_el_amp2d = new TH2F("urw_f125_el_amp2d","uRW-TRD Amp for Electrons ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  urw_f125_pi_amp2d = new TH2F("urw_f125_pi_amp2d","uRW-TRD Amp for Pions ; Time Response (8ns) ; Channel ",250,0.5,250.5,240,0.5,240.5);
  f125_el_clu2d = new TH2F("f125_el_clu2d","GEM-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);
  f125_pi_clu2d = new TH2F("f125_pi_clu2d","GEM-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);
  mmg1_f125_el_clu2d = new TH2F("mmg1_f125_el_clu2d","MMG1-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);
  mmg1_f125_pi_clu2d = new TH2F("mmg1_f125_pi_clu2d","MMG1-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);
  mmg2_f125_el_clu2d = new TH2F("mmg2_f125_el_clu2d","MMG2-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);
  mmg2_f125_pi_clu2d = new TH2F("mmg2_f125_pi_clu2d","MMG2-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);
  urw_f125_el_clu2d = new TH2F("urw_f125_el_clu2d","uRW-TRD Amp for Electrons (Clusters)",200,0.5,200.5,240,0.5,240.5);
  urw_f125_pi_clu2d = new TH2F("urw_f125_pi_clu2d","uRW-TRD Amp for Pions (Clusters)",200,0.5,200.5,240,0.5,240.5);

//=========================================

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  if (MaxEvt>0) nentries=MaxEvt;  //-- limit number of events for test
  int MAX_PRINT=2; //--- debug printing ---

//================ Begin Event Loop ==============

  int N_trk_el=0;
  int N_trk_pi=0;

  Long64_t jentry=0;
  for (jentry=0; jentry<nentries; jentry++) {
	
  	  //printf("------ Next event %lld ---\n",jentry );
	  Long64_t ientry = LoadTree(jentry);
	  if (ientry < 0) break;
	  nb = fChain->GetEntry(jentry);   nbytes += nb;
	  // if (Cut(ientry) < 0) continue;
	  
      if (jentry<MAX_PRINT || !(jentry%NPRT))
	  printf("------- evt=%llu  f125_raw_count=%llu f125_pulse_count=%llu f250_wraw_count=%llu, srs_raw_count=%llu \n"
		 ,jentry,f125_wraw_count, f125_pulse_count, f250_wraw_count,srs_raw_count);
	
//==================================================================================================
//                    Show Event
//==================================================================================================
	/*
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
	        printf("F125:raw: i=%lldd  sl=%d, ch=%d, idx=%d, cnt=%d \n"
		   ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
		   ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
	    }
	}
*/
//==================================================================================================
//                    Process Fa250  Pulse data
//==================================================================================================

	//#define USE_250_PULSE
	#ifdef USE_250_PULSE
/*
	printf("-------------------- Pulse  250  n=%lld ---------------------------\n",f250_pulse_count);
	for (int i=0;i<f250_pulse_count; i++) {
	  printf("F250:: i=%d  sl=%d, ch=%d,  npk=%d  time=%d amp=%d ped=%f \n"
		 ,i,f250_pulse_slot->at(i),f250_pulse_channel->at(i),f250_pulse_pulse_number->at(i)
		 ,f250_pulse_course_time->at(i),f250_pulse_pulse_peak->at(i),f250_pulse_pedestal->at(i)/4.);
	} */
	#endif
 
//==================================================================================================
//                    Process Fa250  RAW data
//==================================================================================================

	//if (jentry<MAX_PRINT) printf("------------------ Fadc250  wraw_count = %llu ---------\n", f250_wraw_count);
	h250_size->Fill(f250_wraw_count);
	
	double CalSum=0;
	double Ch_u=0;
	double Ch_in=0;
	double Ch_out=0;
	bool electron=false;
	
	for (ULong64_t i=0; i<f250_wraw_count; i++) { // --- fadc250 channels loop
	  /*  if (jentry<MAX_PRINT) printf("F250:: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
				   ,i,f250_wraw_slot->at(i),f250_wraw_channel->at(i)
				   ,f250_wraw_samples_index->at(i),f250_wraw_samples_count->at(i));
	    */
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
	        hCal_adc[fadc_chan]->Fill(amax);
	        CalSum+=amax;
	    } else { // Cherenkov
	        if (fadc_chan==13) { hCher_u_adc->Fill(amax);   hCher_u_time->Fill(tmax); Ch_u=amax; }
	        if (fadc_chan==14) { hCher_din_adc->Fill(amax);  hCher_din_time->Fill(tmax); Ch_in=amax; }
	        if (fadc_chan==15) { if(amax>300)electron=true; hCher_dout_adc->Fill(amax);  hCher_dout_time->Fill(tmax);Ch_out=amax; }
	    }
	} // -- end of fadc250 channels loop
	
//=======================  End Fa250 RAW  process Loop  =====================================================
	
	if(electron){
		hCal_sum_el->Fill(CalSum);
	} else {
		hCal_sum_pi->Fill(CalSum);
	}
	
//==================================================================================================
//                    Process Fa125  Pulse  data
//==================================================================================================
	if (!(jentry%NPRT)) {
		if(electron) {
			f125_el_evt->Reset();
			f125_el_raw->Reset();
		} else {
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
	
	//------  Track FIT -----	
	int gem_chan_max=-1;
	int mmg_chan_max=-1;
	int urw_chan_max=-1;
	float f125_amp_max=0.;
	float mmg_f125_amp_max=0.;
	float urw_f125_amp_max=0.;
	
	for (ULong64_t i=0;i<f125_pulse_count; i++) {
	/*    if (jentry<MAX_PRINT) printf("F125:: i=%lld  sl=%d, ch=%d, npk=%d time=%d amp=%d ped=%d \n"
				   ,i,f125_pulse_slot->at(i),f125_pulse_channel->at(i),f125_pulse_npk->at(i)
				   ,f125_pulse_peak_time->at(i),f125_pulse_peak_amp->at(i),f125_pulse_pedestal->at(i));
	    //cout<<" ++++++++++++++++++++ f125_pulse_npk= "<<f125_pulse_npk->at(i)<<endl;
	*/	
		
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
        if(electron){
            if (gemChan>-1) {
                #ifdef USE_TRK
                //f125_el_amp2ds->Fill(time,gemChan,amp);
                f125_el_fit->Fill(time,gemChan,amp);
                //#else
                //f125_el_amp2d->Fill(time,gemChan,amp);
                #endif
			} else {
   		        if (gemChan>-1) {
	                 #ifdef USE_TRK
    	             //f125_pi_amp2ds->Fill(time,gemChan,amp);
        	         f125_pi_fit->Fill(time,gemChan,amp);
            	     //#else
                	 //f125_pi_amp2d->Fill(time,gemChan,amp);
                 	#endif
				}
			}
		}
	} //-- End f125 pulse loop --
	
	//------  Track FIT -----
	bool isSingleTrack=false;
	
	#ifdef USE_TRK
	double chi2_max=3000;
	if (electron) {
	    double chi2el = TrkFit(f125_el_fit, fx, "fx");
	    f125_el_chi2->Fill(chi2el);
	    //if ( 1 || chi2el>0. && chi2el<chi2_max) {
	    if (chi2el>0. && chi2el<chi2_max) {
	    	//f125_el_amp2d->Add(f125_el_amp2ds);
	        Count("n_trk_el");
	        N_trk_el++;
	        isSingleTrack=true;
	    }
	} else {
		double chi2pi = TrkFit(f125_pi_fit,fx,"fx");
	    f125_pi_chi2->Fill(chi2pi);
	    //if ( 1 || chi2pi>0. && chi2pi<chi2_max) {
	    if (chi2pi>0. && chi2pi<chi2_max) {
	        //f125_pi_amp2d->Add(f125_pi_amp2ds);
	        Count("n_trk_pi");
	        N_trk_pi++;
	        isSingleTrack=true;
	    }
	}
	#endif
		
	//----- Single Track Event Hist Filling -----
	if (isSingleTrack) {
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
	    	if(electron){
	        	if (gemChan>-1) {
					//#ifdef USE_TRK
        			//f125_el_amp2ds->Fill(time,gemChan,amp);
					//f125_el_fit->Fill(time,gemChan,amp);
					//#else
	            	f125_el_amp2d->Fill(time,gemChan,amp);
					//#endif
					if (!(jentry%NPRT)) f125_el_evt->Fill(time,gemChan,amp);
	            	f125_el->Fill(amp);
	            	f125_el_clu2d->Fill(time,gemChan,1.);
	        	}
	        	if (amp>MM_THR && mmg1Chan>-1) {
	            	mmg1_f125_el_amp2d->Fill(time,mmg1Chan,amp);
	            	mmg1_f125_el->Fill(amp);
	            	mmg1_f125_el_clu2d->Fill(time,mmg1Chan,1.);
	        	}
	        	if (mmg2Chan>-1) {
		        	mmg2_f125_el_amp2d->Fill(time,mmg2Chan,amp);
		        	mmg2_f125_el->Fill(amp);
		   		    mmg2_f125_el_clu2d->Fill(time,mmg2Chan,1.);
	        	}
	        	if (rwellChan>-1) {
	            	urw_f125_el_amp2d->Fill(time,rwellChan,amp);
	            	urw_f125_el->Fill(amp);
	            	urw_f125_el_clu2d->Fill(time,rwellChan,1.);
	        	}
	    	} else {
	        	if (gemChan>-1) {
					//#ifdef USE_TRK
        			//f125_pi_amp2ds->Fill(time,gemChan,amp);
					//f125_pi_fit->Fill(time,gemChan,amp);
					//#else
	            	f125_pi_amp2d->Fill(time,gemChan,amp);
					//#endif
					if (!(jentry%NPRT)) f125_pi_evt->Fill(time,gemChan,amp);
	            	f125_pi->Fill(amp);
	            	f125_pi_clu2d->Fill(time,gemChan,1.);
	        	}
	        	if (amp>MM_THR && mmg1Chan>-1) {
	            	mmg1_f125_pi_amp2d->Fill(time,mmg1Chan,amp);
	            	mmg1_f125_pi_clu2d->Fill(time,mmg1Chan,1.);
	            	mmg1_f125_pi->Fill(amp);
	        	}
	        	if (mmg2Chan>-1) {
		        	mmg2_f125_pi_amp2d->Fill(time,mmg2Chan,amp);
		        	mmg2_f125_pi->Fill(amp);
		        	mmg2_f125_pi_clu2d->Fill(time,mmg2Chan,1.);
	        	}
	        	if (rwellChan>-1) {
	            	urw_f125_pi_amp2d->Fill(time,rwellChan,amp);
	            	urw_f125_pi_clu2d->Fill(time,rwellChan,1.);
	            	urw_f125_pi->Fill(amp);
	        	}
	    	}
			
	    	//if (peak_amp-ped>f125_amp_max) {
	    	//    f125_amp_max=peak_amp-ped;
	    	//    gem_chan_max = fADCChan;
	    	//}
			
	    	hCal_sum->Fill(CalSum/7.);
	    	hCCor_ud->Fill(Ch_u,Ch_out);
	    	hCCCor_u->Fill(Ch_u,CalSum/7.);
	    	hCCCor_dout->Fill(Ch_out,CalSum/7.);
			
		} //--- end Fa125 Pulse Loop ---
	} //-- end single track condition --
	
//======================= End Process Fa125 Pulse data ================================

//==================================================================================================
//                    Process Fa125  RAW data
//==================================================================================================

	//#define USE_125_RAW
    #ifdef USE_125_RAW
	//if (jentry<MAX_PRINT) printf("------------------ Fadc125  wraw_count = %llu ---------\n", f125_wraw_count);
	
    for (ULong64_t i=0;i<f125_wraw_count; i++) { // --- fadc125 channels loop
    /*  if (jentry<MAX_PRINT) printf("F125:RAW: i=%lld  sl=%d, ch=%d, idx=%d, cnt=%d \n"
				   ,i,f125_wraw_slot->at(i),f125_wraw_channel->at(i)
				   ,f125_wraw_samples_index->at(i),f125_wraw_samples_count->at(i));
	*/
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
				}else {
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
	  			} else {
	    			double cc = f125_pi_raw->GetBinContent(ii, jj);
	    			if (cc == 0.) f125_pi_raw->Fill(ii,jj,pedestal);
	  			}
			}
		}
	}

//=======================  End Fa125 RAW  process Loop  =====================================================

/*    if(electron){
	  f125_el->Fill(f125_amp_max);
	  //if((slot<6)||(slot==7&&chan<24))f125_el->Fill(f125_amp_max);
	  //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_el->Fill(f125_amp_max);
	} else {
	  f125_pi->Fill(f125_amp_max);
	  //if((slot<6)||(slot==7&&chan<24))f125_pi->Fill(f125_amp_max);
	  //if((slot==6&&chan>23)||(slot>6&&slot<9)||(slot==9&&chan<48))mmg_f125_pi->Fill(f125_amp_max);
	}*/

    #endif
	  
//=====================================================================================
//===                     Event Display                                            ====
//=====================================================================================
	/*
	if (jentry<MAX_PRINT || !(jentry%NPRT)) {
	  c0->cd(1); f125_el_amp2d->Draw("colz");
	  c0->cd(4); f125_pi_amp2d->Draw("colz");
	  c0->cd(2); f125_el_evt->Draw("colz");         double chi2e = TrkFit(f125_el_evt,fx1,"fx1"); if (chi2e>0. && chi2e < 1000. ) fx1.Draw("same");
	  c0->cd(5); f125_pi_evt->Draw("colz");         double chi2p = TrkFit(f125_pi_evt,fx2,"fx2"); if (chi2p>0. && chi2p < 1000. ) fx2.Draw("same");
	printf("========================>>>  Chi2 e=%f p=%f \n",chi2e,chi2p);
      //c0->cd(3); f125_el_chi2->Draw("colz");
      //c0->cd(6); f125_pi_chi2->Draw("colz");
      c0->cd(3); f125_el_raw->Draw("colz");  f125_el_evt->Draw("same");
      c0->cd(6); f125_pi_raw->Draw("colz");  f125_pi_evt->Draw("same");

      c0->Modified();   c0->Update();
      if (NPRT<1000) sleep(1);
	}
	*/
  } // -- end of event loop
   
  cout<<" Total events= "<<jentry<< "  N_trk_el=" << N_trk_el++ << " N_trk_pi=" << N_trk_pi <<endl;
   
  //=====================================================================================
  //===                 S A V E   H I S T O G R A M S                                ====
  //=====================================================================================

  //open
  TFile* fOut;
  //char rootFileName[256]; sprintf(rootFileName, "FNAL_JANA2/RunOutput/Run_%06d_SingleTrkOutput.root", RunNum);
  char rootFileName[256]; sprintf(rootFileName, "FNAL_JANA2/RunOutput/singleTrkEventsTree_Run%06d_Output_%06dEntries.root", RunNum, jentry);
  fOut = new TFile(rootFileName, "RECREATE");
  fOut->cd();
  //write hists
  hcount->Write(0, TObject::kOverwrite);
  h250_size->Write(0, TObject::kOverwrite);
  hCher_u_adc->Write(0, TObject::kOverwrite);
  hCher_din_adc->Write(0, TObject::kOverwrite);
  hCher_dout_adc->Write(0, TObject::kOverwrite);
  hCher_u_time->Write(0, TObject::kOverwrite);
  hCher_din_time->Write(0, TObject::kOverwrite);
  hCher_dout_time->Write(0, TObject::kOverwrite);
  hCCCor_u->Write(0, TObject::kOverwrite);
  hCCor_ud->Write(0, TObject::kOverwrite);
  hCCCor_dout->Write(0, TObject::kOverwrite);
  hCal_occ->Write(0, TObject::kOverwrite);
  hCal_sum_el->Write(0, TObject::kOverwrite);
  hCal_sum_pi->Write(0, TObject::kOverwrite);
  f125_el_chi2->Write(0, TObject::kOverwrite);
  f125_pi_chi2->Write(0, TObject::kOverwrite);
  f125_el->Write(0, TObject::kOverwrite);
  f125_pi->Write(0, TObject::kOverwrite);
  f125_el_amp2d->Write(0, TObject::kOverwrite);
  f125_pi_amp2d->Write(0, TObject::kOverwrite);
  f125_el_clu2d->Write(0, TObject::kOverwrite);
  f125_pi_clu2d->Write(0, TObject::kOverwrite);
  mmg1_f125_el->Write(0, TObject::kOverwrite);
  mmg1_f125_pi->Write(0, TObject::kOverwrite);
  mmg1_f125_el_amp2d->Write(0, TObject::kOverwrite);
  mmg1_f125_pi_amp2d->Write(0, TObject::kOverwrite);
  mmg1_f125_el_clu2d->Write(0, TObject::kOverwrite);
  mmg1_f125_pi_clu2d->Write(0, TObject::kOverwrite);
  mmg2_f125_el->Write(0, TObject::kOverwrite);
  mmg2_f125_pi->Write(0, TObject::kOverwrite);
  mmg2_f125_el_amp2d->Write(0, TObject::kOverwrite);
  mmg2_f125_pi_amp2d->Write(0, TObject::kOverwrite);
  mmg2_f125_el_clu2d->Write(0, TObject::kOverwrite);
  mmg2_f125_pi_clu2d->Write(0, TObject::kOverwrite);
  urw_f125_el->Write(0, TObject::kOverwrite);
  urw_f125_pi->Write(0, TObject::kOverwrite);
  urw_f125_el_amp2d->Write(0, TObject::kOverwrite);
  urw_f125_pi_amp2d->Write(0, TObject::kOverwrite);
  urw_f125_el_clu2d->Write(0, TObject::kOverwrite);
  urw_f125_pi_clu2d->Write(0, TObject::kOverwrite);
  fOut->Write(0, TObject::kOverwrite);
	
  //close
  fOut->Close();
  delete fOut;

//=====================================================================================
//===                 P L O T  H I S T O G R A M S                                  ===
//=====================================================================================


  //---  global parameters --
  const char *OutputDir="FNAL_JANA2/RunOutput";
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
  //char pngname[120];  sprintf(pngname,"%s_evt.png","FNAL_JANA2/RunOutput/PNGs");  c0->Print(pngname);
  //char pdfname[120];  sprintf(pdfname,"%s_evt.pdf","FNAL_JANA2/RunOutput/PDFs");  c0->Print(pdfname);

  //---------------------  page 1 --------------------
  htitle(" Calorimeter ");   // if (!COMPACT) cc=NextPlot(0,0);

  cc=NextPlot(nxd,nyd);                   hCal_occ->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[6]->Draw();
  cc=NextPlot(nxd,nyd);                   hCal_sum->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[3]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[4]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[5]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[0]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[1]->Draw();
  cc=NextPlot(nxd,nyd);  gPad->SetLogy(); hCal_adc[2]->Draw();

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

  cc=NextPlot(nxd,nyd);  hcount->Draw();
  cc=NextPlot(nxd,nyd);  f125_el_chi2->Draw("colz");
  cc=NextPlot(nxd,nyd);  f125_pi_chi2->Draw("colz");

  //--- close PDF file ----
  cc=NextPlot(-1,-1);
  //--- the end ---
}
//===============================================================
