//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Jun 17 22:57:15 2023 by ROOT version 6.14/04
// from TTree gem_hits/GEM TTree with single track hit info
// found on file: trd_singleTrackHits_Run_003200.root
//////////////////////////////////////////////////////////


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

   // Declaration of leaf types
   Int_t           event_num;
   Int_t           gem_nhit;
   vector<int>     *xpos;
   vector<float>   *zpos;
   vector<float>   *dedx;
   vector<bool>    *parID;
   vector<float>   *zHist;

   // List of branches
   TBranch        *b_event_num;   //!
   TBranch        *b_gem_nhit;   //!
   TBranch        *b_xpos;   //!
   TBranch        *b_zpos;   //!
   TBranch        *b_dedx;   //!
   TBranch        *b_parID;   //!
   TBranch        *b_zHist;   //!

