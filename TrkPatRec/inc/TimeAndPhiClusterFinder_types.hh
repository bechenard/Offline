#ifndef TimeAndPhiClusterFinderTypes_hh
#define TimeAndPhiClusterFinderTypes_hh

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "TTree.h"

namespace art {class Event;};

namespace mu2e {

  class ComboHitCollection;
  namespace TimeAndPhiClusterFinderTypes {

    struct Config{
        fhicl::Atom<std::string>   tool_type{             fhicl::Name("tool_type"),             fhicl::Comment("Needed fot backward compatibility"), "" };
        fhicl::Atom<bool>          mcDiag{                fhicl::Name("MCDiag"),                fhicl::Comment("Switch to perform MC diag"), true };
        fhicl::Atom<art::InputTag> strawDigiMCCollection{ fhicl::Name("StrawDigiMCCollection"), fhicl::Comment("StrawDigiMC collection name"), "makeSD" };
    };


    struct Data_t {

        enum  {kMaxHits = 8192};

        const art::Event*         event_;
        const ComboHitCollection* chcol_;

        Int_t   iev_;
        Int_t   Nch_,chSel_[kMaxHits],chNhit_[kMaxHits],chPdg_[kMaxHits],chCrCode_[kMaxHits],chSimId_[kMaxHits],chUId_[kMaxHits];
        Float_t chTime_[kMaxHits],chPhi_[kMaxHits],chRad_[kMaxHits],chX_[kMaxHits],chY_[kMaxHits],chZ_[kMaxHits];
        Float_t chTerr_[kMaxHits],chWerr_[kMaxHits],chWDX_[kMaxHits],chWDY_[kMaxHits];

        Int_t   nhit1_,iclu1_[kMaxHits],hitIdx1_[kMaxHits];
        Int_t   nhit2_,iclu2_[kMaxHits],hitIdx2_[kMaxHits];
        Int_t   nclu2_;
        Float_t calo2X_[kMaxHits],calo2Y_[kMaxHits],calo2Z_[kMaxHits],calo2T_[kMaxHits],calo2E_[kMaxHits];

        void  reset() {Nch_=nhit1_=nhit2_=nclu2_=0;}
    };
  }
}

#endif
