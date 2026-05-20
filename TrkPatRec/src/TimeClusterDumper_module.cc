// ============================================================
//  TimeClusterDumper_module.cc
//
//  art EDAnalyzer module that reads a mu2e::TimeClusterCollection
//  and writes every cluster's content into a ROOT TTree.
//
//  Author: B. Echenard
// ============================================================

// ---- art / fhicl headers ----
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ---- Mu2e data products ----
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"

// ---- ROOT headers ----
#include "TTree.h"
#include "TBranch.h"

// ---- STL ----
#include <string>
#include <vector>
#include <iostream>


namespace mu2e {

  class TimeClusterDumper : public art::EDAnalyzer {

  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   chTag {Name("ComboHitCollection"),    Comment("InputTag for ComboHit collection")};
      fhicl::Atom<art::InputTag>   tcTag {Name("TimeClusterCollection"), Comment("InputTag for TimeCluster collection")};
      fhicl::Atom<art::InputTag>   mcTag {Name("MCCollection"),          Comment("InputTag for MC collection")};
      fhicl::Sequence<std::string> hsel  {Name("HitSelectionBits"),      Comment("HitSelectionBits") };
      fhicl::Sequence<std::string> hbkg  {Name("HitBackgroundBits"),     Comment("HitBackgroundBits") };
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit TimeClusterDumper(Parameters const& conf);

    void beginJob()                            override;
    void analyze(art::Event const& evt)        override;


  private:
    art::InputTag _chTag;
    art::InputTag _tcTag;
    art::InputTag _mcTag;
    StrawHitFlag  _hsel;
    StrawHitFlag  _hbkg;


    TTree* _tree {nullptr};
    Int_t   _run,_subRun,_event;

    Int_t   _nClusters,_nsh[1024],_hasCalo[1024];
    Float_t _t0[1024],_t0Err[1024],_caloTime[1024],_caloEnergy[1024];
    std::vector<std::vector<Int_t>> _hitIdx;

    Int_t   Nch_,chSel_[8192],chNhit_[8192],chPdg_[8192],chCrCode_[8192],chSimId_[8192],chUId_[8192];
    Float_t chTime_[8192], chPhi_[8192],chRad_[8192],chX_[8192],chY_[8192],chZ_[8192],chDX_[8192],chDY_[8192],chDZ_[8192],chDQ_[8192];
    Float_t chTerr_[8192],chWerr_[8192],chWDX_[8192],chWDY_[8192],chMomX_[8192],chMomY_[8192],chMomZ_[8192],chEdep_[8192];

    void fillTree(art::Event const& evt,TimeClusterCollection const& tcc);
  };

  // ----------------------------------------------------------------
  TimeClusterDumper::TimeClusterDumper(Parameters const& conf)
    : art::EDAnalyzer(conf)
    , _chTag (conf().chTag())
    , _tcTag (conf().tcTag())
    , _mcTag (conf().mcTag())
    , _hsel  (conf().hsel())
    , _hbkg  (conf().hbkg())
  {}

  // ----------------------------------------------------------------
  void TimeClusterDumper::beginJob() {

    art::ServiceHandle<art::TFileService> tfs;

    _tree = tfs->make<TTree>("tcDump", "TimeCluster dump");

    // event-level branches
    _tree->Branch("run",         &_run,       "run/I");
    _tree->Branch("subRun",      &_subRun,    "subRun/I");
    _tree->Branch("event",       &_event,     "event/I");

    _tree->Branch("Nch",         &Nch_,       "Nch/I");
    _tree->Branch("chSel",       &chSel_,     "chSel[Nch]/I");
    _tree->Branch("chNhit",      &chNhit_,    "chNhit[Nch]/I");
    _tree->Branch("chPdg",       &chPdg_,     "chPdg[Nch]/I");
    _tree->Branch("chCrCode",    &chCrCode_,  "chCrCode[Nch]/I");
    _tree->Branch("chSimId",     &chSimId_,   "chSimId[Nch]/I");
    _tree->Branch("chTime",      &chTime_,    "chTime[Nch]/F");
    _tree->Branch("chPhi",       &chPhi_,     "chPhi[Nch]/F");
    _tree->Branch("chRad",       &chRad_,     "chRad[Nch]/F");
    _tree->Branch("chX",         &chX_,       "chX[Nch]/F");
    _tree->Branch("chY",         &chY_,       "chY[Nch]/F");
    _tree->Branch("chZ",         &chZ_,       "chZ[Nch]/F");
    _tree->Branch("chDX",        &chDX_,      "chDX[Nch]/F");
    _tree->Branch("chDY",        &chDY_,      "chDY[Nch]/F");
    _tree->Branch("chDZ",        &chDZ_,      "chDZ[Nch]/F");
    _tree->Branch("chDQ",        &chDQ_,      "chDQ[Nch]/F");
    _tree->Branch("chMomX",      &chMomX_,    "chMomX[Nch]/F");
    _tree->Branch("chMomY",      &chMomY_,    "chMomY[Nch]/F");
    _tree->Branch("chMomZ",      &chMomZ_,    "chMomZ[Nch]/F");
    _tree->Branch("chEdep",      &chEdep_,    "chEdep[Nch]/F");
    _tree->Branch("chUId",       &chUId_,     "chUId[Nch]/I");
    _tree->Branch("chTerr",      &chTerr_,    "chTerr[Nch]/F");
    _tree->Branch("chWerr",      &chWerr_,    "chWerr[Nch]/F");
    _tree->Branch("chWDX",       &chWDX_,     "chWDX[Nch]/F");
    _tree->Branch("chWDY",       &chWDY_,     "chWDY[Nch]/F");

    _tree->Branch("nClusters",   &_nClusters,   "nClusters/I");
    _tree->Branch("t0",          &_t0,          "t0[nClusters]/F");
    _tree->Branch("t0Err",       &_t0Err,       "t0Err[nClusters]/F");
    _tree->Branch("nsh",         &_nsh,         "nsh[nClusters]/I");
    _tree->Branch("hasCalo",     &_hasCalo,     "hasCalo[nClusters]/I");
    _tree->Branch("caloTime",    &_caloTime,    "caloTime[nClusters]/F");
    _tree->Branch("caloEnergy",  &_caloEnergy,  "caloEnergy[nClusters]/F");
    _tree->Branch("hitIdx",      &_hitIdx);

  }

  // ----------------------------------------------------------------
  void TimeClusterDumper::analyze(art::Event const& evt) {

    _run    = static_cast<Int_t>(evt.run());
    _subRun = static_cast<Int_t>(evt.subRun());
    _event  = static_cast<Int_t>(evt.event());

    // ---- retrieve collection ----
    auto const& tcHandle = evt.getValidHandle<TimeClusterCollection>(_tcTag);
    auto const& tcc = *tcHandle;

    auto const& chHandle = evt.getValidHandle<ComboHitCollection>(_chTag);
    auto const& chc = *chHandle;

    auto  mcdigisTag_(_mcTag);
    auto const& mcdigis = evt.getValidHandle<StrawDigiMCCollection>(_mcTag);


    _hitIdx.clear();
    _nClusters = static_cast<Int_t>(tcc.size());
    for (int ic = 0; ic < _nClusters; ++ic) {
      TimeCluster const& tc = tcc.at(static_cast<std::size_t>(ic));
      _t0[ic]     = static_cast<Float_t>(tc._t0.t0());
      _t0Err[ic]  = static_cast<Float_t>(tc._t0.t0Err());
      _nsh[ic]    = static_cast<Int_t>(tc._strawHitIdxs.size());

      std::vector<int> list;
      for (size_t ih = 0; ih < tc._strawHitIdxs.size(); ++ih) list.push_back(tc._strawHitIdxs.at(ih));
      _hitIdx.push_back(list);

      if (tc._caloCluster.isNonnull()) {
        _hasCalo[ic]    = 1;
        _caloTime[ic]   = static_cast<Float_t>(tc._caloCluster->time());
        _caloEnergy[ic] = static_cast<Float_t>(tc._caloCluster->energyDep());
      } else {
        _hasCalo[ic]    = 0;
        _caloTime[ic]   = 0.f;
        _caloEnergy[ic] = 0.f;
      }
    }

    Nch_= chc.size();
    for (unsigned ich=0; ich<chc.size();++ich)
    {
      const auto& ch = chc.at(ich);

      int selFlag  = ch.flag().hasAllProperties(_hsel) && !ch.flag().hasAnyProperty(_hbkg) ? 1 : 0;
      chSel_[ich]  = selFlag;
      chTime_[ich] = ch.correctedTime();
      chX_[ich]    = ch.pos().x();
      chY_[ich]    = ch.pos().y();
      chZ_[ich]    = ch.pos().z();
      chDX_[ich]   = ch.hDir().x();
      chDY_[ich]   = ch.hDir().y();
      chDZ_[ich]   = ch.hDir().z();
      chDQ_[ich]   = ch.qual();
      chRad_[ich]  = sqrt(ch.pos().x()*ch.pos().x()+ch.pos().y()*ch.pos().y());
      chPhi_[ich]  = ch.pos().phi();
      chNhit_[ich] = ch.nStrawHits();
      chUId_[ich]  = ch.strawId().uniquePanel();
      chTerr_[ich] = ch.posRes(StrawHitPosition::trans);
      chWerr_[ich] = ch.posRes(StrawHitPosition::wire);
      chWDX_[ich]  = ch.vDir().x();
      chWDY_[ich]  = ch.vDir().y();
      chEdep_[ich] = ch.energyDep();

      std::vector<StrawDigiIndex> dids;
      chc.fillStrawDigiIndices(ich,dids);
      const StrawDigiMC& mcdigi        = mcdigis->at(dids[0]);// taking 1st digi: is there a better idea??
      const art::Ptr<SimParticle>& spp = mcdigi.earlyStrawGasStep()->simParticle();
      chPdg_[ich]    = spp->pdgId();
      chCrCode_[ich] = spp->creationCode();
      chSimId_[ich]  = spp->id().asInt();
      chMomX_[ich]   = mcdigi.earlyStrawGasStep()->momentum().x();
      chMomY_[ich]   = mcdigi.earlyStrawGasStep()->momentum().y();
      chMomZ_[ich]   = mcdigi.earlyStrawGasStep()->momentum().z();
    }

    _tree->Fill();

  }

}

DEFINE_ART_MODULE(mu2e::TimeClusterDumper)
