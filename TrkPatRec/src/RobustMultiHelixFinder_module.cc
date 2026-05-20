#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/polySinCos.hh"
#include "Offline/TrkPatRec/inc/RobustMultiHelixFinder_types.hh"

#include <numeric>
#include <stdint.h>

namespace {

  // ------------------------------------------------------------------
  // Pre-transformed calorimeter cluster position in tracker frame
  constexpr unsigned MAXCALOHIT = std::numeric_limits<unsigned int>::max();

  struct CaloCluHit {
     CaloCluHit(float x, float y, float z, float t, float e, float xyerr) :
       x_(x),y_(y),z_(z),t_(t),e_(e),xyerr2_(xyerr*xyerr) {}

     float x_,y_,z_,t_,e_,xyerr2_;
  };
  using CaloCluHits = std::vector<CaloCluHit>;

  // ------------------------------------------------------------------
  // Flat cache of ComboHit fields for vectorization-friendly access.
  struct HitFlat {
    std::vector<float>    x, y, z, t, edep, werr, terr, vdx, vdy, phiC;
    std::vector<unsigned> uid,nsh;

    size_t size() const {return x.size();}
  };

  // ------------------------------------------------------------------
  // Helix candidate structure
  struct CandHelix {
    explicit CandHelix(unsigned ccWeight) :
      nStrawHits_(0), ccWeight_(ccWeight), caloIdx_(MAXCALOHIT), x_(0), y_(0), r_(0), dpdz_(0), phi0_(0),
      dzdt_(0), t0_(0), dzdtErr_(0), t0Err_(0), chisq_(9999)
    {hits_.reserve(32);}

    unsigned nhitsTot()const {return caloIdx_ < MAXCALOHIT ? nStrawHits_+ccWeight_ : nStrawHits_;}
    float    pOffset() const {return dpdz_ > 0 ? -3.14159f : 3.14159f;}

    void clear() {
      hits_.clear();
      nStrawHits_ = 0;
      x_ = y_ = r_ = dpdz_ = phi0_ = dzdt_ = t0_ = dzdtErr_ = t0Err_ =0.0f;
      chisq_ = 9999.0f;
      caloIdx_ = MAXCALOHIT;
    }

    unsigned            nStrawHits_, ccWeight_, caloIdx_;
    float               x_, y_, r_, dpdz_, phi0_, dzdt_, t0_, dzdtErr_, t0Err_, chisq_;
    std::vector<size_t> hits_;
  };

  // ------------------------------------------------------------------
  class LSFitter {
    public:
      LSFitter() : sn_(0),sx_(0),sx2_(0),sy_(0),sxy_(0){};

      void  add(float x, float y, float w=1.0) {sx_+=x*w; sx2_+=x*x*w; sy_+=y*w; sxy_+=x*y*w; sn_+=w;}
      void  clear()                            {sn_=sx_=sx2_=sy_=sxy_=0.0;}

      float fa()   {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sn_*sxy_-sx_*sy_) /(sn_*sx2_-sx_*sx_) : 0.0;}
      float fb()   {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sy_*sx2_-sx_*sxy_)/(sn_*sx2_-sx_*sx_) : 0.0;}
      float xbar() {return sx_/sn_;}

    private:
      float sn_,sx_,sx2_,sy_,sxy_;
  };

}



namespace mu2e {

  class RobustMultiHelixFinder : public art::EDProducer
  {
    public:
      using mapHelix        = std::map<Helicity,std::unique_ptr<HelixSeedCollection>>;
      using strawHitIndices = std::vector<StrawHitIndex>;
      using Config_types    = RobustMultiHelixFinderTypes::Config ;
      using Data_types      = RobustMultiHelixFinderTypes::Data_t;

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> comboHitCollection    {Name("ComboHitCollection"),     Comment("ComboHit collection Name")    };
        fhicl::Atom<art::InputTag> timeClusterCollection {Name("TimeClusterCollection"),  Comment("TimeCluster collection Name") };
        fhicl::Atom<art::InputTag> caloClusterCollection {Name("CaloClusterCollection"),  Comment("CaloCluster collection Name") };
        fhicl::Sequence<int>       helicities            {Name("Helicities"),             Comment("Helicity values") };
        fhicl::Atom<float>         maxEdepHit            {Name("MaxEdepHit"),             Comment("Maximum eDep for hit - circle stage") };
        fhicl::Atom<unsigned>      minHitStride          {Name("MinHitStride"),           Comment("Minimum number of hits before skipping triplets") };
        fhicl::Atom<float>         minRadCircle          {Name("MinRadCircle"),           Comment("Minimum circle radius") };
        fhicl::Atom<float>         maxRadCircle          {Name("MaxRadCircle"),           Comment("Maximum circle radius") };
        fhicl::Atom<float>         minDZCircle           {Name("MinDZCircle"),            Comment("Minimum Z distance between hits for circle fit") };
        fhicl::Atom<float>         minDXY2Circle         {Name("MinDXY2Circle"),          Comment("Minimum XY distance between hits for circle fit") };
        fhicl::Atom<float>         maxDXY2Circle         {Name("MaxDXY2Circle"),          Comment("Maximum XY distance between hits for circle fit") };
        fhicl::Atom<float>         maxDRCircle           {Name("MaxDRCircle"),            Comment("Maximum radial distance between circle and hit - circle stage") };
        fhicl::Atom<float>         resCRad               {Name("ResCRad"),                Comment("Radial resolution of track circle") };
        fhicl::Atom<float>         resCPerp              {Name("ResCPerp"),               Comment("Perpendicular resolution of track circle") };
        fhicl::Atom<float>         maxDphi               {Name("MaxDphi"),                Comment("Maximum phi difference between hit and dz/dphi line init fit") };
        fhicl::Atom<float>         maxDphiC              {Name("MaxDphiC"),               Comment("Maximum phiC difference between consecutive hits") };
        fhicl::Atom<float>         minDZTrk              {Name("MinDZTrk"),               Comment("Minimum DZ track length") };
        fhicl::Atom<float>         maxDZGap              {Name("MaxDZGap"),               Comment("Maximum DZ gap between hits") };
        fhicl::Atom<float>         selHitRatio           {Name("SelHitRatio"),            Comment("Ratio of hits selected after dpdz fit") };
        fhicl::Atom<float>         minDpdz               {Name("MinDpdz"),                Comment("Minimum dp/dz helix candidate") };
        fhicl::Atom<float>         maxDpdz               {Name("MaxDpdz"),                Comment("Maximum dp/dz helix candidate") };
        fhicl::Atom<float>         maxChi2Filter         {Name("MaxChi2Filter"),          Comment("Maximum chi2 to filter hits") };
        fhicl::Atom<float>         maxChi2Recover        {Name("MaxChi2Recover"),         Comment("Maximum chi2 to recover hits") };
        fhicl::Atom<unsigned>      minStrawHits          {Name("MinStrawHits"),           Comment("Minimum number of Straw hits for a helix candidate") };
        fhicl::Atom<unsigned>      minnTotHits           {Name("MinnTotHits"),            Comment("Minimum number of Straw+Calo hits for a helix candidate") };
        fhicl::Atom<float>         ccMinEnergy           {Name("CaloClusterMinE"),        Comment("Minimum calo cluster energy") };
        fhicl::Atom<unsigned>      ccWeight              {Name("CaloClusterWeight"),      Comment("Calo cluster weight") };
        fhicl::Atom<float>         ccmaxDT               {Name("CaloClustermaxDT"),       Comment("Maximum time difference calo time and hits time") };
        fhicl::Atom<float>         ccmaxChi2             {Name("CaloClusterChi2"),        Comment("Maximum chi2 for calo hit to be associated to trk") };
        fhicl::Atom<float>         ccXYerr               {Name("CaloClusterXYerr"),       Comment("XY position error of calo cluster") };
        fhicl::Atom<float>         ccRmin                {Name("CaloClusterRmin"),        Comment("Effective calorimeter inner radius") };
        fhicl::Atom<unsigned>      nMaxTrkIter           {Name("NMaxTrkIter"),            Comment("Number of track finding iterations") };
        fhicl::Atom<int>           diagLevel             {Name("DiagLevel"),              Comment("Diag level"), 0 };
        fhicl::Table<Config_types> diagPlugin            {Name("DiagPlugin"),             Comment("Diag Plugin config")};
      };
      explicit RobustMultiHelixFinder(const art::EDProducer::Table<Config>& config);
      virtual void produce(art::Event& event);
      virtual void beginJob();


    private:
      const art::ProductToken<ComboHitCollection>    chToken_;
      const art::ProductToken<TimeClusterCollection> tcToken_;
      const art::ProductToken<CaloClusterCollection> ccToken_;
      float                                          maxEdepHit_;
      unsigned                                       minHitStride_;
      float                                          minRadCircle_;
      float                                          maxRadCircle_;
      float                                          minDZCircle_;
      float                                          minDXY2Circle_;
      float                                          maxDXY2Circle_;
      float                                          maxDRCircle_;
      float                                          resCRad_;
      float                                          resCPerp_;
      float                                          minDpdz_;
      float                                          maxDpdz_;
      float                                          maxDphi_;
      float                                          maxDphiC_;
      float                                          minDZTrk_;
      float                                          maxDZGap_;
      float                                          selHitRatio_;
      float                                          maxChi2Filter_;
      float                                          maxChi2Recover_;
      unsigned                                       minStrawHits_;
      unsigned                                       minnTotHits_;
      float                                          ccMinEnergy_;
      unsigned                                       ccWeight_;
      float                                          ccmaxDT_;
      float                                          ccmaxChi2_;
      float                                          ccXYerr_;
      float                                          ccRmin2_;
      unsigned                                       nMaxTrkIter_;
      const Calorimeter*                             cal_;
      int                                            iev_;
      int                                            diag_;
      std::unique_ptr<ModuleHistToolBase>            diagTool_;
      Data_types                                     data_;
      std::vector<Helicity>                          hels_;



      void      findAllHelices        (art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
                                       const art::ValidHandle<ComboHitCollection>& chH, const art::ValidHandle<CaloClusterCollection>& ccH);
      void      filterDuplicateHelices(HelixSeedCollection& helices);
      void      findHelicesInTC       (mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc,
                                       const ComboHitCollection& chcol, const CaloClusterCollection& cccol);

      CandHelix findHelixCandidate    (const ComboHitCollection& chcol, const strawHitIndices& hits, const CaloCluHits& caloHits);
      void      assignHitsCircle      (const HitFlat& hitsf,            CandHelix& helix,            const CaloCluHits& caloHits);
      void      checkCaloHitChi2      (const CaloCluHits& caloHits,     CandHelix& helix);

      //std::vector<float> chi2XYCircleFast(const HitFlat& hitsf, float centerX, float centerY, float radius);
      std::vector<float> chi2Helix    (const HitFlat& hitsf, const CandHelix& helix);

      //void      initDpdzScan          (const HitFlat& hitsf, CandHelix& helix, const std::vector<size_t>& hitsScan);
      void      initDpdzTriplet       (const HitFlat& hitsf, CandHelix& helix, size_t is, size_t js, size_t ks);
      void      fillPhiZFitter        (const HitFlat& hitsf, const CandHelix& helix, LSFitter& zphiFitter);
      void      cleanDpdzHits         (const HitFlat& hitsf, CandHelix& helix);
      bool      checkHelixQual        (const HitFlat& hitsf, const CandHelix& helix, int nshCircle);
      void      calcTimeProp          (const HitFlat& hitsf, CandHelix& helix);

      void      helixFilterChi2       (const HitFlat& hitsf, CandHelix& helix);
      void      fitCircle             (const HitFlat& hitsf, CandHelix& helix);
      void      refineHelix              (      HitFlat& hitsf, CandHelix& helix);
      void      helixRecoverChi2      (const HitFlat& hitsf, CandHelix& helix);
      float     helixCalcChi2         (const HitFlat& hitsf, CandHelix& helix);

      void      printHelix            (const HelixSeed& helix);
      void      fillDiag              (mapHelix& helcols, const ComboHitCollection& chcol);
  };


  RobustMultiHelixFinder::RobustMultiHelixFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer   {config},
    chToken_          {consumes<ComboHitCollection>   (config().comboHitCollection())},
    tcToken_          {consumes<TimeClusterCollection>(config().timeClusterCollection())},
    ccToken_          {consumes<CaloClusterCollection>(config().caloClusterCollection())},
    maxEdepHit_       (config().maxEdepHit()),
    minHitStride_     (config().minHitStride()),
    minRadCircle_     (config().minRadCircle()),
    maxRadCircle_     (config().maxRadCircle()),
    minDZCircle_      (config().minDZCircle()),
    minDXY2Circle_    (config().minDXY2Circle()),
    maxDXY2Circle_    (config().maxDXY2Circle()),
    maxDRCircle_      (config().maxDRCircle()),
    resCRad_          (config().resCRad()),
    resCPerp_         (config().resCPerp()),
    minDpdz_          (config().minDpdz()),
    maxDpdz_          (config().maxDpdz()),
    maxDphi_          (config().maxDphi()),
    maxDphiC_         (config().maxDphiC()),
    minDZTrk_         (config().minDZTrk()),
    maxDZGap_         (config().maxDZGap()),
    selHitRatio_      (config().selHitRatio()),
    maxChi2Filter_    (config().maxChi2Filter()),
    maxChi2Recover_   (config().maxChi2Recover()),
    minStrawHits_     (std::max(config().minStrawHits(),3u)),
    minnTotHits_      (std::max(config().minnTotHits(),3u)),
    ccMinEnergy_      (config().ccMinEnergy()),
    ccWeight_         (config().ccWeight()),
    ccmaxDT_          (config().ccmaxDT()),
    ccmaxChi2_        (config().ccmaxChi2()),
    ccXYerr_          (config().ccXYerr()),
    ccRmin2_          (config().ccRmin()*config().ccRmin()),
    nMaxTrkIter_      (config().nMaxTrkIter()),
    iev_              (0),
    diag_             (config().diagLevel()),
    diagTool_(),
    data_()
  {
    std::vector<int> helvals = config().helicities();
    for (const auto& hv : helvals) {
      Helicity hel(hv);
      hels_.emplace_back(hel);
      produces<HelixSeedCollection>(Helicity::name(hel));
    }
    if (diag_) diagTool_ = art::make_tool<ModuleHistToolBase>(config().diagPlugin," ");
  }


  //--------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::beginJob(){
    if (diag_){
       art::ServiceHandle<art::TFileService> tfs;
       diagTool_->bookHistograms(tfs);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::produce(art::Event& event )
  {
    if (diag_>0) std::cout<<"Event "<<event.id().event()<<std::endl;
    iev_ = event.id().event();

    const auto& tcH = event.getValidHandle(tcToken_);
    const auto& chH = event.getValidHandle(chToken_);
    const auto& ccH = event.getValidHandle(ccToken_);

    mapHelix helcols;
    for (const auto& hel : hels_) helcols[hel] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
    findAllHelices(event, helcols,tcH, chH, ccH);
    for (const auto& hel : hels_) event.put(std::move(helcols[hel]),Helicity::name(hel));
  }




  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in an event (hopefully), some timeCluster implementations have overlapping content, so filter duplicates
  void RobustMultiHelixFinder::findAllHelices(art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
                                              const art::ValidHandle<ComboHitCollection>& chH,
                                              const art::ValidHandle<CaloClusterCollection>& ccH)
  {
    cal_ = &(*GeomHandle<Calorimeter>());

    const TimeClusterCollection& tccol(*tcH);
    const ComboHitCollection&    chcol(*chH);
    const CaloClusterCollection& cccol(*ccH);

    if (diag_) {data_.reset(); data_.event_=&event; data_.chcol_ = &chcol;}

    for (size_t index=0;index<tccol.size();++index) {
      const auto tcArtPtr = art::Ptr<TimeCluster>(tcH,index);
      const auto& tc = tccol[index];
      findHelicesInTC(helcols,tcArtPtr,tc,chcol,cccol);
    }
    filterDuplicateHelices(*helcols[Helicity::poshel]);
    filterDuplicateHelices(*helcols[Helicity::neghel]);

    if (diag_) {
      fillDiag(helcols,chcol);
      diagTool_->fillHistograms(&data_);
      if (diag_>2) for (const auto& hel : *helcols[Helicity::poshel]) printHelix(hel);
      if (diag_>2) for (const auto& hel : *helcols[Helicity::neghel]) printHelix(hel);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Filter helices which have the same set of hits in common
  void RobustMultiHelixFinder::filterDuplicateHelices(HelixSeedCollection& helices)
  {
    for (auto first = helices.begin(); first != helices.end(); ++first)
    {
      if (first->hits().empty()) continue;
      for (auto second = std::next(first); second != helices.end(); ++second)
      {
        if (second->hits().empty()) continue;

        size_t result = 0;
        auto helix1 = first->hits();
        auto helix2 = second->hits();
        auto iter1  = helix1.begin();
        auto iter2  = helix2.begin();
        while (iter1 != helix1.end() && iter2 != helix2.end()) {
          if      (iter1->index(0) < iter2->index(0)) ++iter1;
          else if (iter2->index(0) < iter1->index(0)) ++iter2;
          else { ++result; ++iter1; ++iter2; }
        }

        if (result != helix1.size() || result != helix2.size()) continue;

        if (first->hits().size() >= second->hits().size()) second->_hhits.clear();
        else                                               first->_hhits.clear();
      }
    }

    auto predEmpty = [](const HelixSeed& hel) { return hel.hits().empty(); };
    helices.erase(std::remove_if(helices.begin(), helices.end(), predEmpty), helices.end());
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in a timeCluster for both helicities. Can find up to nMaxTrkIter helices for each helicity
  void RobustMultiHelixFinder::findHelicesInTC(mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc,
                                                const ComboHitCollection& chcol, const CaloClusterCollection& cccol)
  {
    //Select calo clusters compatible with time cluster
    CaloCluHits caloHits;
    float timeCluster(tc.t0().t0());
    for (const auto& ccalo : cccol)
    {
       if (fabs(ccalo.time() - timeCluster) > ccmaxDT_) continue;
       auto ccPosTrk = cal_->geomUtil().mu2eToTracker(cal_->geomUtil().diskFFToMu2e(ccalo.diskID(),ccalo.cog3Vector()));
       caloHits.emplace_back(CaloCluHit(ccPosTrk.x(),ccPosTrk.y(),ccPosTrk.z(),ccalo.time(),ccalo.energyDep(),ccXYerr_));
    }

    auto hitsToProcess = tc.hits();
    auto pred = [&chcol](const auto& i, const auto& j) {return chcol[i].strawId().uniquePanel()<chcol[j].strawId().uniquePanel();};
    sort(hitsToProcess.begin(),hitsToProcess.end(),pred);

    for (size_t it=0;it<nMaxTrkIter_;++it)
    {
      unsigned strawHitsToGo = std::accumulate(hitsToProcess.begin(), hitsToProcess.end(), 0,
                               [&](unsigned acc, auto ich){return acc + chcol[ich].nStrawHits();});

      if (strawHitsToGo < minStrawHits_) break;

      CandHelix bestHelix = findHelixCandidate(chcol, hitsToProcess, caloHits);
      if (bestHelix.nStrawHits_ < minStrawHits_) break;

      // Remove elements from HitsToProcess if they are found in
      auto predFound = [&bestHelix](auto x) {return std::find(bestHelix.hits_.begin(),
                                             bestHelix.hits_.end(), x) != bestHelix.hits_.end();};
      hitsToProcess.erase(std::remove_if(hitsToProcess.begin(), hitsToProcess.end(), predFound),hitsToProcess.end());

      //Save helix in data structure
      Helicity bestHelicity = bestHelix.dpdz_ > 0 ? Helicity::poshel: Helicity::neghel;
      float Rcent  = sqrt(bestHelix.x_*bestHelix.x_+bestHelix.y_*bestHelix.y_);
      float Fcent  = polyAtan2(bestHelix.y_,bestHelix.x_);


      HelixSeed hseed;
      hseed._helix = RobustHelix(Rcent,Fcent,bestHelix.r_,1.0/bestHelix.dpdz_,bestHelix.phi0_);
      hseed._hhits.setParent(chcol.parent());
      hseed._helix._helicity  = bestHelicity;
      hseed._helix._chi2dXY   = bestHelix.dzdt_;  // Dirty hack to save particle propagation direction
      hseed._helix._chi2dZPhi = bestHelix.chisq_;
      hseed._t0 = TrkT0(bestHelix.t0_,bestHelix.t0Err_);
      for (const auto& ich : bestHelix.hits_) hseed._hhits.emplace_back(chcol[ich]);
      hseed._status.merge(TrkFitFlag::MPRHelix);
      hseed._status.merge(TrkFitFlag::helixOK);
      hseed._timeCluster = tcArtPtr;

      helcols[bestHelicity]->push_back(std::move(hseed));
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Loop over triplets to find helices. Find circle first, then look at full helix
  CandHelix RobustMultiHelixFinder::findHelixCandidate(const ComboHitCollection& chcol, const strawHitIndices& hitsInit, const CaloCluHits& caloHits)
  {
    //-- Start loop over all hit triplets
    CandHelix helix(ccWeight_), bestHelix(ccWeight_);

    //Flag proton tracks here
    /*
    int nFlatStr(0),nFlagPRot(0);
    for (auto idx : hitsInit) {
      const auto& ch = chcol[idx];
      nFlatStr += ch.nStrawHits();
      if (ch.energyDep()>maxEdepHit_) nFlagPRot += ch.nStrawHits();
    }
    if (float(nFlagPRot)/float(nFlatStr) > 0.5) return helix;
    */

    //Flat hit structure to improve vectorization
    HitFlat hitsf;
    for (auto idx : hitsInit) {
      const auto& ch = chcol[idx];
      hitsf.x.emplace_back(ch.pos().x());
      hitsf.y.emplace_back(ch.pos().y());
      hitsf.z.emplace_back(ch.pos().z());
      hitsf.t.emplace_back(ch.correctedTime());
      hitsf.edep.emplace_back(ch.energyDep());
      hitsf.werr.emplace_back(ch.posRes(StrawHitPosition::wire));
      hitsf.terr.emplace_back(ch.posRes(StrawHitPosition::trans));
      hitsf.vdx.emplace_back(ch.vDir().x());
      hitsf.vdy.emplace_back(ch.vDir().y());
      hitsf.phiC.emplace_back(0.0);
      hitsf.uid.emplace_back(ch.strawId().uniquePanel());
      hitsf.nsh.emplace_back(ch.nStrawHits());
    }

    size_t kstride = hitsInit.size() < minHitStride_ ? 1 : 2;
    for (size_t i=0; i+2<hitsf.size(); ++i)
    {
      float rad2x1y1  = hitsf.x[i]*hitsf.x[i] + hitsf.y[i]*hitsf.y[i];
      if (hitsf.edep[i] > maxEdepHit_) continue;

      for (size_t j=i+1; j+1<hitsf.size(); ++j)
      {
        float x2        = hitsf.x[j] - hitsf.x[i];
        float y2        = hitsf.y[j] - hitsf.y[i];
        float d12       = x2*x2 + y2*y2;
        float rad2x2y2  = hitsf.x[j]*hitsf.x[j] + hitsf.y[j]*hitsf.y[j];
        if (hitsf.edep[j] > maxEdepHit_)                  continue;
        if (fabs(hitsf.z[j] - hitsf.z[i]) < minDZCircle_) continue;
        if (d12 < minDXY2Circle_ || d12 > maxDXY2Circle_) continue;

        for (size_t k=j+1; k<hitsf.size(); k+=kstride)
        {
          float x3        = hitsf.x[k] - hitsf.x[i];
          float y3        = hitsf.y[k] - hitsf.y[i];
          float d13       = x3*x3+y3*y3;
          float rad2x3y3  = hitsf.x[k]*hitsf.x[k]+hitsf.y[k]*hitsf.y[k];

          if (hitsf.edep[k] > maxEdepHit_)                  continue;
          if (fabs(hitsf.z[k] - hitsf.z[j]) < minDZCircle_) continue;
          if (d13 < minDXY2Circle_ || d13 > maxDXY2Circle_) continue;

          float x4  = hitsf.x[k] - x2;
          float y4  = hitsf.y[k] - y2;
          float d23 = x4*x4 + y4*y4;
          if (d23 < minDXY2Circle_ || d23 > maxDXY2Circle_) continue;

          float denom      = 2.0f*(x2*y3 - x3*y2);
          float numX       = d12*y3 - d13*y2;
          float numY       = d13*x2 - d12*x3;
          float centerX    = numX/denom + hitsf.x[i];
          float centerY    = numY/denom + hitsf.y[i];
          float radCenter2 = centerX*centerX + centerY*centerY;

          float dcx        = hitsf.x[i] - centerX;
          float dcy        = hitsf.y[i] - centerY;
          float radius     = sqrtf(dcx*dcx + dcy*dcy);
          if (radius < minRadCircle_ || radius > maxRadCircle_ ) continue;
          if (radCenter2 > rad2x1y1  || radCenter2 > rad2x2y2 || radCenter2 > rad2x3y3) continue;

          //-- Set the helix candidate
          helix.clear();
          helix.x_ = centerX;
          helix.y_ = centerY;
          helix.r_ = radius;

          //-- Select all hits close to the circle
          assignHitsCircle(hitsf, helix, caloHits);
          if (helix.nhitsTot() < minnTotHits_) continue;
          auto nshCircle = helix.nStrawHits_;

          //-- Estimate dp/dz - slower but possibly more efficient
          //for (size_t n = 0; n < hitsf.size(); ++n) {
          //  hitsf.phiC[n] = polyAtan2(hitsf.y[n] - helix.y_, hitsf.x[n] - helix.x_);
          //}
          //std::vector<size_t> hitsToScan{i,j,k};
          //initDpdzScan(hitsf, helix, hitsToScan);
          //if (fabs(helix.dpdz_) < minDpdz_ || fabs(helix.dpdz_) > maxDpdz_) continue;

          //-- Estimate dp/dz and recalculate phiC for hits
          initDpdzTriplet(hitsf,helix,i,j,k);
          if (fabs(helix.dpdz_) < minDpdz_ || fabs(helix.dpdz_) > maxDpdz_) continue;
          for (size_t n = 0; n < hitsf.size(); ++n) {
            hitsf.phiC[n] = polyAtan2(hitsf.y[n] - helix.y_, hitsf.x[n] - helix.x_);
          }

          //-- Remove incompatible hits
          cleanDpdzHits(hitsf,helix);
          if (helix.nStrawHits_ < minStrawHits_) continue;

          //-- Check helix quality after initial helix reco
          bool passHelix = checkHelixQual(hitsf, helix, nshCircle);
          if (!passHelix) continue;

          //-- Filter hits
          helixFilterChi2(hitsf,helix);

          //-- Refine helix xy position and recalculate dpdz/phi0
          fitCircle(hitsf,helix);
          if (helix.r_ < minRadCircle_ || helix.r_ > maxRadCircle_ ) continue;

          for (size_t n = 0; n < hitsf.size(); ++n) {
            hitsf.phiC[n] = atan2f(hitsf.y[n] - helix.y_, hitsf.x[n] - helix.x_);
          }
          LSFitter zphiFitter;
          fillPhiZFitter(hitsf,helix,zphiFitter);
          helix.dpdz_ = zphiFitter.fa();
          helix.phi0_ = zphiFitter.fb();
          if (fabs(helix.dpdz_) < minDpdz_ || fabs(helix.dpdz_) > maxDpdz_) continue;

          //-- Find all compatible hits and calculate dzdt and t0
          helixRecoverChi2(hitsf,helix);
          calcTimeProp(hitsf,helix);

          //-- Recheck calo cluster agreement in space and time
          checkCaloHitChi2(caloHits, helix);

          //-- Select best helix
          if (helix.nhitsTot() > bestHelix.nhitsTot() || (helix.nhitsTot() == bestHelix.nhitsTot() && helix.chisq_<bestHelix.chisq_) ) {
            refineHelix(hitsf,helix);
            helix.chisq_ = helixCalcChi2(hitsf,helix);
            bestHelix = helix;

            //One could potentially check for a early stopping condition for protons, do not forget to add
            //for (auto& val : bestHelix.hits_) val = hitsInit[val];
          }
        }
      }
    }

    if (bestHelix.nhitsTot() < minnTotHits_ ) {bestHelix.clear(); return bestHelix;}

    //-- Do not forget to restore the original indices!
    for (auto& val : bestHelix.hits_) val = hitsInit[val];

    return bestHelix;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  //Vectorization friendly assignation of hits to the helix
  void RobustMultiHelixFinder::assignHitsCircle(const HitFlat& hitsf, CandHelix& helix, const CaloCluHits& caloHits)
  {
     float invDiameter = 0.5f / helix.r_ ;
     float radius2     = helix.r_*helix.r_ ;

     std::vector<float> val(hitsf.size());
     for (size_t n = 0; n < hitsf.size(); ++n) {
         float dx = hitsf.x[n] - helix.x_;
         float dy = hitsf.y[n] - helix.y_;
         val[n] = fabsf((dx*dx + dy*dy - radius2)*invDiameter);
     }

     float previousValue = 1e6f;
     for (size_t n = 0; n < hitsf.size(); ++n) {
        if (val[n] > maxDRCircle_) continue;

        if (helix.hits_.empty() || hitsf.uid[n] != hitsf.uid[helix.hits_.back()]) {
            previousValue      = val[n];
            helix.nStrawHits_ += hitsf.nsh[n];
            helix.hits_.emplace_back(n);
        } else if (val[n] < previousValue) {
            previousValue       = val[n];
            helix.nStrawHits_  += hitsf.nsh[n] - hitsf.nsh[helix.hits_.back()];
            helix.hits_.back()  = n;
        }
     }

     //-- Add the calo cluster based on radial distance with circle
     //   take most energetic cluster if there are many
     float caloE(0);
     for (size_t icalo=0;icalo<caloHits.size();++icalo) {
       const auto& caloHit = caloHits[icalo];
       float dx = caloHit.x_ - helix.x_;
       float dy = caloHit.y_ - helix.y_;
       float dr = fabs((dx*dx+dy*dy-radius2)*invDiameter);
       if (dr < maxDRCircle_ && caloHit.e_ > caloE)
       {
          caloE = caloHit.e_;
          helix.caloIdx_ = icalo;
       }
     }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::checkCaloHitChi2(const CaloCluHits& caloHits, CandHelix& helix)
  {
    const int   nMAXITER(10);

    if (helix.caloIdx_==MAXCALOHIT) return;
    const auto& caloHit = caloHits[helix.caloIdx_];

    float deltaTime = helix.dzdt_*caloHits[helix.caloIdx_].z_ + helix.t0_ - caloHits[helix.caloIdx_].t_;
    if (fabs( helix.dzdt_) > 1e-3 && fabs(deltaTime)>ccmaxDT_)  return;

    float phiAtZ = caloHit.z_*helix.dpdz_ + helix.phi0_;
    float s, c;
    polySinCos(phiAtZ, s, c);
    float cx     = helix.x_ + helix.r_*c;
    float cy     = helix.y_ + helix.r_*s;
    float cr2    = cx*cx+cy*cy;

    int niter(0);
    while (cr2 < ccRmin2_ && niter<nMAXITER)
    {
      phiAtZ += 20*helix.dpdz_;
      polySinCos(phiAtZ, s, c);
      cx     = helix.x_ + helix.r_*c;
      cy     = helix.y_ + helix.r_*s;
      cr2    = cx*cx+cy*cy;
      ++niter;
    }

    float dhx  = caloHit.x_ - (helix.x_ + helix.r_*c);
    float dhy  = caloHit.y_ - (helix.y_ + helix.r_*s);
    float chi2 = dhx*dhx/caloHit.xyerr2_ + dhy*dhy/caloHit.xyerr2_;
    if (chi2 > ccmaxChi2_) helix.caloIdx_ = MAXCALOHIT;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 based on circle fit - this might be needed in the future, so keep it
  /*
  float RobustMultiHelixFinder::chi2XYCircleFast(const HitFlat& hitsf, float centerX, float centerY, float radius)
  {
    std::vector<float> chisqvec;
    chisqvec.reserve(hitsf.size());

    for (size_t k=0;k<hitsf.size();++k)
    {
      float dx     = hitsf.x[k] - centerX;
      float dy     = hitsf.y[k] - centerY;
      float r      = sqrtf(dx*dx+dy*dy);
      float rwdot  = (hitsf.vdx[k]*dx + hitsf.vdy[k]*dy)/r;
      float dr     = r - radius;
      float rwdot2 = rwdot*rwdot;
      float werr   = hitsf.werr[k];
      float terr   = hitsf.terr[k];
      float rres2  = werr*werr*rwdot2 + terr*terr*(1.0-rwdot2);
      float chisq  = dr*dr/rres2;
      chisqvec.emplace_back(chisq);
    }
    return chisqvec;
  */

  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 with full helix hypothesis - vectorization friendly
  std::vector<float> RobustMultiHelixFinder::chi2Helix(const HitFlat& hitsf, const CandHelix& helix)
  {
    const size_t N = hitsf.size();
    std::vector<float> chisqvec(N);

    const float cx       = helix.x_;
    const float cy       = helix.y_;
    const float r        = helix.r_;
    const float dpdz     = helix.dpdz_;
    const float phi0     = helix.phi0_;
    const float cradres  = resCRad_  * r;
    const float cperpres = resCPerp_ * r;
    const float crad2    = cradres  * cradres;
    const float cperp2   = cperpres * cperpres;

    // Hoisting loop-invariants and using pointer-based access removes
    // repeated bounds-checking and exposes a simple stride-1 pattern.
    const float* __restrict__ xp    = hitsf.x.data();
    const float* __restrict__ yp    = hitsf.y.data();
    const float* __restrict__ zp    = hitsf.z.data();
    const float* __restrict__ vdxp  = hitsf.vdx.data();
    const float* __restrict__ vdyp  = hitsf.vdy.data();
    const float* __restrict__ werrp = hitsf.werr.data();
    const float* __restrict__ terrp = hitsf.terr.data();
    float*       __restrict__ outp  = chisqvec.data();

    for (size_t k = 0; k < N; ++k)
    {
      const float dx   = xp[k] - cx;
      const float dy   = yp[k] - cy;
      const float r2   = dx*dx + dy*dy;
      const float cdd2 =  dx*vdxp[k] + dy*vdyp[k];
      const float cpd2 = -dy*vdxp[k] + dx*vdyp[k];

      float s, c;
      polySinCos(zp[k]*dpdz + phi0, s, c);

      const float dhx    = dx - r*c;
      const float dhy    = dy - r*s;
      const float dwire  =  vdxp[k]*dhx + vdyp[k]*dhy;
      const float dtrans = -vdyp[k]*dhx + vdxp[k]*dhy;

      const float r2inv  = 1.0f / r2;
      const float wres2  = werrp[k]*werrp[k] + (crad2*cdd2*cdd2 + cperp2*cpd2*cpd2)*r2inv;
      const float tres2  = terrp[k]*terrp[k] + (crad2*cpd2*cpd2 + cperp2*cdd2*cdd2)*r2inv;

      outp[k] = dwire*dwire/wres2 + dtrans*dtrans/tres2;
    }
    return chisqvec;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::initDpdzTriplet(const HitFlat& hitsf, CandHelix& helix, size_t is, size_t js, size_t ks)
  {
    const float zi   = hitsf.z[is];
    const float zj   = hitsf.z[js];
    const float zk   = hitsf.z[ks];
    const float phii = polyAtan2(hitsf.y[is] - helix.y_, hitsf.x[is] - helix.x_);
    const float phij = polyAtan2(hitsf.y[js] - helix.y_, hitsf.x[js] - helix.x_);
    const float phik = polyAtan2(hitsf.y[ks] - helix.y_, hitsf.x[ks] - helix.x_);

    int ilmax  = int((zj-zi)*maxDpdz_/6.2831f - (phij-phii)/6.2831f);
    int ilmin  = int((zi-zj)*maxDpdz_/6.2831f - (phij-phii)/6.2831f);

    float mindPhi(maxDphi_);
    for (int il=ilmin; il<=ilmax; ++il)
    {
      float dpdz = (phij-phii + il*6.2831f)/(zj-zi);
      if (fabs(dpdz) > maxDpdz_ || fabs(dpdz)<minDpdz_) continue;

      float pOffset = dpdz > 0 ? -3.1415f : 3.1415f;
      float phi0    = phii - zi*dpdz;
      float n       = std::trunc((zk*dpdz + phi0 - pOffset)/6.2831f);
      float dPhi    = fabs(zk*dpdz + phi0 - n*6.2831f - phik);
      if (dPhi > mindPhi) continue;
      mindPhi = dPhi;

      float pi    = phii;
      float pj    = phij + il*6.2831f;
      float pk    = phik + n*6.2831f;
      float a     = (3.0*(zi*pi+zj*pj+zk*pk)-(pi+pj+pk)*(zi+zj+zk))/(3.0*(zi*zi+zj*zj+zk*zk)-(zi+zj+zk)*(zi+zj+zk));
      float b     = ((pi+pj+pk) - a*(zi+zj+zk))/3.0;
      helix.dpdz_ = a;
      helix.phi0_ = b;
    }
  }


  /*
  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::initDpdzScan(const HitFlat& hitsf, CandHelix& helix, const std::vector<size_t>& hitsScan)
  {
    unsigned nshBest(0);
    for (const auto& is : hitsScan){
      for (const auto js : hitsScan){
        if (is>=js) continue;

        float dp = fabs(hitsf.phiC[js] - hitsf.phiC[is]);
        float dz = fabs(hitsf.z[js] - hitsf.z[is]);

        int ilmax = int( dz*maxDpdz_/6.2831f - dp/6.2831f);
        int ilmin = int(-dz*maxDpdz_/6.2831f - dp/6.2831f);
        for (int il=ilmin;il<=ilmax;++il){
           float dpdz = (dp + il*6.2831f)/dz;
           if (fabs(dpdz) > maxDpdz_ || fabs(dpdz) < minDpdz_) continue;

           unsigned nsh  = 0;
           float phi0    = hitsf.phiC[is] - hitsf.z[is]*dpdz;
           float pOffset = dpdz > 0 ? -3.1415f : 3.1415f;
           for (const auto& idx : helix.hits_){
              float dPhi = fmodf(fabsf(hitsf.z[idx]*dpdz +phi0 - hitsf.phiC[idx]),6.2831f);
              if (dPhi>3.1415f) dPhi = 6.2831f - dPhi;
              if (dPhi < maxDPhi_) nsh += hitsf.nsh[idx];
           }

           if (nsh <= nshBest ) continue;
           nshBest     = nsh;
           helix.dpdz_ = dpdz;
           helix.phi0_ = phi0;

           if (nshBest == helix.nStrawHits_) return;
        }
      }
    }
  }
  */

  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::cleanDpdzHits(const HitFlat& hitsf, CandHelix& helix)
  {
    auto pred = [&](size_t idx)
    {
       float dPhi = fmodf(fabsf(hitsf.z[idx]*helix.dpdz_ + helix.phi0_ - hitsf.phiC[idx]),6.2831f);
       return dPhi > maxDphi_ && 6.2831f-dPhi > maxDphi_;
    };

    helix.hits_.erase(std::remove_if(helix.hits_.begin(), helix.hits_.end(), pred),helix.hits_.end());
    helix.nStrawHits_=0;
    for (const auto& i : helix.hits_) helix.nStrawHits_ +=  hitsf.nsh[i];
  }


 //---------------------------------------------------------------------------------------------------------------------------
 void RobustMultiHelixFinder::fillPhiZFitter(const HitFlat& hitsf, const CandHelix& helix, LSFitter& zphiFitter)
 {
    zphiFitter.clear();
    for (const auto& ih : helix.hits_)
    {
      float base = hitsf.z[ih]*helix.dpdz_ + helix.phi0_;
      int   n    = int((base - helix.pOffset())/6.2831f);
      float dPhi = base - n*6.2831f - hitsf.phiC[ih];
      if      (dPhi > 3.1415f)  {dPhi -= 6.2831f; ++n;}
      else if (dPhi < -3.1415f) {dPhi += 6.2831f; --n;}
      if (fabs(dPhi) < maxDphi_) zphiFitter.add(hitsf.z[ih], hitsf.phiC[ih] + n*6.2831f, hitsf.nsh[ih]);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  bool RobustMultiHelixFinder::checkHelixQual(const HitFlat& hitsf, const CandHelix& helix, int nshCircle)
  {
     //-- Check  the ratio of surviving hits / circle hits, z entend, and phiC distribution
     float selHitRatio = float(helix.nStrawHits_)/float(nshCircle);
     if (selHitRatio < selHitRatio_) return false;

     //-- Check z length and gaps
     float zLength = hitsf.z[helix.hits_.back()] - hitsf.z[helix.hits_.front()];
     if (zLength < minDZTrk_) return false;

     float zgap(0);
     for (size_t i=1;i<helix.hits_.size();++i) zgap = std::max(zgap,hitsf.z[helix.hits_[i]] - hitsf.z[helix.hits_[i-1]]);
     if (zgap > maxDZGap_) return false;

     //-- Check gaps in phi
     std::vector<float> phiCvec;
     for (const auto k : helix.hits_) phiCvec.push_back(hitsf.phiC[k]);
     std::sort(phiCvec.begin(), phiCvec.end());

     std::vector<float> diffPhiC;
     for (size_t ii = 0; ii < phiCvec.size(); ++ii)
     {
       float dang  = std::abs(phiCvec[(ii+1) % phiCvec.size()] - phiCvec[ii]);
       diffPhiC.emplace_back(std::min(dang, 6.2831f - dang));
     }
     if (diffPhiC.size() > 2) {
       std::nth_element(diffPhiC.begin(), diffPhiC.begin()+1, diffPhiC.end(),std::greater<float>());
       if (diffPhiC[1] > maxDphiC_) return false;
     }

     return true;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::helixFilterChi2(const HitFlat& hitsf, CandHelix& helix)
  {
    auto chisqvec = chi2Helix(hitsf,helix);
    auto pred   = [&](size_t k) {return chisqvec[k] < maxChi2Filter_;};
    auto newEnd = std::stable_partition(helix.hits_.begin(), helix.hits_.end(), pred);
    helix.hits_.erase(newEnd, helix.hits_.end());

    helix.nStrawHits_= 0;
    for (const auto ih : helix.hits_) helix.nStrawHits_ += hitsf.nsh[ih];
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::helixRecoverChi2(const HitFlat& hitsf, CandHelix& helix)
  {
    helix.hits_.clear();
    float previousValue(1e6);

    auto chisqvec = chi2Helix(hitsf,helix);
    for (size_t k=0;k<hitsf.size();++k)
    {
       if (chisqvec[k] > maxChi2Recover_) continue;

       if (helix.hits_.empty() || hitsf.uid[k] != hitsf.uid[helix.hits_.back()]) {
          helix.hits_.emplace_back(k);
          previousValue = chisqvec[k];
          continue;
       }
       if (chisqvec[k] < previousValue){
          previousValue = chisqvec[k];
          helix.hits_.back() = k;
       }
    }

    helix.nStrawHits_= 0;
    helix.chisq_     = 0;
    for (const auto ih : helix.hits_) {
      helix.nStrawHits_ += hitsf.nsh[ih];
      helix.chisq_      += chisqvec[ih]*hitsf.nsh[ih];
    }
    helix.chisq_ /= helix.nStrawHits_;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::fitCircle(const HitFlat& hitsf, CandHelix& helix)
  {
    float w_(0),x_(0),x2_(0),x3_(0),y_(0),y2_(0),y3_(0),xy_(0),x2y_(0),xy2_(0);

    for (const auto ih : helix.hits_)
    {
      float w = hitsf.nsh[ih];
      float x = hitsf.x[ih];
      float y = hitsf.y[ih];
      x_   += x*w;
      x2_  += x*x*w;
      x3_  += x*x*x*w;
      y_   += y*w;
      y2_  += y*y*w;
      y3_  += y*y*y*w;
      xy_  += x*y*w;
      x2y_ += x*x*y*w;
      xy2_ += x*y*y*w;
      w_   += w;
    }

    float d11 = w_*xy_  - x_*y_;
    float d20 = w_*x2_  - x_*x_;
    float d30 = w_*x3_  - x2_*x_;
    float d21 = w_*x2y_ - x2_*y_;
    float d02 = w_*y2_  - y_*y_;
    float d03 = w_*y3_  - y2_*y_;
    float d12 = w_*xy2_ - x_*y2_;
    float den = 2*(d20*d02-d11*d11);
    if (fabs(den) < 1e-9) return;
    float cx  = ((d30+d12)*d02-(d03+d21)*d11)/den;
    float cy  = ((d03+d21)*d20-(d30+d12)*d11)/den;
    float rad = (x2_+y2_-2*cx*x_-2*cy*y_)/w_+cx*cx+cy*cy;

    helix.x_ = cx;
    helix.y_ = cy;
    helix.r_ = sqrtf(rad);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  float RobustMultiHelixFinder::helixCalcChi2(const HitFlat& hitsf, CandHelix& helix)
  {
    auto chisqvec = chi2Helix(hitsf,helix);

    float sum(0);
    for (auto k : helix.hits_) {sum += chisqvec[k]*hitsf.nsh[k]; }
    return sum/helix.nStrawHits_;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // One-step gradient descent to find best helix parameters
  void RobustMultiHelixFinder::refineHelix(HitFlat& hitsf, CandHelix& helix)
  {
    std::vector<float> val;
    std::vector<float> dispX{0.0f, -0.1f, 0.1f,  0.0f, 0.0f, 0.1f};
    std::vector<float> dispY{0.0f,  0.0f, 0.0f, -0.1f, 0.1f, 0.1f};

    CandHelix thelix(helix);
    for (size_t i=0;i<dispX.size();++i){
       float centerX = helix.x_ + dispX[i];
       float centerY = helix.y_ + dispY[i];

       float radius(0);
       for (const auto k : helix.hits_){
         float dx = (hitsf.x[k]-centerX);
         float dy = (hitsf.y[k]-centerY);
         radius += (dx*dx+dy*dy);
       }
       radius = sqrtf(radius/helix.hits_.size());

       thelix.x_ = centerX;
       thelix.y_ = centerY;
       thelix.r_ = radius;

       // There is no need to recalculate phiC as the values are essentially constant,
       // the Hessian is well approximated by just changing the radius and position of the helices
       val.push_back(helixCalcChi2(hitsf,thelix));
    }

    //Calculate Hessian entries
    float d2fx  = (val[2]+val[1]-2*val[0])/0.01;
    float d2fy  = (val[4]+val[3]-2*val[0])/0.01;
    float d2fxy = (val[5]-val[2]-val[4]+val[0])/0.01; //fast approximation
    float dfx   = (val[2]-val[0])/0.1;
    float dfy   = (val[4]-val[0])/0.1;
    float det   = d2fx*d2fy-d2fxy*d2fxy;
    if (fabs(det) <1e-9) return;

    float h11 = d2fy/det;
    float h22 = d2fx/det;
    float h12 = -d2fxy/det;
    float dx = h11*dfx+h12*dfy;
    float dy = h12*dfx+h22*dfy;
    if (fabs(dx) > helix.r_ || fabs(dy) > helix.r_) return;

    helix.x_ -= dx;
    helix.y_ -= dy;

    float radius(0);
    for (const auto idx : helix.hits_){
      float dx = (hitsf.x[idx]-helix.x_);
      float dy = (hitsf.y[idx]-helix.y_);
      radius  += (dx*dx+dy*dy);
    }
    helix.r_ = sqrtf(radius/helix.hits_.size());

    for (size_t n=0; n<hitsf.size(); ++n) hitsf.phiC[n] = polyAtan2(hitsf.y[n] - helix.y_, hitsf.x[n] - helix.x_);
    LSFitter zphiFitter;
    fillPhiZFitter(hitsf,helix,zphiFitter);
    helix.dpdz_ = zphiFitter.fa();
    helix.phi0_ = zphiFitter.fb();
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Estimate time parameters dzdt and t0 with error
  void RobustMultiHelixFinder::calcTimeProp(const HitFlat& hitsf, CandHelix& helix)
  {
    LSFitter ztFitter;
    for (const auto& k : helix.hits_) ztFitter.add(hitsf.z[k], hitsf.t[k],hitsf.nsh[k]);
    helix.dzdt_ = ztFitter.fa();
    helix.t0_   = ztFitter.fb();

    float s2(0),szb2(0),sn(0),zbar(ztFitter.xbar());
    for (const auto& k : helix.hits_)
    {
      float dz = hitsf.z[k] - zbar;
      float s  = helix.t0_ + hitsf.z[k]*helix.dzdt_ - hitsf.t[k];
      s2      += s*s*hitsf.nsh[k];
      szb2    += dz*dz*hitsf.nsh[k];
      sn      += hitsf.nsh[k];
    }
    if (sn<2) return;
    s2 /= (sn-2);

    helix.t0Err_= sqrtf((zbar*zbar/szb2+1.0/sn)*s2);
    helix.dzdtErr_ = sqrtf(s2/szb2);
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Print me!
  void RobustMultiHelixFinder::printHelix(const HelixSeed& helix)
  {
    std::cout<<"Helicity   "<<Helicity::name(helix.helix().helicity())<<std::endl;
    std::cout<<"Radius     "<<helix.helix().radius()<<std::endl;
    std::cout<<"Rcent      "<<helix.helix().rcent()<<std::endl;
    std::cout<<"Fcent      "<<helix.helix().fcent()<<std::endl;
    std::cout<<"Lambda     "<<helix.helix().lambda()<<std::endl;
    std::cout<<"fz0        "<<helix.helix().fz0()<<std::endl;
    std::cout<<"Time       "<<helix.t0().t0()<<"  +- "<<helix.t0().t0Err()<<std::endl;
    std::cout<<"Hits       "<<helix.hits().size()<<std::endl;
    std::cout<<"SH indx(0) ";
    for (const auto& ch : helix.hits()) std::cout<<ch.index(0)<<" ";
    std::cout<<std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Diagnosis
  void RobustMultiHelixFinder::fillDiag(mapHelix& helcols, const ComboHitCollection& chcol)
  {
    data_.iev_ = iev_;

    data_.Nch_ = chcol.size();
    for (unsigned ich=0; ich<chcol.size();++ich)
    {
       data_.chTime_[ich] = chcol[ich].correctedTime();
       data_.chX_[ich]    = chcol[ich].pos().x();
       data_.chY_[ich]    = chcol[ich].pos().y();
       data_.chZ_[ich]    = chcol[ich].pos().z();
       data_.chRad_[ich]  = sqrtf(data_.chX_[ich]*data_.chX_[ich]+data_.chY_[ich]*data_.chY_[ich]);
       data_.chPhi_[ich]  = chcol[ich].pos().phi();
       data_.chNhit_[ich] = chcol[ich].nStrawHits();
       data_.chUId_[ich]  = chcol[ich].strawId().uniquePanel();
       data_.chTerr_[ich] = chcol[ich].posRes(StrawHitPosition::trans);
       data_.chWerr_[ich] = chcol[ich].posRes(StrawHitPosition::wire);
       data_.chWDX_[ich]  = chcol[ich].uDir2D().x();
       data_.chWDY_[ich]  = chcol[ich].uDir2D().y();

       std::vector<StrawHitIndex> strawHitIdxs;
       chcol.fillStrawHitIndices(ich,strawHitIdxs);
       data_.chStrawId_.push_back(strawHitIdxs);
    }

    for (const auto& helix : *helcols[Helicity::poshel]){
      std::vector<int> hhits;
      for (auto hi : helix.hits()){
        for (size_t j=0;j<chcol.size();++j){
          if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
        }
      }
      data_.helhel_[data_.Nhel_]  = 1;
      data_.helrad_[data_.Nhel_]  = helix.helix().radius();
      data_.helrcen_[data_.Nhel_] = helix.helix().rcent();
      data_.helfcen_[data_.Nhel_] = helix.helix().fcent();
      data_.hellam_[data_.Nhel_]  = helix.helix().lambda();
      data_.helfz0_[data_.Nhel_]  = helix.helix().fz0();
      data_.heldzdt_[data_.Nhel_] = helix.helix().chi2dXY();
      data_.helnhi_[data_.Nhel_]  = hhits.size();
      data_.helchi2_[data_.Nhel_] = helix.helix().chi2dZPhi();
      data_.helhits_.push_back(hhits);
      ++data_.Nhel_;
   }

   for (const auto& helix : *helcols[Helicity::neghel]){
      std::vector<int> hhits;
      for (auto hi : helix.hits()){
        for (size_t j=0;j<chcol.size();++j){
          if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
        }
      }
      data_.helhel_[data_.Nhel_]  = -1;
      data_.helrad_[data_.Nhel_]  = helix.helix().radius();
      data_.helrcen_[data_.Nhel_] = helix.helix().rcent();
      data_.helfcen_[data_.Nhel_] = helix.helix().fcent();
      data_.hellam_[data_.Nhel_]  = helix.helix().lambda();
      data_.helfz0_[data_.Nhel_]  = helix.helix().fz0();
      data_.heldzdt_[data_.Nhel_] = helix.helix().chi2dXY();
      data_.helnhi_[data_.Nhel_]  = hhits.size();
      data_.helchi2_[data_.Nhel_] = helix.helix().chi2dZPhi();
      data_.helhits_.push_back(hhits);
      ++data_.Nhel_;
    }
  }


}

DEFINE_ART_MODULE(mu2e::RobustMultiHelixFinder);
