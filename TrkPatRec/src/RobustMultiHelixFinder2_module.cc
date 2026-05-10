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
#include "Offline/TrkPatRec/inc/RobustMultiHelixFinder_types.hh"

#include <numeric>
#include <stdint.h>
#include "/usr/include/valgrind/callgrind.h"


namespace {

  struct CaloCluHit {
    CaloCluHit() : x_(0),y_(0),z_(0),t_(0),e_(0),xyerr2_(0) {}
    CaloCluHit(float x, float y, float z, float t, float e, float xyerr) :
      x_(x),y_(y),z_(z),t_(t),e_(e),xyerr2_(xyerr*xyerr) {}

    float x_,y_,z_,t_,e_,xyerr2_;
  };
  using CaloCluHits = std::vector<CaloCluHit>;

  struct CandHelix {
    using artCptr = art::Ptr<mu2e::CaloCluster>;
    using hitsIdx = std::vector<size_t>;

    CandHelix() :
      nStrawHits_(0),nCaloCluHits_(0),x_(0),y_(0),r_(0),dpdz_(0),phi0_(0),dzdt_(0),t0_(0),caloIdx_(-1)
    {hits_.reserve(32);} ;

//FIXME THE 10 HERE NEEDS TO BE A VARIABLE AND SET THAT WAY
    const unsigned nhits()   const {return caloIdx_>-1? nStrawHits_ + 10 : nStrawHits_;}
    const float    pOffset() const {return dpdz_ > 0 ? -3.14159 : 3.14159;}

    void           clear()   {hits_.clear(); nStrawHits_ = nCaloCluHits_ = 0;x_=y_=r_=dpdz_=phi0_=dzdt_=t0_=0.0;}

    unsigned  nStrawHits_,nCaloCluHits_;
    float     x_,y_,r_,dpdz_,phi0_,dzdt_,t0_;
    hitsIdx   hits_;
    int       caloIdx_;
  };

  class LSFitter {
    public:
      LSFitter() : sn_(0),sx_(0),sx2_(0),sy_(0),sxy_(0){};

      void  add(float x, float y, float w=1.0) {sx_+=x*w; sx2_+=x*x*w; sy_+=y*w; sxy_+=x*y*w; sn_+=w;}
      void  clear()                            {sn_=sx_=sx2_=sy_=sxy_=0.0;}

      float fa()      {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sn_*sxy_-sx_*sy_) /(sn_*sx2_-sx_*sx_) : 0.0;}
      float fb()      {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sy_*sx2_-sx_*sxy_)/(sn_*sx2_-sx_*sx_) : 0.0;}
      bool  isValid() {return fabs(sn_)>1e-6;}

    private:
      float sn_,sx_,sx2_,sy_,sxy_;
  };
}



namespace mu2e {

  class RobustMultiHelixFinder2 : public art::EDProducer
  {
    public:
      using mapHelix        = std::map<Helicity,std::unique_ptr<HelixSeedCollection>>;
      using strawHitIndices = std::vector<StrawHitIndex>;
      using Config_types    = RobustMultiHelixFinderTypes::Config ;
      using Data_types      = RobustMultiHelixFinderTypes::Data_t;

      enum circleFitter {NoChoice, HyperFit, ChiSquared};

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> comboHitCollection    {Name("ComboHitCollection"),     Comment("ComboHit collection Name")    };
        fhicl::Atom<art::InputTag> timeClusterCollection {Name("TimeClusterCollection"),  Comment("TimeCluster collection Name") };
        fhicl::Atom<art::InputTag> caloClusterCollection {Name("CaloClusterCollection"),  Comment("CaloCluster collection Name") };
        fhicl::Sequence<int>       helicities            {Name("Helicities"),             Comment("Helicity values") };
        fhicl::Atom<float>         maxEdepHit            {Name("MaxEdepHit"),             Comment("Maximum eDep for hit - circle stage") };
        fhicl::Atom<float>         minRadCircle          {Name("MinRadCircle"),           Comment("Minimum circle radius") };
        fhicl::Atom<float>         maxRadCircle          {Name("MaxRadCircle"),           Comment("Maximum circle radius") };
        fhicl::Atom<float>         minDZCircle           {Name("MinDZCircle"),            Comment("Minimum Z distance between hits for circle fit") };
        fhicl::Atom<float>         minDXY2Circle         {Name("MinDXY2Circle"),          Comment("Minimum XY distance between hits for circle fit") };
        fhicl::Atom<float>         maxDXY2Circle         {Name("MaxDXY2Circle"),          Comment("Maximum XY distance between hits for circle fit") };
        fhicl::Atom<float>         maxDRCircle           {Name("MaxDRCircle"),            Comment("Maximum radial distance between circle and hit - circle stage") };
        fhicl::Atom<float>         maxChi2Circle         {Name("MaxChi2Circle"),          Comment("Maximum chi2 between circle and hit - circle stage") };
        fhicl::Atom<float>         resCRad               {Name("ResCRad"),                Comment("Radial resolution of track circle") };
        fhicl::Atom<float>         resCPerp              {Name("ResCPerp"),               Comment("Perpendicular resolution of track circle") };
        fhicl::Atom<float>         maxDphi               {Name("MaxDphi"),                Comment("Maximum phi difference between hit and dz/dphi line init fit") };
        fhicl::Atom<float>         maxDphiC              {Name("MaxDphiC"),               Comment("Maximum phiC difference between consecutive hits") };
        fhicl::Atom<float>         minDZTrk              {Name("MinDZTrk"),               Comment("Minimum DZ track length") };
        fhicl::Atom<float>         selHitRatio           {Name("SelHitRatio"),            Comment("Ratio of hits selected after dpdz fit") };
        fhicl::Atom<float>         minDpdz               {Name("MinDpdz"),                Comment("Minimum dp/dz helix candidate") };
        fhicl::Atom<float>         maxDpdz               {Name("MaxDpdz"),                Comment("Maximum dp/dz helix candidate") };
        fhicl::Atom<bool>          useDist2              {Name("UseDist2"),               Comment("Use distance based hit filtering and recovery") };
        fhicl::Atom<float>         maxDist2Filter        {Name("MaxDist2Filter"),         Comment("Maximum distance to filter hits") };
        fhicl::Atom<float>         maxDist2Recover       {Name("MaxDist2Recover"),        Comment("Maximum distance to recover hits") };
        fhicl::Atom<float>         maxChi2Filter         {Name("MaxChi2Filter"),          Comment("Maximum chi2 to filter hits") };
        fhicl::Atom<float>         maxChi2Recover        {Name("MaxChi2Recover"),         Comment("Maximum chi2 to recover hits") };
        fhicl::Atom<unsigned>      minStrawHits          {Name("MinStrawHits"),           Comment("Minimum number of Straw hits for a helix candidate") };
        fhicl::Atom<unsigned>      minnTotHits           {Name("MinnTotHits"),            Comment("Minimum number of Straw+Calo hits for a helix candidate") };
        fhicl::Atom<float>         ccMinEnergy           {Name("CaloClusterMinE"),        Comment("Minimum calo cluster energy") };
        fhicl::Atom<int>           ccWeight              {Name("CaloClusterWeight"),      Comment("Calo cluster weight ") };
        fhicl::Atom<float>         maxEDepAvg            {Name("MaxEDepAvg"),             Comment("Maximum EDep average")};
        fhicl::Atom<std::string>   fitCircleStr          {Name("FitCircleStrategy"),      Comment("Fit Circle algorhithm HyperFit or ChiSquared") };
        fhicl::Atom<unsigned>      nMaxTrkIter           {Name("NMaxTrkIter"),            Comment("Number of track finding iterations ") };
        fhicl::Atom<int>           diagLevel             {Name("DiagLevel"),              Comment("Diag level"), 0 };
        fhicl::Table<Config_types> diagPlugin            {Name("DiagPlugin"),             Comment("Diag Plugin config")};
      };
      explicit RobustMultiHelixFinder2(const art::EDProducer::Table<Config>& config);
      virtual void produce(art::Event& event);
      virtual void beginJob();


    private:
      const art::ProductToken<ComboHitCollection>    chToken_;
      const art::ProductToken<TimeClusterCollection> tcToken_;
      const art::ProductToken<CaloClusterCollection> ccToken_;
      float                                          maxEdepHit_;
      float                                          minRadCircle_;
      float                                          maxRadCircle_;
      float                                          minDZCircle_;
      float                                          minDXY2Circle_;
      float                                          maxDXY2Circle_;
      float                                          maxDRCircle_;
      float                                          maxChi2Circle_;
      float                                          resCRad_;
      float                                          resCPerp_;
      float                                          minDpdz_;
      float                                          maxDpdz_;
      float                                          maxDphi_;
      float                                          maxDphiC_;
      float                                          minDZTrk_;
      float                                          selHitRatio_;
      bool                                           useDist2_;
      float                                          maxDist2Filter_;
      float                                          maxDist2Recover_;
      float                                          maxChi2Filter_;
      float                                          maxChi2Recover_;
      unsigned                                       minStrawHits_;
      unsigned                                       minnTotHits_;
      float                                          ccMinEnergy_;
      int                                            ccWeight_;
      float                                          maxEDepAvg_;
      circleFitter                                   fitCircleStrategy_;
      unsigned                                       nMaxTrkIter_;
      const Calorimeter*                             cal_;
      int                                            iev_;
      int                                            diag_;
      std::unique_ptr<ModuleHistToolBase>            diagTool_;
      Data_types                                     data_;
      std::vector<Helicity>                          hels_;



      void      findAllHelices        (art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
                                       const art::ValidHandle<ComboHitCollection>& chH, const art::ValidHandle<CaloClusterCollection>& ccH);
      void      findHelicesInTC       (mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc,
                                       const ComboHitCollection& chcol, const CaloClusterCollection& cccol);
      CandHelix findHelixCandidate    (const ComboHitCollection& chcol, const strawHitIndices& hits, const CaloCluHits& caloHits);
      void      fitCircleAlg          (const CandHelix& circle, const ComboHitCollection& chcol, float& centerX, float& centerY, float& radius);
      void      fitCircleChi2         (const CandHelix& circle, const ComboHitCollection& chcol, float& centerX, float& centerY, float& radius);
      float     weight                (const ComboHit& ch,      const CandHelix& circle);
      void      checkCaloCluHitChi2   (const CaloCluHits& caloHits, CandHelix& helix);

      float     chi2XYCircleFast      (const ComboHit& ch, float centerX, float centerY, float radius);
      float     chi2XYCircle          (const ComboHit& ch, float centerX, float centerY, float radius);
      float     chi2XYHelix           (const ComboHit& ch, const CandHelix& helix);
      float     dist2Helix            (const ComboHit& ch, const CandHelix& helix);

      float     mindPhiTriplet        (const ComboHitCollection& chcol, int i, int j, int k, float centerX, float centerY);
      void      initDpdzTriplet       (const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitIdx);
      void      initDpdzScan          (const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitIdx);
      void      fillPhiZFitter        (const ComboHitCollection& chcol, const std::vector<float>& phiCcol, const CandHelix& helix, LSFitter& zphiFitter);
      void      dpdzFit               (const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix);
      void      cleanDpdzHits         (const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix);
      float     timeUncertainty       (const ComboHitCollection& chcol, const CandHelix& helix);

      void      helixRecoverDist      (const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit);
      void      helixRecoverChi2      (const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit);
      void      helixRecover          (const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit);

      void      helixFilterDist       (const ComboHitCollection& chcol, CandHelix& helix);
      void      helixFilterChi2       (const ComboHitCollection& chcol, CandHelix& helix);
      void      helixFilter           (const ComboHitCollection& chcol, CandHelix& helix);

      float     helixSumDist2         (const ComboHitCollection& chcol, CandHelix& helix);
      float     helixSumChi2          (const ComboHitCollection& chcol, CandHelix& helix);
      float     helixSumResid         (const ComboHitCollection& chcol, CandHelix& helix);

      void      refineXY              (const ComboHitCollection& chcol, std::vector<float>& phiCcol, CandHelix& helix);

      void      filterDuplicateHelices(HelixSeedCollection& helices);

      void      printHelix            (const HelixSeed& helix);
      void      fillDiag              (mapHelix& helcols, const ComboHitCollection& chcol);


      void printVec(const std::vector<float>& vec) {for (const auto& v : vec) {std::cout<<v<<" ";} std::cout<<"\n";}
      void printVec(const std::vector<int>& vec) {for (const auto& v : vec) {std::cout<<v<<" ";}std::cout<<"\n";}
      void printVec(const std::vector<unsigned>& vec) {for (const auto& v : vec) {std::cout<<v<<" ";} std::cout<<"\n";}
      void printVec(const strawHitIndices& vec) {for (const auto& v : vec) {std::cout<<v<<" ";} std::cout<<"\n";}

  };


  RobustMultiHelixFinder2::RobustMultiHelixFinder2(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    chToken_       {consumes<ComboHitCollection>   (config().comboHitCollection())},
    tcToken_       {consumes<TimeClusterCollection>(config().timeClusterCollection())},
    ccToken_       {consumes<CaloClusterCollection>(config().caloClusterCollection())},
    maxEdepHit_       (config().maxEdepHit()),
    minRadCircle_     (config().minRadCircle()),
    maxRadCircle_     (config().maxRadCircle()),
    minDZCircle_      (config().minDZCircle()),
    minDXY2Circle_    (config().minDXY2Circle()),
    maxDXY2Circle_    (config().maxDXY2Circle()),
    maxDRCircle_      (config().maxDRCircle()),
    maxChi2Circle_    (config().maxChi2Circle()),
    resCRad_          (config().resCRad()),
    resCPerp_         (config().resCPerp()),
    minDpdz_          (config().minDpdz()),
    maxDpdz_          (config().maxDpdz()),
    maxDphi_          (config().maxDphi()),
    maxDphiC_         (config().maxDphiC()),
    minDZTrk_         (config().minDZTrk()),
    selHitRatio_      (config().selHitRatio()),
    useDist2_         (config().useDist2()),
    maxDist2Filter_   (config().maxDist2Filter()),
    maxDist2Recover_  (config().maxDist2Recover()),
    maxChi2Filter_    (config().maxChi2Filter()),
    maxChi2Recover_   (config().maxChi2Recover()),
    minStrawHits_     (std::max(config().minStrawHits(),3u)),
    minnTotHits_      (std::max(config().minnTotHits(),3u)),
    ccMinEnergy_      (config().ccMinEnergy()),
    ccWeight_         (config().ccWeight()),
    maxEDepAvg_       (config().maxEDepAvg()),
    fitCircleStrategy_(circleFitter::NoChoice),
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

    if      (config().fitCircleStr()=="HyperFit")   fitCircleStrategy_ = circleFitter::HyperFit;
    else if (config().fitCircleStr()=="ChiSquared") fitCircleStrategy_ = circleFitter::ChiSquared;
    else    throw cet::exception("CATEGORY")<< "RobustMultiHelixFinder2: unrecognixed FitCirclestrategy specified";
  }


  //--------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::beginJob(){
    if (diag_){
       art::ServiceHandle<art::TFileService> tfs;
       diagTool_->bookHistograms(tfs);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::produce(art::Event& event )
  {
    if (diag_>0) std::cout<<"Event "<<event.id().event()<<std::endl;
    iev_ = event.id().event();

    const auto& tcH = event.getValidHandle(tcToken_);
    const auto& chH = event.getValidHandle(chToken_);
    const auto& ccH = event.getValidHandle(ccToken_);

CALLGRIND_START_INSTRUMENTATION;
    mapHelix helcols;
    for (const auto& hel : hels_) helcols[hel] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());

    findAllHelices(event, helcols,tcH, chH, ccH);
CALLGRIND_STOP_INSTRUMENTATION;

    for (const auto& hel : hels_) event.put(std::move(helcols[hel]),Helicity::name(hel));
  }




  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in an event (hopefully), some timeCluster implementations have overlapping content, so filter duplicates
  void RobustMultiHelixFinder2::findAllHelices(art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
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
  void RobustMultiHelixFinder2::filterDuplicateHelices(HelixSeedCollection& helices)
  {
    for (auto first = helices.begin(); first != helices.end(); ++first){
      for (auto second = std::next(first); second != helices.end(); ++second){

        size_t result(0);
        auto helix1 = first->hits();
        auto helix2 = second->hits();
        auto iter1 = helix1.begin();
        auto iter2 = helix2.begin();
        while (iter1!=helix1.end() && iter2!=helix2.end()) {
           if (iter1->index(0)<iter2->index(0)) ++iter1;
           else if (iter2->index(0)<iter1->index(0)) ++iter2;
           else {++result; ++iter1; ++iter2;}
        }

        if (result != helix1.size() && result != helix2.size()) continue;

        if (first->hits().size()>second->hits().size()) second->_hhits.clear();
        else first->_hhits.clear();
      }
    }

    auto predEmpty = [](const HelixSeed& hel) {return hel.hits().empty();};
    helices.erase(std::remove_if(helices.begin(),helices.end(),predEmpty),helices.end());
  }








  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in a timeCluster for both helicities. Can find up to nMaxTrkIter helices for each helicity
  void RobustMultiHelixFinder2::findHelicesInTC(mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc,
                                               const ComboHitCollection& chcol, const CaloClusterCollection& cccol)
  {

if (chcol.empty()){
CandHelix h;
std::vector<float> ph;
std::vector<StrawHitIndex> u;
float v(0),vv(0),vvv(0);
dpdzFit(chcol,ph,h);
fitCircleChi2(h,chcol, v,vv,vvv);
fitCircleAlg(h,chcol,v,vv,vvv);
initDpdzScan(chcol,ph,h,u);
initDpdzTriplet(chcol,ph,h,u);
}


    //for (auto i : tc.hits()) std::cout<<i<<" ";
    //std::cout<<"\n";
    //Calculate average time in tccol Hits and create CaloCluHit collection
    float taverage(0);
    for (auto& ich : tc.hits()) taverage += chcol[ich].correctedTime()/tc.hits().size();

    CaloCluHits caloHits;
    for (const auto& ccalo : cccol){
      if (fabs(ccalo.time()-taverage) > 40) continue;
      auto ccPosTrk = cal_->geomUtil().mu2eToTracker(cal_->geomUtil().diskFFToMu2e(ccalo.diskID(),ccalo.cog3Vector()));
      caloHits.emplace_back(CaloCluHit(ccPosTrk.x(),ccPosTrk.y(),ccPosTrk.z(),ccalo.time(),ccalo.energyDep(),100));
    }


    std::vector<StrawHitIndex> usedHits;
    for (size_t it=0;it<nMaxTrkIter_;++it){

      size_t strawHitsToGo(0);
      strawHitIndices hitsToProcess;
      hitsToProcess.reserve(64);
      for (const auto& ich : tc.hits()) {
        if (find(usedHits.begin(),usedHits.end(),ich) != usedHits.end()) continue;
        hitsToProcess.emplace_back(ich);
        strawHitsToGo += chcol[ich].nStrawHits();
      }
      if (strawHitsToGo < minStrawHits_) break;


      auto pred = [&chcol](const auto& i, const auto& j) {return chcol[i].strawId().uniquePanel()<chcol[j].strawId().uniquePanel();};
      sort(hitsToProcess.begin(),hitsToProcess.end(),pred);

      CandHelix bestHelix = findHelixCandidate(chcol, hitsToProcess, caloHits);
      if (bestHelix.nStrawHits_ < minStrawHits_) break;


      const Helicity bestHelicity = bestHelix.dpdz_>0 ? Helicity::poshel: Helicity::neghel;
      for (const auto& ich : bestHelix.hits_) usedHits.push_back(ich);

      //FIX THESE, THEY CAN BE SAVED IN THE HELIX SINCE THE CALCULATION IS MADE EARLIER
      float chi2dXY(0), chi2dZPhi(0);
      for (const auto& ich : bestHelix.hits_) {
         chi2dXY   += chi2XYCircle(chcol[ich],bestHelix.x_,bestHelix.y_,bestHelix.r_)*chcol[ich].nStrawHits();
         chi2dZPhi += chi2XYHelix(chcol[ich],bestHelix)*chcol[ich].nStrawHits();
      }
      chi2dXY   /= bestHelix.nStrawHits_;
      chi2dZPhi /= bestHelix.nStrawHits_;
//Dirty hack to save particle propagation direction, will be gone when we have updated the data products
chi2dXY = bestHelix.dzdt_;


      float Rcent  = sqrt(bestHelix.x_*bestHelix.x_+bestHelix.y_*bestHelix.y_);
      float Fcent  = polyAtan2(bestHelix.y_,bestHelix.x_);
      float t0err  = timeUncertainty(chcol,bestHelix);

      HelixSeed hseed;
      hseed._helix = RobustHelix(Rcent,Fcent,bestHelix.r_,bestHelix.dpdz_,bestHelix.phi0_);
      hseed._hhits.setParent(chcol.parent());
      hseed._helix._helicity  = bestHelicity;
      hseed._helix._chi2dXY   = chi2dXY;
      hseed._helix._chi2dZPhi = chi2dZPhi;
      hseed._t0 = TrkT0(bestHelix.t0_,t0err);
      for (const auto& ich : bestHelix.hits_) hseed._hhits.emplace_back(chcol[ich]);
      hseed._status.merge(TrkFitFlag::MPRHelix);
      hseed._status.merge(TrkFitFlag::helixOK);
      hseed._timeCluster = tcArtPtr;
      hseed._eDepAvg = hseed._hhits.eDepAvg();
      if (hseed._eDepAvg > maxEDepAvg_) continue;

      helcols[bestHelicity]->push_back(std::move(hseed));
    }



}


  //---------------------------------------------------------------------------------------------------------------------------
  // Loop over triplets and find the circle with the maximum hits using the distance between hit and circle as (fast) metric
  // or chi2. Exclude compton hits in triplet search, but include them in circle assignment. Add calo cluster if present
  CandHelix RobustMultiHelixFinder2::findHelixCandidate(const ComboHitCollection& chcol, const strawHitIndices& hits, const CaloCluHits& caloHits)
  {

    float chi2Best(1e6);
    CandHelix helix, bestHelix;
    std::vector<float> phiCcol(chcol.size(),-999);

    //Find the best circle
    size_t stride(1);
    for (size_t i=0; i+2<hits.size(); i+=stride){
      const auto& chi  = chcol[hits[i]];
      float rad2x1y1   = chi.pos().x()*chi.pos().x()+chi.pos().y()*chi.pos().y();
      float x1         = chi.pos().x();
      float y1         = chi.pos().y();
      if (chi.energyDep() > maxEdepHit_) continue;

      for (size_t j=i+1; j<hits.size(); j+=stride){
        const auto& chj    = chcol[hits[j]];
        float rad2x2y2     = chj.pos().x()*chj.pos().x()+chj.pos().y()*chj.pos().y();
        float x2           = chj.pos().x()-chi.pos().x();
        float y2           = chj.pos().y()-chi.pos().y();
        float x12x12y12y12 = x2*x2+y2*y2;
        if (chj.energyDep() > maxEdepHit_)     continue;
        if (chj.pos().z()-chi.pos().z() < minDZCircle_) continue;
        if (x12x12y12y12 < minDXY2Circle_ || x12x12y12y12 > maxDXY2Circle_) continue;

        for (size_t k=j+1; k<hits.size();k+=stride){
          const auto& chk    = chcol[hits[k]];
          float rad2x3y3     = chk.pos().x()*chk.pos().x()+chk.pos().y()*chk.pos().y();
          float x3           = chk.pos().x()-chi.pos().x();
          float y3           = chk.pos().y()-chi.pos().y();
          float x13x13y13y13 = x3*x3+y3*y3;
          if (chk.energyDep() > maxEdepHit_)     continue;
          if (chk.pos().z()-chj.pos().z() < minDZCircle_) continue;
          if (x13x13y13y13 < minDXY2Circle_ || x13x13y13y13 > maxDXY2Circle_) continue;

          float x4           = chk.pos().x()-x2;
          float y4           = chk.pos().y()-y2;
          float x23x23y23y23 = x4*x4+y4*y4;
          if (x23x23y23y23 < minDXY2Circle_ || x23x23y23y23 > maxDXY2Circle_) continue;

          float denominator = 2*(x2*y3 - x3*y2);
          float numeratorX  = x12x12y12y12*y3 - x13x13y13y13*y2;
          float numeratorY  = x13x13y13y13*x2 - x12x12y12y12*x3;

          float centerX     = numeratorX/denominator+x1;
          float centerY     = numeratorY/denominator+y1;
          float radius      = sqrtf((x1-centerX)*(x1-centerX) + (y1-centerY)*(y1-centerY));
          float radCenter2  = centerX*centerX+centerY*centerY;

          if (radius < minRadCircle_ || radius > maxRadCircle_ ) continue;
          if (radCenter2 > rad2x1y1  || radCenter2 > rad2x2y2 || radCenter2 > rad2x3y3) continue;

          //Check if the triplets roughly align here
          //float minDP = mindPhiTriplet(chcol, hits[i], hits[j], hits[k], centerX, centerY);
          //if (minDP > 0.8) continue;

          // select hits close enough to the circle
          helix.clear();
          helix.x_ = centerX;
          helix.y_ = centerY;
          helix.r_ = radius;

          float previousValue(1e6);
          for (const auto ih : hits){
            const ComboHit& ch = chcol[ih];
            float value        = (ih==hits[i] || ih==hits[j] || ih==hits[k]) ? 0 : chi2XYCircleFast(ch,centerX,centerY,radius);
            //float dx           = ch.pos().x() - centerX;
            //float dy           = ch.pos().y() - centerY;
            //float value        = (ih==hits[i] || ih==hits[j] || ih==hits[k]) ? 0 : fabs(sqrtf(dx*dx+dy*dy)-radius);

            //if (value > maxDRCircle_) continue;
            if (value > maxChi2Circle_) continue;

            // if the previous hit as the same uid, then we need to keep only the best one
            if (helix.hits_.empty() || ch.strawId().uniquePanel() != chcol[helix.hits_.back()].strawId().uniquePanel()){
               previousValue = value;
               helix.nStrawHits_ += chcol[ih].nStrawHits();
               helix.hits_.emplace_back(ih);
            }
            else if (value < previousValue){
               previousValue = value;
               helix.nStrawHits_ -= chcol[helix.hits_.back()].nStrawHits();
               helix.nStrawHits_ += chcol[ih].nStrawHits();
               helix.hits_.back() = ih;
            }
          }


          // add the calorimeter here, based on radius only
          double caloDR(50);
          for (size_t icalo=0;icalo<caloHits.size();++icalo) {
            const auto& ccalo = caloHits[icalo];
            double dr = fabs(sqrtf((ccalo.x_ - centerX)*(ccalo.x_ - centerX) + (ccalo.y_ - centerY)*(ccalo.y_ - centerY)) - radius);
            if (dr < caloDR)
            {
               caloDR = dr;
               helix.caloIdx_ = icalo;
            }
          }

          if (helix.nStrawHits_ < minStrawHits_) continue;
          //if (helix.nhits() < minnTotHits_ || helix.nhits() < bestHelix.nhits()) continue;
          auto nshCircle = helix.nStrawHits_;




          //Estimate dp/dz and remove incompatible hits
          for (auto ich : helix.hits_) phiCcol[ich] = polyAtan2(chcol[ich].pos().y() - helix.y_, chcol[ich].pos().x() - helix.x_);

          std::vector<StrawHitIndex> hitsScan{hits[i], hits[j], hits[k]};
          initDpdzScan(chcol, phiCcol, helix, hitsScan);
          //initDpdzTriplet(chcol, phiCcol, helix, hitsScan);
          if (fabs(helix.dpdz_) < minDpdz_ || fabs(helix.dpdz_) > maxDpdz_) continue;

          cleanDpdzHits(chcol,phiCcol,helix);
          checkCaloCluHitChi2(caloHits,helix);

          if (helix.nStrawHits_ < minStrawHits_) continue;
          if (helix.nhits() < minnTotHits_ || helix.nhits() < bestHelix.nhits()) continue;


          //-- Perform full dp/dz fit
          //fitCircle(chcol,helix);
          //for (auto ich : helix.hits_) phiCcol[ich] = polyAtan2(chcol[ich].y - helix.y_, chcol[ich].x - helix.x_);
          //dpdzFit(chcol,phiCcol, helix);


          //plotTrip(chcol,helix.hits_,helix.x_,helix.y_,helix.r_,helix.dpdz_,helix.phi0_);
          //cout<<"Fit dpdz "<<helix.nsh_<<endl;



          // checks on the ratio of surviving hits / circle hits, z entend, and phiC distribution
          //FIXME
          float selHitRatio = float(helix.nStrawHits_)/float(nshCircle);
          if (selHitRatio < selHitRatio_) continue;

          // checks the z length
          float zLength = chcol[helix.hits_.back()].pos().z() - chcol[helix.hits_.front()].pos().z();
          if (zLength < minDZTrk_) continue;

          float zgap(0);
          for (size_t i=1;i<helix.hits_.size();++i) zgap = std::max(zgap,chcol[helix.hits_[i]].pos().z() -chcol[helix.hits_[i-1]].pos().z());
          if (zgap > 1500) continue;


          //check if there are holes in phiC
          std::vector<float> phiCvec,diffPhiC;
          for (const auto& k : helix.hits_) phiCvec.emplace_back(phiCcol[k]);
          std::sort(phiCvec.begin(),phiCvec.end());

          for (size_t i=0;i<phiCvec.size();++i){
            float dang = std::abs(phiCvec[(i+1)%phiCvec.size()]-phiCvec[i]);
            float delta = std::min(dang, 6.2831f-dang);
            diffPhiC.push_back(delta);
          }

          std::sort(diffPhiC.rbegin(),diffPhiC.rend());
          if (diffPhiC.size()>2 && diffPhiC[1]>maxDphiC_) continue;

//CBE cout<<"Passed all cuts\n";
          //Check if the clorimeter time is roughly consistent with the hits
          //FIXME
          if (helix.caloIdx_>-1) {
            LSFitter ztFitter;
              for (const auto& k : helix.hits_) ztFitter.add(chcol[k].pos().z(), chcol[k].correctedTime(),chcol[k].nStrawHits());

            float deltaTime= ztFitter.fa()*caloHits[helix.caloIdx_].z_ + ztFitter.fb() - caloHits[helix.caloIdx_].t_;
            if (fabs(deltaTime)>30)  continue;
          }


          // Simplest solution, filter, refit and associate hits
                //-- Filter hits, refine helic parameters, and select final hit list
          helixFilter(chcol,helix);

          //--Refine helix
          refineXY(chcol,phiCcol,helix);
          //try this to go faster
          //fitCircle(chcol,helix);
          //for (auto ich : helix.hits_) phiCcol[ich] = polyAtan2(chcol[ich].pos().y() - helix.y_, chcol[ich].pos().x() - helix.x_);
          //dpdzFit(chcol,phiCcol,helix);

          if (fabs(helix.dpdz_) < minDpdz_ || fabs(helix.dpdz_) > maxDpdz_) continue;
          if (helix.r_ < minRadCircle_ || helix.r_ > maxRadCircle_ ) continue;


          //Add all compatible hits
          helixRecover(chcol, helix, hits);
          checkCaloCluHitChi2(caloHits, helix);

          float chindf = helixSumChi2(chcol,helix)/helix.hits_.size();
          if (helix.nhits() > bestHelix.nhits() || (helix.nhits() == bestHelix.nhits() && chindf < chi2Best)){
            chi2Best = chindf;
            bestHelix = helix;
          }
        }
      }
    }

    if (bestHelix.nhits() < minnTotHits_ ) {bestHelix.clear(); return bestHelix;}

    // recalculate dz/dt and t0 for the best helix
    // TODO: also calculate the error

    LSFitter ztFitter;
    for (const auto& k : bestHelix.hits_) ztFitter.add(chcol[k].pos().z(), chcol[k].correctedTime(),chcol[k].nStrawHits());
    bestHelix.dzdt_ = ztFitter.fa();
    bestHelix.t0_   = ztFitter.fb();

    return bestHelix;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::checkCaloCluHitChi2(const CaloCluHits& caloHits, CandHelix& helix)
  {
     helix.caloIdx_ = -1;

     //FIXME
     float chi2min(5);
     //float chi2min(maxChi2Calo);
     for (size_t icalo=0;icalo<caloHits.size();++icalo){
        const auto caloHit = caloHits[icalo];
        float phiAtZ = caloHit.z_*helix.dpdz_ + helix.phi0_;
        float dhx    = caloHit.x_ - (helix.x_ + helix.r_*cos(phiAtZ));
        float dhy    = caloHit.y_ - (helix.y_ + helix.r_*sin(phiAtZ));
        float chi2 = dhx*dhx/caloHit.xyerr2_ + dhy*dhy/caloHit.xyerr2_;
        if (chi2 < chi2min) {chi2min = chi2; helix.caloIdx_ = icalo;}
     }
  }







  //---------------------------------------------------------------------------------------------------------------------------
  // Algebraic weighted circle fit using weights from hit resolution (adapted from HyperFit from N. Chernov)
  void RobustMultiHelixFinder2::fitCircleAlg(const CandHelix& circle, const ComboHitCollection& chcol, float& centerX,
                                            float& centerY, float& radius)
  {
    std::vector<float> hitWeight(chcol.size(),0);
    for (const auto& ich : circle.hits_) hitWeight[ich] = weight(chcol[ich],circle);//*chcol[ich].nStrawHits();


    float sumWeights(0),xmean(0),ymean(0);
    for (const auto& ich : circle.hits_) {
      xmean      += chcol[ich].pos().x()*hitWeight[ich];
      ymean      += chcol[ich].pos().y()*hitWeight[ich];
      sumWeights += hitWeight[ich];
    }
    xmean /= sumWeights;
    ymean /= sumWeights;

    float Mxx(0),Mxy(0),Myy(0),Mxz(0),Myz(0),Mzz(0);
    for (const auto& ich : circle.hits_){
        float xx = chcol[ich].pos().x()-xmean;
        float yy = chcol[ich].pos().y()-ymean;
        float zz = xx*xx+yy*yy;
        float w  = hitWeight[ich];
        Mxy += xx*yy*w;
        Mxx += xx*xx*w;
        Myy += yy*yy*w;
        Mxz += xx*zz*w;
        Myz += yy*zz*w;
        Mzz += zz*zz*w;
    }
    Mxx /= sumWeights;
    Myy /= sumWeights;
    Mxy /= sumWeights;
    Mxz /= sumWeights;
    Myz /= sumWeights;
    Mzz /= sumWeights;

    float Mz     = Mxx + Myy;
    float Cov_xy = Mxx*Myy - Mxy*Mxy;
    float Var_z  = Mzz - Mz*Mz;
    float A2     = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    float A1     = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    float A0     = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    float A22    = A2+A2;

    unsigned iterMAX = 10;
    float x(0.0),y(A0);
    for (unsigned iter=0; iter<iterMAX; ++iter)  {
        float Dy = A1 + x*(A22 + 16.*x*x);
        float xnew = x - y/Dy;
        if (abs(xnew-x) < 1e-4 || abs(x) > 1e10) break;
        float ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
        if (abs(ynew)> abs(y)) break;
        x = xnew;  y = ynew;
    }

    float DET = x*x - x*Mz + Cov_xy;
    float centerXMean = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
    float centerYMean = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;

    centerX = centerXMean + xmean;
    centerY = centerYMean + ymean;
    radius  = sqrt(centerXMean*centerXMean + centerYMean*centerYMean + Mz - x - x);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Circle fit minimizing chi2 = (delta R2)^2
  void RobustMultiHelixFinder2::fitCircleChi2(const CandHelix& circle, const ComboHitCollection& chcol, float& centerX,
                                             float& centerY, float& radius)
  {
    float w_(0),x_(0),x2_(0),x3_(0),y_(0),y2_(0),y3_(0),xy_(0),x2y_(0),xy2_(0);

    for (const auto& ich : circle.hits_) {
      float w = weight(chcol[ich],circle);
      float x = chcol[ich].pos().x();
      float y = chcol[ich].pos().y();
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
    centerX  = ((d30+d12)*d02-(d03+d21)*d11)/den;
    centerY  = ((d03+d21)*d20-(d30+d12)*d11)/den;
    float c  = (x2_+y2_-2*centerX*x_-2*centerY*y_)/w_;
    radius   = sqrt(c+centerX*centerX+centerY*centerY);
  }







  //---------------------------------------------------------------------------------------------------------------------------
  // Simple transverse distance to helix calculation
  float RobustMultiHelixFinder2::dist2Helix(const ComboHit& ch, const CandHelix& helix)
  {
    float phiPred = ch.pos().z()*helix.dpdz_ + helix.phi0_;
    float xPred   = helix.x_ + helix.r_*cos(phiPred);
    float yPred   = helix.y_ + helix.r_*sin(phiPred);
    float dx      = ch.pos().x() - xPred;
    float dy      = ch.pos().y() - yPred;
    return dx*dx+dy*dy;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 based on circle fit
  float RobustMultiHelixFinder2::chi2XYCircleFast(const ComboHit& ch, float centerX, float centerY, float radius)
  {
    float dx     = ch.pos().x() - centerX;
    float dy     = ch.pos().y() - centerY;
    float r      = sqrt(dx*dx+dy*dy);
    float dr     = r - radius;
    float rwdot  = (ch.vDir().x()*dx + ch.vDir().y()*dy)/r;
    float rwdot2 = rwdot*rwdot;
    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float rres2  = werr*werr*rwdot2 + terr*terr*(1.0-rwdot2);
    float chisq  = dr*dr/rres2;
    return chisq;
  }

  // This version includes the uncertainty along the circle
  float RobustMultiHelixFinder2::chi2XYCircle(const ComboHit& ch, float centerX, float centerY, float radius)
  {
    float cradres_(resCRad_*radius),cperpres_(resCPerp_*radius);

    float dx     = ch.pos().x() - centerX;
    float dy     = ch.pos().y() - centerY;
    float dnorm  = sqrtf(dx*dx+dy*dy);
    float cdd    = (dx*ch.vDir().x()+dy*ch.vDir().y())/dnorm;
    float cpd    = (-dy*ch.vDir().x()+dx*ch.vDir().y())/dnorm;

    float phiAtZ = polyAtan2(ch.pos().y()-centerY,ch.pos().x()-centerX);
    float dhx    = ch.pos().x() - (centerX + radius*cos(phiAtZ));
    float dhy    = ch.pos().y() - (centerY + radius*sin(phiAtZ));
    float dtrans = fabs(-ch.vDir().y()*dhx + ch.vDir().x()*dhy);
    float dwire  = fabs(ch.vDir().x()*dhx + ch.vDir().y()*dhy);

    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float wres2 = werr*werr + cradres_*cradres_*cdd*cdd + cperpres_*cperpres_*cpd*cpd;
    float tres2 = terr*terr + cradres_*cradres_*cpd*cpd + cperpres_*cperpres_*cdd*cdd;
    float chisq = dwire*dwire/wres2 + dtrans*dtrans/tres2;

    return chisq;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 with full helix hypothesis
  float RobustMultiHelixFinder2::chi2XYHelix(const ComboHit& ch, const CandHelix& helix)
  {
    float dx      = ch.pos().x()-helix.x_;
    float dy      = ch.pos().y()-helix.y_;
    float dnorm   = sqrtf(dx*dx+dy*dy);
    float cdd     = (dx*ch.vDir().x()+dy*ch.vDir().y())/dnorm;
    float cpd     = (-dy*ch.vDir().x()+dx*ch.vDir().y())/dnorm;

    float phiAtZ  = ch.pos().z()*helix.dpdz_+helix.phi0_;
    float dhx     = ch.pos().x() - (helix.x_ + helix.r_*cos(phiAtZ));
    float dhy     = ch.pos().y() - (helix.y_ + helix.r_*sin(phiAtZ));
    float dtrans  = fabs(-ch.vDir().y()*dhx + ch.vDir().x()*dhy);
    float dwire   = fabs(ch.vDir().x()*dhx + ch.vDir().y()*dhy);

    float werr    = ch.posRes(StrawHitPosition::wire);
    float terr    = ch.posRes(StrawHitPosition::trans);
    float cradres = resCRad_*helix.r_;
    float cperpres= resCPerp_*helix.r_;
    float wres2   = werr*werr + cradres*cradres*cdd*cdd + cperpres*cperpres*cpd*cpd;
    float tres2   = terr*terr + cradres*cradres*cpd*cpd + cperpres*cperpres*cdd*cdd;
    return dwire*dwire/wres2 + dtrans*dtrans/tres2;
}




  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate weight based on circle fit and straw resolution
  float RobustMultiHelixFinder2::weight(const ComboHit& ch, const CandHelix& circle)
  {
    float dx     = ch.pos().x() - circle.x_;
    float dy     = ch.pos().y() - circle.y_;
    float rwdot  = (ch.uDir2D().x()*dx + ch.uDir2D().y()*dy);
    float costh2 = rwdot*rwdot/(dx*dx+dy*dy);
    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float res2   = werr*werr*costh2 + terr*terr*(1.0-costh2);
    return 1.0/res2;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  float RobustMultiHelixFinder2::mindPhiTriplet(const ComboHitCollection& chcol, int i, int j, int k, float centerX, float centerY)
  {
    float zi   = chcol[i].pos().z();
    float zj   = chcol[j].pos().z();
    float zk   = chcol[k].pos().z();
    float phii = polyAtan2(chcol[i].pos().y()-centerY,chcol[i].pos().x()-centerX);
    float phij = polyAtan2(chcol[j].pos().y()-centerY,chcol[j].pos().x()-centerX);
    float phik = polyAtan2(chcol[k].pos().y()-centerY,chcol[k].pos().x()-centerX);
    int ilmax  = int((zj-zi)*maxDpdz_/6.2831 - (phij-phii)/6.2831);
    int ilmin  = int((zi-zj)*maxDpdz_/6.2831 - (phij-phii)/6.2831);

    float mindPhi(1e6);
    for (int il=ilmin; il<=ilmax; ++il){
       float dpdz  = (phij-phii + il*6.2831)/(zj-zi);
       if (fabs(dpdz) > maxDpdz_ || fabs(dpdz)<minDpdz_) continue;

       float pOffset = dpdz > 0 ? -3.14159 : 3.14159;
       float phi0    = phii - zi*dpdz;
       int   n       = int((zk*dpdz + phi0 - pOffset)/6.2831);
       float dPhi    = fabs(zk*dpdz + phi0 - n*6.2831 - phik);
       mindPhi = std::min(dPhi,mindPhi);
    }
    return mindPhi;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::initDpdzTriplet(const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitIdx)
  {
     const float& zi   = chcol[hitIdx[0]].pos().z();
     const float& zj   = chcol[hitIdx[1]].pos().z();
     const float& zk   = chcol[hitIdx[2]].pos().z();
     const float& phii = phiCcol[hitIdx[0]];
     const float& phij = phiCcol[hitIdx[1]];
     const float& phik = phiCcol[hitIdx[2]];
     int ilmax  = int((zj-zi)*maxDpdz_/6.2831 - (phij-phii)/6.2831);
     int ilmin  = int((zi-zj)*maxDpdz_/6.2831 - (phij-phii)/6.2831);

     float mindPhi(maxDphi_);
     for (int il=ilmin; il<=ilmax; ++il){

        float dpdz = (phij-phii + il*6.2831)/(zj-zi);
        if (fabs(dpdz) > maxDpdz_ || fabs(dpdz)<minDpdz_) continue;

        float pOffset = dpdz > 0 ? -3.14159 : 3.14159;
        float phi0    = phii - zi*dpdz;
        int   n       = int((zk*dpdz + phi0 - pOffset)/6.2831);
        float dPhi    = fabs(zk*dpdz + phi0 - n*6.2831 - phik);
        if (dPhi > mindPhi) continue;
        mindPhi = dPhi;

        float pi    = phii;
        float pj    = phij + il*6.2831;
        float pk    = phik + n*6.2831;
        float a     = (3.0*(zi*pi+zj*pj+zk*pk)-(pi+pj+pk)*(zi+zj+zk))/(3.0*(zi*zi+zj*zj+zk*zk)-(zi+zj+zk)*(zi+zj+zk));
        float b     = ((pi+pj+pk) - a*(zi+zj+zk))/3.0;
        helix.dpdz_ = a;
        helix.phi0_ = b;
     }
  }



  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::initDpdzScan(const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitIdx)
  {
    float    dpdzBest(0),phi0Best(0);
    unsigned nshBest(0);

    size_t stride(1);
    for (size_t is=0; is < hitIdx.size(); is+=stride){
      for (size_t js=is+1;js < hitIdx.size(); js+=stride){
        auto  i    = hitIdx[is];
        auto  j    = hitIdx[js];

        float dp   = phiCcol[j] - phiCcol[i];
        float dz   = chcol[j].pos().z() - chcol[i].pos().z();
        if (fabs(dp) < 0.1) continue;
        if (fabs(dz) < 500) continue;

        int ilmax = int( dz*maxDpdz_/6.2831 - dp/6.2831);
        int ilmin = int(-dz*maxDpdz_/6.2831 - dp/6.2831);
        for (int il=ilmin;il<=ilmax;++il){
           float dpdz    = (dp + il*6.2831)/dz;
           float pOffset = dpdz > 0 ? -3.14159 : 3.14159;
           if (fabs(dpdz) > maxDpdz_ || fabs(dpdz) < minDpdz_) continue;

           unsigned nsh(0);
           float    phi0(phiCcol[i] - chcol[i].pos().z()*dpdz);
           for (auto k : helix.hits_){
              int   n    = int((chcol[k].pos().z()*dpdz + phi0 - pOffset)/6.2831);
              float dPhi = fabs(chcol[k].pos().z()*dpdz + phi0 - n*6.2831 - phiCcol[k]);
              if (dPhi < maxDphi_) {nsh+=chcol[k].nStrawHits();}
           }

           if (nsh < nshBest ) continue;
           dpdzBest  = dpdz;
           phi0Best  = phi0;
           nshBest   = nsh;

           if (nshBest == helix.nStrawHits_) break;
        }
      }
    }

    helix.dpdz_ = dpdzBest;
    helix.phi0_ = phi0Best;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::cleanDpdzHits(const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix)
  {
    for (auto it = helix.hits_.begin(); it != helix.hits_.end(); )
    {
      int   n      = int((chcol[*it].pos().z()*helix.dpdz_ + helix.phi0_ - helix.pOffset())/6.2831);
      float dPhi   = fabs(chcol[*it].pos().z()*helix.dpdz_ + helix.phi0_ - n*6.2831 - phiCcol[*it]);
      if (dPhi > maxDphi_) {helix.nStrawHits_ -= chcol[*it].nStrawHits();  it = helix.hits_.erase(it); }
      else                 {++it;}
    }
  }


 //---------------------------------------------------------------------------------------------------------------------------
 void RobustMultiHelixFinder2::fillPhiZFitter(const ComboHitCollection& chcol, const std::vector<float>& phiCcol, const CandHelix& helix, LSFitter& zphiFitter){

    zphiFitter.clear();
    for (const auto& ih : helix.hits_){
       const auto& ch = chcol[ih];
       int   n        = int((ch.pos().z()*helix.dpdz_ + helix.phi0_ - helix.pOffset())/6.2831);
       float dPhi     = ch.pos().z()*helix.dpdz_ + helix.phi0_ - n*6.2831 - phiCcol[ih];

       if      (dPhi > 3.1415)  {dPhi -= 6.2831; ++n;}
       else if (dPhi < -3.1415) {dPhi += 6.2831; --n;}

       if (fabs(dPhi) < maxDphi_) {zphiFitter.add(ch.pos().z(), phiCcol[ih] + n*6.2831,ch.nStrawHits());}
    }
 }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::dpdzFit(const ComboHitCollection& chcol, const std::vector<float>& phiCcol, CandHelix& helix)
  {
    LSFitter zphiFitter;
    fillPhiZFitter(chcol, phiCcol, helix, zphiFitter);

    helix.dpdz_ = zphiFitter.fa();
    helix.phi0_ = zphiFitter.fb();
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixRecoverDist(const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit){

    helix.hits_.clear();
    float previousValue(1e6);

    for (auto k : hitsInit){
        const auto& ch = chcol[k];
        float dist2 = dist2Helix(ch,helix);
        if (dist2 > maxDist2Recover_) continue;

        if (helix.hits_.empty() || ch.strawId().uniquePanel() != chcol[helix.hits_.back()].strawId().uniquePanel()) {
           helix.hits_.emplace_back(k);
           previousValue = dist2;
           continue;
        }

        if (dist2 < previousValue){
             previousValue = dist2;
             helix.hits_.back() = k;
        }
    }

    helix.nStrawHits_=0;
    for (const auto& ih : helix.hits_) helix.nStrawHits_ += chcol[ih].nStrawHits();
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixRecoverChi2(const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit){


    helix.hits_.clear();
    float previousValue(1e6);

    for (auto k : hitsInit){
        const auto& ch = chcol[k];
        float chi2 = chi2XYHelix(ch,helix);
        if (chi2 > maxChi2Recover_) continue;

        if (helix.hits_.empty() || ch.strawId().uniquePanel() != chcol[helix.hits_.back()].strawId().uniquePanel()) {
           helix.hits_.emplace_back(k);
           previousValue = chi2;
           continue;
        }

        if (chi2 < previousValue){
             previousValue = chi2;
             helix.hits_.back() = k;
        }
    }

    helix.nStrawHits_=0;
    for (const auto& ih : helix.hits_) helix.nStrawHits_ += chcol[ih].nStrawHits();
  }

  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixRecover(const ComboHitCollection& chcol, CandHelix& helix, const std::vector<StrawHitIndex>& hitsInit)
  {
    if (useDist2_) helixRecoverDist(chcol,helix,hitsInit);
    else          helixRecoverChi2(chcol,helix,hitsInit);
  }




  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixFilterDist(const ComboHitCollection& chcol, CandHelix& helix)
  {
    for (auto it = helix.hits_.begin(); it != helix.hits_.end(); )
    {
      if (dist2Helix(chcol[*it],helix) <  maxDist2Filter_) {++it;}
      else {helix.nStrawHits_ -= chcol[*it].nStrawHits();  it = helix.hits_.erase(it); }
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixFilterChi2(const ComboHitCollection& chcol, CandHelix& helix)
  {
    for (auto it = helix.hits_.begin(); it != helix.hits_.end(); )
    {
      if (chi2XYHelix(chcol[*it],helix) <  maxChi2Filter_) {++it;}
      else {helix.nStrawHits_ -= chcol[*it].nStrawHits();  it = helix.hits_.erase(it); }
    }
  }

  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::helixFilter(const ComboHitCollection& chcol, CandHelix& helix)
  {
    if (useDist2_) return helixFilterDist(chcol,helix);
    return helixFilterChi2(chcol,helix);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  float RobustMultiHelixFinder2::helixSumDist2(const ComboHitCollection& chcol, CandHelix& helix)
  {
    float sum(0);
    for (const auto& k : helix.hits_) sum += dist2Helix(chcol[k],helix)*chcol[k].nStrawHits();
    return sum;
  }

  //---------------------------------------------------------------------------------------------------------------------------
  float RobustMultiHelixFinder2::helixSumChi2(const ComboHitCollection& chcol, CandHelix& helix)
  {
    float sum(0);
    for (const auto k : helix.hits_) sum += chi2XYHelix(chcol[k],helix)*chcol[k].nStrawHits();
    return sum;
  }

  //---------------------------------------------------------------------------------------------------------------------------
  float RobustMultiHelixFinder2::helixSumResid(const ComboHitCollection& chcol, CandHelix& helix)
  {
    if (useDist2_) return helixSumDist2(chcol,helix);
    return helixSumChi2(chcol,helix);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder2::refineXY(const ComboHitCollection& chcol, std::vector<float>& phiCcol, CandHelix& helix)
  {
    std::vector<float> val;
    std::vector<float> dispX{0.0f, -0.1f, 0.1f,  0.0f, 0.0f, 0.1f};
    std::vector<float> dispY{0.0f,  0.0f, 0.0f, -0.1f, 0.1f, 0.1f};

    LSFitter zphiFitter;
    float sumWeights(helix.hits_.size());

    for (size_t i=0;i<dispX.size();++i){
       float centerX = helix.x_ + dispX[i];
       float centerY = helix.y_ + dispY[i];

       float radius(0);
       for (const auto k : helix.hits_){
         float dx = (chcol[k].pos().x()-centerX);
         float dy = (chcol[k].pos().y()-centerY);
         radius += (dx*dx+dy*dy);
       }
       radius = sqrtf(radius/sumWeights);



       CandHelix thelix(helix);
       thelix.x_ = centerX;
       thelix.y_ = centerY;
       thelix.r_ = radius;


       for (const auto& ich : thelix.hits_) phiCcol[ich] = polyAtan2(chcol[ich].pos().y() - thelix.y_, chcol[ich].pos().x() - thelix.x_);
       fillPhiZFitter(chcol,phiCcol,helix,zphiFitter);
       if (zphiFitter.isValid()) {
         float fa = zphiFitter.fa();
         float fb = zphiFitter.fb();
         thelix.dpdz_ = fa;
         thelix.phi0_ = fb;
         val.push_back(helixSumResid(chcol,thelix));
       } else {
         val.push_back(0);
       }

       //Do it without recalculating dpdz to speed up?
      //val.push_back(helixSumResid(chcol,thelix));
    }


    //Calculate Hessian entries
    float d2fx  = (val[2]+val[1]-2*val[0])/0.01;
    float d2fy  = (val[4]+val[3]-2*val[0])/0.01;
    float d2fxy = (val[5]-val[2]-val[4]+val[0])/0.01;
    float dfx   = (val[2]-val[0])/0.1;
    float dfy   = (val[4]-val[0])/0.1;
    float det   = d2fx*d2fy-d2fxy*d2fxy;
    if (fabs(det) <1e-6) return;

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
      float dx = (chcol[idx].pos().x()-helix.x_);
      float dy = (chcol[idx].pos().y()-helix.y_);
      radius  += (dx*dx+dy*dy);
    }
    helix.r_ = sqrtf(radius/sumWeights);

    for (auto ich : helix.hits_) phiCcol[ich] = polyAtan2(chcol[ich].pos().y() - helix.y_, chcol[ich].pos().x() - helix.x_);
    fillPhiZFitter(chcol,phiCcol,helix,zphiFitter);
    helix.dpdz_ = zphiFitter.fa();
    helix.phi0_ = zphiFitter.fb();

    return;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Estimate time uncertainty - variance of hit with simple nHit weight
  float RobustMultiHelixFinder2::timeUncertainty(const ComboHitCollection& chcol, const CandHelix& helix)
  {
    return 0;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Print me!
  void RobustMultiHelixFinder2::printHelix(const HelixSeed& helix)
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
  void RobustMultiHelixFinder2::fillDiag(mapHelix& helcols, const ComboHitCollection& chcol)
  {
    data_.iev_ = iev_;

    data_.Nch_ = chcol.size();
    for (unsigned ich=0; ich<chcol.size();++ich)
    {
       data_.chTime_[ich] = chcol[ich].correctedTime();
       data_.chX_[ich]    = chcol[ich].pos().x();
       data_.chY_[ich]    = chcol[ich].pos().y();
       data_.chZ_[ich]    = chcol[ich].pos().z();
       data_.chRad_[ich]  = sqrt(data_.chX_[ich]*data_.chX_[ich]+data_.chY_[ich]*data_.chY_[ich]);
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
      data_.helchi2_[data_.Nhel_] = helix.helix().chi2dZPhi();
      data_.heldzdt_[data_.Nhel_] = helix.helix().chi2dXY();
      data_.helnhi_[data_.Nhel_]  = hhits.size();
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

DEFINE_ART_MODULE(mu2e::RobustMultiHelixFinder2);
