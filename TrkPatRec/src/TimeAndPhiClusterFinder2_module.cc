//
// Tracker time / phi / z cluster finder
//
// Original author B. Echenard
//

#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"

#include <algorithm>
#include <numeric>


namespace
{
   struct TimePhiCandidate
   {
       TimePhiCandidate() {strawIdx_.reserve(64);}

       TimePhiCandidate(unsigned nsh, const std::vector<StrawHitIndex>& strawIdx, const art::Ptr<mu2e::CaloCluster>& caloCluster) :
         nsh_(nsh), strawIdx_(strawIdx), caloCluster_(caloCluster)
       {}

       unsigned                    nsh_      = 0;
       std::vector<StrawHitIndex>  strawIdx_ = {};
       art::Ptr<mu2e::CaloCluster> caloCluster_;
       float                       a_ = 0;
       float                       b_ = 0;
   };

   typedef std::vector<TimePhiCandidate> TimePhiCandidateCollection;
}






namespace mu2e {

  class TimeAndPhiClusterFinder2 : public art::EDProducer
  {
    public:

      using Config_types = TimeAndPhiClusterFinderTypes::Config ;
      using Data_types   = TimeAndPhiClusterFinderTypes::Data_t ;
      struct Config
      {
          using Name    = fhicl::Name;
          using Comment = fhicl::Comment;
          fhicl::Atom<art::InputTag>              comboHitCollection     {Name("ComboHitCollection"),     Comment("ComboHit collection {Name") };
          fhicl::Atom<art::InputTag>              caloClusterCollection  {Name("CaloClusterCollection"),  Comment("Calo cluster collection {Name") };
          fhicl::Table<MVATools::Config>          MVATime                {Name("MVATime"),                Comment("MVA for time cluster cleaning") };
          fhicl::Sequence<std::string>            hsel                   {Name("HitSelectionBits"),       Comment("HitSelectionBits") };
          fhicl::Sequence<std::string>            hbkg                   {Name("HitBackgroundBits"),      Comment("HitBackgroundBits") };
          fhicl::Atom<bool>                       usecc                  {Name("UseCaloCluster"),         Comment("Use calorimeter cluster") };
          fhicl::Atom<unsigned>                   minNSHits              {Name("MinNSHits"),              Comment("Minimum number of hits for cluster") };
          fhicl::Atom<float>                      tbin                   {Name("Tbin"),                   Comment("Time histogram bin width") };
          fhicl::Atom<unsigned>                   minTimeYbin            {Name("MinTimeYbin"),            Comment("Minimum number of bins to start recording max for scanning algo") };
          fhicl::Atom<float>                      maxTimeDT              {Name("MaxTimeDT"),              Comment("Max time difference for hits in cluster") };
          fhicl::Atom<float>                      maxFitDT               {Name("MaxFitDT"),               Comment("Max time difference for hits in cluster with linear regression") };
          fhicl::Atom<float>                      maxDeltaPhi            {Name("MaxDeltaPhi"),            Comment("Max delta phi between consecutive hits in cluster") };
          fhicl::Atom<int>                        diag                   {Name("Diag"),                   Comment("Diag level"), 0 };
          fhicl::Table<Config_types>              diagPlugin             {Name("DiagPlugin"),             Comment("Diag Plugin config")};
      };


      explicit TimeAndPhiClusterFinder2(const art::EDProducer::Table<Config>& config);
      ~TimeAndPhiClusterFinder2() = default;

      void beginJob() override;
      void produce(art::Event& e) override;


    private:
      int                                             iev_;
      const art::ProductToken<ComboHitCollection>     chToken_;
      const art::ProductToken<CaloClusterCollection>  ccToken_;
      MVATools                                        MVATime_;
      StrawHitFlag                                    hsel_;
      StrawHitFlag                                    hbkg_;
      bool                                            usecc_;
      unsigned                                        minNSHits_;
      float                                           tbin_;
      unsigned                                        minTimeYbin_;
      float                                           maxTimeDT_;
      float                                           maxFitDT_;
      float                                           maxDeltaPhi_;
      const Calorimeter*                              cal_;
      int                                             diag_;
      std::unique_ptr<ModuleHistToolBase>             diagTool_;
      Data_types                                      data_;


      void findClusters     (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                             TimeClusterCollection& tccol);
      void findTimePeaks    (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                             TimePhiCandidateCollection& timeCandidates);
      void dzdt_fit         (const ComboHitCollection& chcol,               TimePhiCandidate& tc);
      void addCaloPtr       (const art::Handle<CaloClusterCollection>& ccH, TimePhiCandidate& tc);
      void findPhiPeaks     (const ComboHitCollection& chcol,               TimePhiCandidateCollection& timeCandidates,
                             TimePhiCandidateCollection& phiCandidates);
      void calculateMean    (const ComboHitCollection& chcol,               TimeCluster& tc);
      void fillTCcol        (const TimePhiCandidateCollection& candidates,  const ComboHitCollection& chcol,
                             TimeClusterCollection& tccol);
      void fillDiag         (const TimePhiCandidateCollection& candidates,  const ComboHitCollection& chcol,
                             const TimeClusterCollection& tccol);

      void timeClusterFromCalo(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol, TimePhiCandidateCollection& tccol);
  };




  TimeAndPhiClusterFinder2::TimeAndPhiClusterFinder2(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    iev_(0),
    chToken_             {consumes<ComboHitCollection>      (config().comboHitCollection())     },
    ccToken_             {mayConsume<CaloClusterCollection> (config().caloClusterCollection())  },
    MVATime_             (config().MVATime()),
    hsel_                (config().hsel()),
    hbkg_                (config().hbkg()),
    usecc_               (config().usecc()),
    minNSHits_           (config().minNSHits()),
    tbin_                (config().tbin()),
    minTimeYbin_         (config().minTimeYbin()),
    maxTimeDT_           (config().maxTimeDT()),
    maxFitDT_            (config().maxFitDT()),
    maxDeltaPhi_         (config().maxDeltaPhi()),
    diag_                (config().diag()),
    diagTool_(),
    data_()
    {
       produces<TimeClusterCollection>();
       if (diag_) diagTool_ = art::make_tool<ModuleHistToolBase>(config().diagPlugin," ");
    }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder2::beginJob()
  {
      MVATime_.initMVA();

      if (diag_){
         art::ServiceHandle<art::TFileService> tfs;
         diagTool_->bookHistograms(tfs);
         MVATime_.showMVA();
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder2::produce(art::Event & event )
  {
      iev_ = event.id().event();
//std::cout<<"Event "<<iev_<<std::endl;
      art::Handle<CaloClusterCollection> ccH{};
      if (usecc_) ccH = event.getHandle<CaloClusterCollection>(ccToken_);

      const auto& chH = event.getValidHandle(chToken_);
      const auto& chcol(*chH);

      if (diag_) {data_.reset(); data_.event_=&event; data_.chcol_ = &chcol;}

      std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
      findClusters(ccH, chcol, *tccol);

      if (diag_) diagTool_->fillHistograms(&data_);
      event.put(std::move(tccol));
  }


  //-------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder2::findClusters(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                             TimeClusterCollection& tccol)
  {
      cal_ = &(*GeomHandle<Calorimeter>()); //yes this is clunky....

      // Find time peaks
      std::vector<TimePhiCandidate> timeCandidates;
      timeCandidates.reserve(64);
      findTimePeaks(ccH, chcol, timeCandidates);


//if (1<0)
 timeClusterFromCalo(ccH,chcol,timeCandidates);


      // Refine time peaks and create split phi clusters
      std::vector<TimePhiCandidate> phiCandidates;
      phiCandidates.reserve(64);
      findPhiPeaks(chcol, timeCandidates, phiCandidates);

      // Finally create the timeClusters
      tccol.reserve(64);
      fillTCcol(phiCandidates,chcol,tccol);

      if (diag_) fillDiag(timeCandidates, chcol, tccol);
  }


  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in time distribution
  void TimeAndPhiClusterFinder2::findTimePeaks(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                               TimePhiCandidateCollection& timeCandidates)
  {
      // Select good hits and sort them by time
      std::vector<unsigned> chGood;
      chGood.reserve(chcol.size());
      for (size_t ich=0; ich<chcol.size();++ich){
        if (chcol[ich].flag().hasAllProperties(hsel_) && !chcol[ich].flag().hasAnyProperty(hbkg_)) chGood.emplace_back(ich);
      }
      auto sortFcn = [&chcol](const auto& i1, const auto& i2){return chcol[i1].correctedTime() < chcol[i2].correctedTime();};
      sort(chGood.begin(),chGood.end(),sortFcn);

      // Fill the time histogram with hit corrected times, add calo if needed
      float tmax     = chcol[chGood.back()].correctedTime() + tbin_;
      float tmin     = std::max(0.0f, chcol[chGood.front()].correctedTime() - tbin_);
      unsigned nbins = unsigned((tmax-tmin)/tbin_);

      std::vector<unsigned> timeHist(nbins,0), timeIdx(nbins,0);
      std::iota(timeIdx.begin(),timeIdx.end(),0);
      for (const auto ich : chGood){
          unsigned ibin = unsigned((chcol[ich].correctedTime()-tmin)/tbin_);
          timeHist[ibin] += chcol[ich].nStrawHits();
      }

      std::sort(timeIdx.begin(),timeIdx.end(),[&timeHist](unsigned i, unsigned j){return timeHist[i] > timeHist[j];});


      //std::set<int> vetoBins;
      std::vector<float> timePeakMin, timePeakMax;
      for (unsigned i=0;i<timeIdx.size();++i){
         if (timeHist[timeIdx[i]] < minTimeYbin_) break;
         float tpeak = tmin + timeIdx[i]*tbin_ + 0.5*tbin_;
         timePeakMin.push_back(tpeak - 0.5*tbin_ - tbin_);
         timePeakMax.push_back(tpeak + 0.5*tbin_ + tbin_);
         //vetoBins.insert(timeIdx[i]-1);
         //vetoBins.insert(timeIdx[i]+1);
      }

      //Create time cluster candidates
      std::vector<bool> usedHit(chcol.size(),false);
      for (size_t i=0;i<timePeakMin.size();++i){
           TimePhiCandidate tc;
           for (const auto& ich : chGood){
               //if (usedHit[ich]) continue;
               float time = chcol[ich].correctedTime();
               if (time < timePeakMin[i]) continue;
               if (time > timePeakMax[i]) break;
               tc.strawIdx_.emplace_back(ich);
               tc.nsh_ += chcol.at(ich).nStrawHits();
           }
           if (tc.strawIdx_.size()<2) continue;

           //calculate dz/dt with linear regression
           dzdt_fit(chcol,tc);

           //reassociate hits based on the fit
           tc.strawIdx_.clear();
           tc.nsh_ = 0;
           for (const auto& ich : chGood){
              if (usedHit[ich]) continue;
              float time = chcol[ich].correctedTime();
              float dt = time - tc.a_*chcol[ich].pos().z()-tc.b_;
              if (abs(dt) > maxFitDT_)   continue;
              if (dt      > 2*maxFitDT_) break;

              tc.nsh_ += chcol.at(ich).nStrawHits();
              tc.strawIdx_.emplace_back(ich);
              usedHit[ich] = true;
          }

          if (usecc_) addCaloPtr(ccH,tc);
          if (tc.nsh_ < minNSHits_) continue;
          timeCandidates.emplace_back(std::move(tc));
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Simple linear t-z regression, weighted by number of hits
  void TimeAndPhiClusterFinder2::dzdt_fit(const ComboHitCollection& chcol, TimePhiCandidate& tc)
  {
      float sz(0),sz2(0),st(0),szt(0),sn(0);
      for (auto ich : tc.strawIdx_){
         const auto& ch = chcol[ich];
         int   nsh = ch.nStrawHits();
         float z   = ch.pos().z();
         float t   = ch.correctedTime();
         sn  += nsh;
         sz  += z*nsh;
         sz2 += z*z*nsh;
         st  += t*nsh;
         szt += t*z*nsh;
      }
      float den = sn*sz2 - sz*sz;
      if (abs(den) < 1e-6) den = 1e-6;
      tc.a_ = (sn*szt -sz*st)/den;
      tc.b_ = (st*sz2 - sz*szt)/den;
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder2::addCaloPtr(const art::Handle<CaloClusterCollection>& ccH, TimePhiCandidate& tc){

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;

      unsigned icalo(cccol.size());
      float maxEcalo(0);
      for (size_t i=0; i < cccol.size(); ++i){
          const auto& cluster = cccol[i];
          const auto trackerPos = cal_->geomUtil().mu2eToTracker(cal_->geomUtil().diskFFToMu2e(cluster.diskID(),cluster.cog3Vector()));
          float dt = abs(cccol[i].time() - tc.a_*trackerPos.z()-tc.b_);
          if (dt > maxTimeDT_ || cccol[i].energyDep() < maxEcalo) continue;
          icalo=i;
          maxEcalo = cccol[i].energyDep();
      }
      if (icalo<cccol.size()) tc.caloCluster_ = art::Ptr<CaloCluster>(ccH,icalo);
   }




  //--------------------------------------------------------------------------------------------------------------
  //Filtering and split clusters in phi
  void TimeAndPhiClusterFinder2::findPhiPeaks(const ComboHitCollection& chcol, TimePhiCandidateCollection& timeCandidates,
                                             TimePhiCandidateCollection& phiCandidates)
  {
       // Cache phi value of comboHits for efficiency, this is recalculated every time the fcn is called
       std::vector<float> chPhi;
       chPhi.reserve(chcol.size());
       for (const auto& ch : chcol) chPhi.emplace_back(polyAtan2(ch.pos().y(),ch.pos().x()));

       for (auto& tc : timeCandidates){

           // Sort hits in ascending azimuthal angle
           auto& hits = tc.strawIdx_;
           sort(hits.begin(),hits.end(),[&chPhi](const int a, const int b){return chPhi[a]<chPhi[b];});

//std::cout<<"Begin dphi"<<std::endl;
//for (auto h : hits) std::cout<<h<<" ";
//std::cout<<std::endl;

           // Start by fast forwarding to first large gap in phi.
           // note: (idx+vector_size)%vector_size maps negative and overflow indices -> valid vector indices
           size_t idx(0),vsize(hits.size());
           while (idx < vsize){
               float deltaPhi = chPhi[hits[idx]] - chPhi[hits[(idx+vsize-1)%vsize]];
               if (deltaPhi < 0) deltaPhi += 6.2832;
               if (deltaPhi > maxDeltaPhi_) break;
               ++idx;
           }
//std::cout<<"Idx start "<<idx<<std::endl;
           // Find sequences of contiguous phi hits - first (last) point separated by an phi angle greater than threshold
           // with previous (next) point - flag them for removal if their number is below a threshold
           // Only save split clusters smaller than the original cluster to avoid duplicates
           size_t idxMax(idx+vsize), nStrawHits(0);
           std::vector<StrawHitIndex> buffer;
           while (idx < idxMax){
               buffer.emplace_back(hits[idx%vsize]);
               nStrawHits += chcol[hits[idx%vsize]].nStrawHits();

               float deltaPhi = chPhi[hits[(idx+1)%vsize]]-chPhi[hits[idx%vsize]];
               if (deltaPhi < 0) deltaPhi += 6.2832;

               if (deltaPhi< maxDeltaPhi_ && idx+1!=idxMax ){++idx;continue;}

//std::cout<<"Empty buffer"<<std::endl;
//for (auto h : buffer) std::cout<<h<<" ";
//std::cout<<std::endl;

               if (nStrawHits >= minNSHits_) {
                   phiCandidates.emplace_back(TimePhiCandidate(nStrawHits,buffer,tc.caloCluster_));
               }

               buffer.clear();
               nStrawHits = 0;
               ++idx;
           }
       }
   }



  //--------------------------------------------------------------------------------------------------------------
  // Caluclate cluster variables, find duplicates and fill collections
  void TimeAndPhiClusterFinder2::calculateMean(const ComboHitCollection& chcol, TimeCluster& tc)
  {
      if (tc._strawHitIdxs.empty()){tc._pos = XYZVectorF(0,0,0);tc._t0._t0=0;tc._t0._t0err=0;tc._nsh=0; return;};

      tc._nsh = 0;
      float tacc(0),tacc2(0),xacc(0),yacc(0),zacc(0),weight(0);
      for (const auto& ish :tc._strawHitIdxs){
          const ComboHit& ch = chcol[ish];
          tc._nsh += ch.nStrawHits();

          float htime = chcol[ish].correctedTime();
          float hwt   = ch.nStrawHits();

          weight += hwt;
          tacc   += htime*hwt;
          tacc2  += htime*htime*hwt;
          xacc   += ch.pos().x()*hwt;
          yacc   += ch.pos().y()*hwt;
          zacc   += ch.pos().z()*hwt;
      }

      tacc/=weight;
      tacc2/=weight;
      xacc/=weight;
      yacc/=weight;
      zacc/=weight;

      tc._t0._t0    = tacc;
      tc._t0._t0err = sqrtf(tacc2-tacc*tacc);
      tc._pos       = XYZVectorF(xacc, yacc, zacc);
  }

  //----------------------------------------------------------ch----------------------------------------------------
  void TimeAndPhiClusterFinder2::fillTCcol(const TimePhiCandidateCollection& candidates, const ComboHitCollection& chcol,
                                          TimeClusterCollection& tccol)
  {
      for (const auto& cand : candidates){
          if (cand.nsh_ < minNSHits_) continue;
          TimeCluster tclu;
          tclu._caloCluster  = cand.caloCluster_;
          tclu._nsh          = cand.nsh_;
          tclu._strawHitIdxs = std::move(cand.strawIdx_);

          calculateMean(chcol, tclu);
          tccol.emplace_back(std::move(tclu));
      }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Diagnosis
  void TimeAndPhiClusterFinder2::fillDiag(const TimePhiCandidateCollection& timeCandidates, const ComboHitCollection& chcol,
                                         const TimeClusterCollection& tccol)
  {
     data_.iev_ = iev_;

     data_.Nch_ = chcol.size();
     for (unsigned ich=0; ich<chcol.size();++ich)
     {
        int   selFlag = (chcol[ich].flag().hasAllProperties(hsel_) && !chcol[ich].flag().hasAnyProperty(hbkg_)) ? 1 : 0;

        data_.chSel_[ich]  = selFlag;
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
        data_.chWDX_[ich]  = chcol[ich].wdir().x();
        data_.chWDY_[ich]  = chcol[ich].wdir().y();
     }

     int iclu1(0),ih1(0);
     for (const auto& tc : timeCandidates)
     {
        for (auto& ish : tc.strawIdx_)
        {
           data_.iclu1_[ih1]  = iclu1;
           data_.hitIdx1_[ih1]= ish;
           ++ih1;
        }
//std::cout<<"Cluster i "<<iclu1<<std::endl;
//for (const auto& h  : tc.strawIdx_) std::cout<<h<<" ";
//std::cout<<std::endl;
        ++iclu1;
     }
     data_.nhit1_=ih1;

     int iclu2(0),ih2(0);
     for (const auto& tc : tccol)
     {
        for (auto& ish : tc._strawHitIdxs)
        {
           data_.iclu2_[ih2]    = iclu2;
           data_.hitIdx2_[ih2]  = ish;
           ++ih2;
        }
//std::cout<<"Cluster i "<<iclu2<<std::endl;
//for (const auto& h  : tc._strawHitIdxs) std::cout<<h<<" ";
//std::cout<<std::endl;
        if (tc.caloCluster()){
          data_.calo2X_[iclu2] = tc.caloCluster()->cog3Vector().x();
          data_.calo2Y_[iclu2] = tc.caloCluster()->cog3Vector().y();
          data_.calo2Z_[iclu2] = tc.caloCluster()->cog3Vector().z();
          data_.calo2E_[iclu2] = tc.caloCluster()->energyDep();
          data_.calo2T_[iclu2] = tc.caloCluster()->time();
        } else {
          data_.calo2X_[iclu2]=data_.calo2Y_[iclu2]=data_.calo2Z_[iclu2]=data_.calo2E_[iclu2]=data_.calo2T_[iclu2]=0.0;
        }
        ++iclu2;
     }
     data_.nhit2_=ih2;
     data_.nclu2_=iclu2;
   }








void TimeAndPhiClusterFinder2::timeClusterFromCalo(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
TimePhiCandidateCollection& timeCandidates)
{

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;


      std::vector<float> chPhi;
      chPhi.reserve(chcol.size());
      for (const auto& ch : chcol) chPhi.emplace_back(polyAtan2(ch.pos().y(),ch.pos().x()));


      std::vector<unsigned> chGood;
      chGood.reserve(chcol.size());
      for (size_t ich=0; ich<chcol.size();++ich){
        if (chcol[ich].flag().hasAllProperties(hsel_)) chGood.emplace_back(ich);
      }
      auto sortFcn = [&chcol](const auto& i1, const auto& i2){return chcol[i1].correctedTime() < chcol[i2].correctedTime();};
      sort(chGood.begin(),chGood.end(),sortFcn);


      for (size_t icalo=0; icalo < cccol.size(); ++icalo){
          const auto& cluster = cccol[icalo];
          if (cluster.energyDep()<60) continue;

          const auto trackerPos = cal_->geomUtil().mu2eToTracker(cal_->geomUtil().diskFFToMu2e(cluster.diskID(),cluster.cog3Vector()));
          float caloZ   = trackerPos.z();
          float caloPhi = polyAtan2(trackerPos.y(),trackerPos.x());
          float caloT   = cluster.time();

/*
          int nsteps(20);

          int nhiMax(0);
          float fabest(-1),sdeltaMax(10000);
          for (int j=0;j<=nsteps;++j){
              int nhi(0);
              float fa(j*0.01/nsteps),sdt(0);
              for (const auto& ich : chGood){
                 float dt2 = abs(chcol[ich].correctedTime()-fa*(chcol[ich].pos().z()-caloZ)-caloT);
                 float dphi = std::min(6.2831f-fabs(chPhi[ich]-caloPhi), fabs(chPhi[ich]-caloPhi));
                 if (dt2<15 && fabs(dphi)<1.57) {nhi += chcol[ich].nStrawHits(); sdt += dt2;}
              }
              if (nhi > nhiMax || (nhi==nhiMax && sdt < sdeltaMax)) {nhiMax=nhi; fabest=fa; sdeltaMax = sdt;}
          }
          if (nhiMax < 15) return;

std::cout<<fabest<<std::endl;
*/

float fabest(0.005);

          TimePhiCandidate tc;
          for (const auto& ich : chGood){
             float dt2 = abs(chcol[ich].correctedTime()-fabest*(chcol[ich].pos().z()-caloZ)-caloT);
             float dphi = std::min(6.2831f-fabs(chPhi[ich]-caloPhi), fabs(chPhi[ich]-caloPhi));
             if (dt2>20 || fabs(dphi)>1.57) continue;
             tc.nsh_ += chcol.at(ich).nStrawHits();
             tc.strawIdx_.emplace_back(ich);
          }

          if (tc.nsh_ < 15) return;

          tc.caloCluster_ = art::Ptr<CaloCluster>(ccH,icalo);
          timeCandidates.emplace_back(std::move(tc));
      }
}


}

DEFINE_ART_MODULE(mu2e::TimeAndPhiClusterFinder2);
