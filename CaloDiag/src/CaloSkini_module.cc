//
// An EDAnalyzer module that reads back the hits created by the calorimeter and produces an ntuple
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"

#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include "TDirectory.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"


namespace
{
   constexpr int ntupLen = 16384;
}


namespace mu2e {

    class CaloSkini : public art::EDAnalyzer {

      public:
      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;
        fhicl::Atom<art::InputTag>     caloHitCollection     { Name("caloHitCollection"),      Comment("Calo Hit collection name") };
        fhicl::Atom<art::InputTag>     caloClusterCollection { Name("caloClusterCollection"),  Comment("Calo cluster collection name") };
        fhicl::Atom<int>               diagLevel             { Name("diagLevel"),              Comment("Diag Level"),0 };
      };

      explicit CaloSkini(const art::EDAnalyzer::Table<Config>& config);
      virtual ~CaloSkini() {}

      virtual void beginJob();
      virtual void endJob() {};
      virtual void analyze(const art::Event& e);


      private:
        art::InputTag         virtualDetectorTag_;
        art::InputTag         caloHitTag_;
        art::InputTag         caloClusterTag_;
        art::InputTag         caloHitTruthTag_;
        art::InputTag         caloClusterTruthTag_;
        int                   diagLevel_;
        int                   nProcess_;


        TH1F *hcryE_,*hcryT_,*hcryX_,*hcryY_,*hcryZ_;
        TH1F *hcluE_,*hcluT_,*hcluX_,*hcluY_,*hcluZ_,*hcluE_1Et,*hcluE_1E9,*hcluE_1E25,*hcluE_F;
        TH2F *hxy_,*hECE,*hCryEEMC_,*hCryTTMC_,*hCluEEMC_,*hCluTTMC_,*hCryEEMC2_;

        TTree* Ntup_;
        int   _evt,_run;

        int   nHits_,cryId_[ntupLen],crySectionId_[ntupLen],crySimIdx_[ntupLen],crySimLen_[ntupLen];
        float cryEtot_,cryTime_[ntupLen],cryEdep_[ntupLen],cryEdepErr_[ntupLen],cryPosX_[ntupLen],cryPosY_[ntupLen],cryPosZ_[ntupLen],_cryLeak[ntupLen];

        int   nCluster_,nCluSim_,cluNcrys_[ntupLen];
        float cluEnergy_[ntupLen],cluEnergyErr_[ntupLen],cluTime_[ntupLen],cluTimeErr_[ntupLen],cluCogX_[ntupLen],cluCogY_[ntupLen],
        cluCogZ_[ntupLen],cluE1_[ntupLen],cluE9_[ntupLen],cluE25_[ntupLen],cluSecMom_[ntupLen];
        int   cluSplit_[ntupLen],cluConv_[ntupLen],cluSimIdx_[ntupLen],cluSimLen_[ntupLen];
        std::vector<std::vector<int> > cluList_;

      };


    CaloSkini::CaloSkini(const art::EDAnalyzer::Table<Config>& config) :
      EDAnalyzer{config},
      caloHitTag_         (config().caloHitCollection()),
      caloClusterTag_     (config().caloClusterCollection()),
      diagLevel_          (config().diagLevel()),
      nProcess_(0),
      Ntup_(0)
      {}

    void CaloSkini::beginJob(){

      art::ServiceHandle<art::TFileService> tfs;

      Ntup_  = tfs->make<TTree>("Calo", "Calo");

      Ntup_->Branch("evt",          &_evt ,         "evt/I");
      Ntup_->Branch("run",          &_run ,         "run/I");
      Ntup_->Branch("cryEtot",      &cryEtot_ ,     "cryEtot/F");

      Ntup_->Branch("nCry",         &nHits_ ,       "nCry/I");
      Ntup_->Branch("cryId",        &cryId_ ,       "cryId[nCry]/I");
      Ntup_->Branch("crySectionId", &crySectionId_, "crySectionId[nCry]/I");
      Ntup_->Branch("cryPosX",      &cryPosX_ ,     "cryPosX[nCry]/F");
      Ntup_->Branch("cryPosY",      &cryPosY_ ,     "cryPosY[nCry]/F");
      Ntup_->Branch("cryPosZ",      &cryPosZ_ ,     "cryPosZ[nCry]/F");
      Ntup_->Branch("cryEdep",      &cryEdep_ ,     "cryEdep[nCry]/F");
      Ntup_->Branch("cryEdepErr",   &cryEdepErr_ ,  "cryEdepErr[nCry]/F");
      Ntup_->Branch("cryTime",      &cryTime_ ,     "cryTime[nCry]/F");

      Ntup_->Branch("nCluster",     &nCluster_ ,    "nCluster/I");
      Ntup_->Branch("cluEnergy",    &cluEnergy_ ,   "cluEnergy[nCluster]/F");
      Ntup_->Branch("cluEnergyErr", &cluEnergyErr_ ,"cluEnergyErr[nCluster]/F");
      Ntup_->Branch("cluTime",      &cluTime_ ,     "cluTime[nCluster]/F");
      Ntup_->Branch("cluTimeErr",   &cluTimeErr_ ,  "cluTimeErr[nCluster]/F");
      Ntup_->Branch("cluCogX",      &cluCogX_ ,     "cluCogX[nCluster]/F");
      Ntup_->Branch("cluCogY",      &cluCogY_ ,     "cluCogY[nCluster]/F");
      Ntup_->Branch("cluCogZ",      &cluCogZ_ ,     "cluCogZ[nCluster]/F");
      Ntup_->Branch("cluNcrys",     &cluNcrys_ ,    "cluNcrys[nCluster]/I");
      Ntup_->Branch("cluE1",        &cluE1_ ,       "cluE1[nCluster]/F");
      Ntup_->Branch("cluE9",        &cluE9_ ,       "cluE9[nCluster]/F");
      Ntup_->Branch("cluE25",       &cluE25_ ,      "cluE25[nCluster]/F");
      Ntup_->Branch("cluSecMom",    &cluSecMom_ ,   "cluSecMom[nCluster]/F");
      Ntup_->Branch("cluSplit",     &cluSplit_ ,    "cluSplit[nCluster]/I");
      Ntup_->Branch("cluList",      &cluList_);

    }

    void CaloSkini::analyze(const art::Event& event){
      ++nProcess_;
      if (nProcess_%10==0 && diagLevel_ > 0) std::cout<<"Processing event from CaloSkini =  "<<nProcess_ <<std::endl;

      //Handle to the calorimeter
      art::ServiceHandle<GeometryService> geom;
      if (!geom->hasElement<Calorimeter>() ) return;
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());

      //Calorimeter crystal hits (average from readouts)
      art::Handle<CaloHitCollection> CaloHitsHandle;
      event.getByLabel(caloHitTag_, CaloHitsHandle);
      const CaloHitCollection& CaloHits(*CaloHitsHandle);

      //Calorimeter clusters
      art::Handle<CaloClusterCollection> caloClustersHandle;
      event.getByLabel(caloClusterTag_, caloClustersHandle);
      const CaloClusterCollection& caloClusters(*caloClustersHandle);


      _evt = event.id().event();
      _run = event.run();

      if (diagLevel_ == 3){std::cout << "processing event in calo_example " << nProcess_ << " run and event  = " << _run << " " << _evt << std::endl;}

      //--------------------------  Do calorimeter hits --------------------------------

      nHits_ = 0;

      for (unsigned int ic=0; ic<CaloHits.size();++ic)
      {
        const CaloHit& hit            = CaloHits.at(ic);
        int diskId                    = cal.crystal(hit.crystalID()).diskID();
        CLHEP::Hep3Vector crystalPos  = cal.geomUtil().mu2eToDiskFF(diskId,cal.crystal(hit.crystalID()).position());  //in disk FF frame

        cryId_[nHits_]        = hit.crystalID();
        crySectionId_[nHits_] = diskId;
        cryEdep_[nHits_]      = hit.energyDep();
        cryEdepErr_[nHits_]   = hit.energyDepErr();
        cryTime_[nHits_]      = hit.time();
        cryPosX_[nHits_]      = crystalPos.x();
        cryPosY_[nHits_]      = crystalPos.y();
        cryPosZ_[nHits_]      = crystalPos.z();
        ++nHits_;
      }

      //--------------------------  Do clusters --------------------------------
      nCluster_ = 0;
      cluList_.clear();
      for (unsigned int ic=0; ic<caloClusters.size();++ic)
      {
        const CaloCluster& cluster = caloClusters.at(ic);

        cluEnergy_[nCluster_]    = cluster.energyDep();
        cluEnergyErr_[nCluster_] = cluster.energyDepErr();
        cluTime_[nCluster_]      = cluster.time();
        cluTimeErr_[nCluster_]   = cluster.timeErr();
        cluNcrys_[nCluster_]     = cluster.size();
        cluCogX_[nCluster_]      = cluster.cog3Vector().x(); //in disk FF frame
        cluCogY_[nCluster_]      = cluster.cog3Vector().y();
        cluCogZ_[nCluster_]      = cluster.cog3Vector().z();
        ++nCluster_;
      }


      Ntup_->Fill();
}



}

DEFINE_ART_MODULE(mu2e::CaloSkini)
