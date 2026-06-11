//
// Create a compressed representation of Calorimeter StepPointMCs
//
// Basic idea: We recollect all StepPointMCS attached to a SimParticle ancestor and collapse them into CaloShowerStep objects.
// We call a SimParticle entering the calorimeter an Ancestore SimParticle. This ancestor will generate a shower of
// SimParticles and StepPointMcs in the crystal, which will be compressed
// At the end of the modules, all StepPointMCs can be dropped, as well as a large fraction of SimParticles
// with negligible loss of information.
//
// Note: if a SimParticle enters the calorimeter, generates a secondary SimParticle that hit another disk of the calorimeter
// (e.g. a photon leaks from the first disk and hits the second disk), then the SimParticle hitting the second disk is considered
// to be an ancestor SimParticle
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CaloMC/inc/ShowerStepUtil.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>
#include <utility>
#include <numeric>



namespace {

   class CaloCompressUtil
   {
       public:
           CaloCompressUtil() : steps_(), sims_() {}

           const std::vector<const mu2e::StepPointMC*>& steps() const {return steps_;}
           const std::set<art::Ptr<mu2e::SimParticle>>& sims()  const {return sims_;}

           void fill(const mu2e::StepPointMC* step, std::vector<art::Ptr<mu2e::SimParticle>> sims)
           {
               steps_.push_back(step);
               for (const auto& sim: sims) sims_.insert(sim);
           }

        private:
           std::vector<const mu2e::StepPointMC*> steps_;
           std::set<art::Ptr<mu2e::SimParticle>> sims_;
   };

   struct diagSummary
   {
       diagSummary() : totalEdep_(0.0),totalStep_(0),totalSim_(0),totalChk_(0),nCompress_(0),ncompressInfo_(0) {};
       void reset() {totalEdep_=0.0;totalStep_=totalSim_=totalChk_=nCompress_=ncompressInfo_=0;}

       float     totalEdep_;
       unsigned  totalStep_,totalSim_,totalChk_,nCompress_,ncompressInfo_;
   };

}



namespace mu2e {

  class CaloShowerStepMaker : public art::EDProducer
  {
      public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Sequence<std::string>  caloStepPointCollection { Name("caloStepPointCollection"), Comment("Calo crystal stepPointMC collection name") };
             fhicl::Atom<unsigned>         numZSlices              { Name("numZSlices"),              Comment("Number of crystal longitudinal slices ") };
             fhicl::Atom<float>            deltaTime               { Name("deltaTime"),               Comment("Max time difference to be inside a ShowerStep") };
             fhicl::Atom<bool>             compressData            { Name("compressData"),            Comment("Compress stepPointMC and SimParticles in crystal") };
             fhicl::Atom<double>           eDepThreshold           { Name("eDepThreshold"),           Comment("Threshold on energy deposited by SimParticle to keep it") };
             fhicl::Atom<int>              diagLevel               { Name("diagLevel"),               Comment("Debug"),0 };
         };

         explicit CaloShowerStepMaker(const art::EDProducer::Table<Config>& config);
         void produce( art::Event& e) override;


      private:
         using HandleVector = std::vector<art::Handle<StepPointMCCollection>>;
         using SimPtr       = art::Ptr<SimParticle>;
         using SimStepMap   = std::map<SimPtr,std::vector<const StepPointMC*>>;

         void makeCompressedHits       (const HandleVector&, CaloShowerStepCollection&, SimParticlePtrCollection&);
         void collectStepBySimAncestor (const Calorimeter&, const HandleVector&, std::map<SimPtr,CaloCompressUtil>&);
         void collectStepBySim         (const HandleVector&, SimStepMap&);
         void compressSteps            (const Calorimeter&, CaloShowerStepCollection&, int, const SimPtr&, std::vector<const StepPointMC*>&);
         void dumpAllInfo              (const HandleVector&, const Calorimeter&);


         std::vector<std::string>                 calorimeterStepPoints_;
         int                                      numZSlices_;
         double                                   deltaTime_;
         bool                                     compressData_;
         double                                   eDepThreshold_;
         int                                      diagLevel_;
         double                                   zSliceSize_;

         diagSummary                              diagSummary_;
  };



  CaloShowerStepMaker::CaloShowerStepMaker(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},
     calorimeterStepPoints_(config().caloStepPointCollection()),
     numZSlices_           (config().numZSlices()),
     deltaTime_            (config().deltaTime()),
     compressData_         (config().compressData()),
     eDepThreshold_        (config().eDepThreshold()),
     diagLevel_            (config().diagLevel()),
     zSliceSize_(0),
     diagSummary_()
     {
         consumesMany<StepPointMCCollection>();
         produces<CaloShowerStepCollection>();
         produces<SimParticlePtrCollection>();
     }



  //------------------------------------------------------------------------------------------------------------
  //
  void CaloShowerStepMaker::produce(art::Event& event)
  {
      diagSummary_.reset();
      if (diagLevel_ > 0) std::cout << "[CaloShowerStepMaker::produce] begin" << std::endl;

      auto caloShowerStepMCs = std::make_unique<CaloShowerStepCollection>();
      auto simsToKeep        = std::make_unique<SimParticlePtrCollection>();

      HandleVector crystalStepsHandles;
      for (const auto& stepPts : calorimeterStepPoints_)
      {
          art::Handle<StepPointMCCollection> hc;
          event.getByLabel(art::InputTag(stepPts), hc);
          crystalStepsHandles.push_back(hc);
      }

      makeCompressedHits(crystalStepsHandles,*caloShowerStepMCs,*simsToKeep);

      event.put(std::move(caloShowerStepMCs));
      event.put(std::move(simsToKeep));

      if (diagLevel_ > 0) std::cout << "[CaloShowerStepMaker::produce] end" << std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::makeCompressedHits(const HandleVector& crystalStepsHandle,
                                               CaloShowerStepCollection& caloShowerStepMCs,SimParticlePtrCollection& simsToKeep)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      zSliceSize_            = cal.G4Info().getDouble("crystalZLength")/float(numZSlices_)+1e-5;


      //-----------------------------------------------------------------
      // Collect the StepPointMC's produced by each SimParticle Ancestor
      std::map<SimPtr,CaloCompressUtil> crystalAncestorsMap;
      collectStepBySimAncestor(cal,crystalStepsHandle,crystalAncestorsMap);

      if (diagLevel_ > 2) dumpAllInfo(crystalStepsHandle,cal);


      //---------------------------------------------------------------------------------------------------------------
      //Loop over ancestor simParticles, check if they are compressible, and produce the corresponding caloShowerStepMC

      std::set<SimPtr> SimsToKeepUnique;
      for (const auto& iter : crystalAncestorsMap )
      {
          const SimPtr&           sim  = iter.first;
          const CaloCompressUtil& info = iter.second;

          diagSummary_.totalSim_ += info.sims().size();

          std::map<unsigned,std::vector<const StepPointMC*>> crystalMap;
          for (const StepPointMC* step : info.steps()) crystalMap[step->volumeId()].push_back(step);

          for (const auto& iterCrystal : crystalMap)
          {
              unsigned crid = iterCrystal.first;
              std::vector<const StepPointMC*> steps = iterCrystal.second;

              //Filter very small energy deposits at this stage
              double eDep(0);
              for (const auto& step : steps) eDep += step->totalEDep();
              if (eDep < eDepThreshold_) continue;

              if (compressData_)
              {
                  SimsToKeepUnique.insert(sim);
                  compressSteps(cal, caloShowerStepMCs, crid, sim, steps);
              }
              else
              {
                  std::map<SimPtr, std::vector<const StepPointMC*>> newSimStepMap;
                  for (const StepPointMC* step : steps) newSimStepMap[step->simParticle()].push_back(step);
                  for (auto& iter : newSimStepMap)
                  {
                      compressSteps(cal, caloShowerStepMCs, crid, iter.first, iter.second);
                      SimsToKeepUnique.insert(iter.first);
                  }
              }
          }
          ++diagSummary_.ncompressInfo_;
          if (compressData_) ++diagSummary_.nCompress_;
      }

      //dump the unique set of SimParticles to keep into final vector
      simsToKeep.assign(SimsToKeepUnique.begin(),SimsToKeepUnique.end());

      //---------------------------------------------------------------------------------------------------------------
      // Final diag info
      if (diagLevel_ > 1)
      {
          std::cout<<"CaloShowerStepMaker summary"<<std::endl;

          std::set<int> volIds{};
          for (auto caloShowerStepMC : caloShowerStepMCs) volIds.insert(caloShowerStepMC.volumeG4ID());

          for (auto volId: volIds)
          {
             std::map<const art::Ptr<SimParticle>, double> simMap;
             for (const auto& caloShowerStepMC : caloShowerStepMCs)
                if (caloShowerStepMC.volumeG4ID()==volId) simMap[caloShowerStepMC.simParticle()] += caloShowerStepMC.energyDepG4();

             for (auto& kv : simMap) std::cout<<"Vol id: "<<volId<<"  Sim id: "<<kv.first.id()<<"   energy="<<kv.second<<std::endl;
          }
      }

      if (diagLevel_ > 0)
        std::cout << "[CaloShowerStepMaker::makeCompressedHits] compressed "<<diagSummary_.nCompress_<<" / "<<diagSummary_.ncompressInfo_<<" incoming SimParticles"<<std::endl
                  << "[CaloShowerStepMaker::makeCompressedHits] keeping "<<simsToKeep.size()<<" SimParticles"<<std::endl
                  << "[CaloShowerStepMaker::makeCompressedHits] Total sims init: " <<diagSummary_.totalSim_<<std::endl
                  << "[CaloShowerStepMaker::makeCompressedHits] Total caloShower steps: " <<caloShowerStepMCs.size()<<std::endl
                  << "[CaloShowerStepMaker::makeCompressedHits] Total energy deposited / number of stepPointMC: " <<diagSummary_.totalEdep_<<" / "<<diagSummary_.totalStep_<<std::endl
                  << "[CaloShowerStepMaker::makeCompressedHits] Total stepPointMCs seen: " <<diagSummary_.totalChk_<<std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::collectStepBySimAncestor(const Calorimeter& cal,
                                                     const HandleVector& stepsHandles,
                                                     std::map<SimPtr,CaloCompressUtil>& ancestorsMap)
  {
     SimParticlePtrCollection inspectedSims;
     std::unordered_map<SimPtr,SimPtr> simToAncestorMap;

     for (HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
     {
         const art::Handle<StepPointMCCollection>& handle(*i);
         const StepPointMCCollection& steps(*handle);

         for (const auto& step : steps )
         {
             SimPtr sim = step.simParticle();

             inspectedSims.clear();
             while (sim->hasParent())
             {
                 const auto& alreadyInspected = simToAncestorMap.find(sim);
                 if (alreadyInspected != simToAncestorMap.end()) {sim = alreadyInspected->second; break;}
                 inspectedSims.push_back(sim);

                 if (!cal.isInsideAnyCrystal(sim->startPosition()))  break;
                 if (!cal.isInsideSameDisk(sim->startPosition(),sim->endPosition()) ) break;

                 sim = sim->parent();
             }

             for (const SimPtr& inspectedSim : inspectedSims) simToAncestorMap[inspectedSim] = sim;
             ancestorsMap[sim].fill(&step,inspectedSims);

             diagSummary_.totalEdep_ += step.totalEDep();
         }
         diagSummary_.totalStep_ += steps.size();
      }
  }


  //-----------------------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::collectStepBySim(const HandleVector& stepsHandles,
                                             std::map<SimPtr,std::vector<const StepPointMC*>>& simStepMap)
  {
      for (HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i)
      {
          const art::Handle<StepPointMCCollection>& handle(*i);
          const StepPointMCCollection& steps(*handle);
          for (const auto& step : steps ) simStepMap[step.simParticle()].push_back(&step);
      }
  }


  //-------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::compressSteps(const Calorimeter& cal, CaloShowerStepCollection &caloShowerStepMCs,
                                          int volId, const SimPtr& sim, std::vector<const StepPointMC*>& steps)
  {
     auto sortFunctor = [](const StepPointMC* a, const StepPointMC* b) {return a->time() < b->time();};
     std::sort(steps.begin(), steps.end(), sortFunctor);

     ShowerStepUtil buffer(numZSlices_, ShowerStepUtil::weight_type::energy );

     for (const StepPointMC* step : steps)
     {
         CLHEP::Hep3Vector pos  = cal.caloUtil().mu2eToCrystal(volId,step->position());
         int               idx  = int(std::max(1e-6,pos.z())/zSliceSize_);

         if (buffer.entries(idx)>0 && (step->time()-buffer.t0(idx) > deltaTime_) )
         {
             if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepMaker::compressSteps] inserted  "; buffer.printBucket(idx);}
             diagSummary_.totalChk_ += buffer.entries(idx);

             caloShowerStepMCs.push_back(CaloShowerStep(volId, sim, buffer.entries(idx), buffer.time(idx), buffer.energyG4(idx),
                                                        buffer.energyVis(idx),buffer.pIn(idx),buffer.pos(idx)));
             buffer.reset(idx);
         }

         buffer.add(idx, step->totalEDep(), step->visibleEDep(), step->time(), step->momentum().mag(), pos);
     }

     //do not forget to flush the final buffer(s) :-)
     for (unsigned i=0;i<buffer.nBuckets();++i)
     {
         if (buffer.entries(i) == 0) continue;

         if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepMaker::compressSteps] inserted ";  buffer.printBucket(i);}
         diagSummary_.totalChk_ += buffer.entries(i);

         caloShowerStepMCs.push_back(CaloShowerStep(volId, sim,  buffer.entries(i), buffer.time(i), buffer.energyG4(i),
                                                    buffer.energyVis(i),buffer.pIn(i),buffer.pos(i)));
     }
  }

  //-------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::dumpAllInfo(const HandleVector& stepsHandles, const Calorimeter& cal)
  {
      std::cout<<"Dumping StepPointMCs  Mu2e / crystal / disk / diskFF frames"<<std::endl;
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
      {
          const art::Handle<StepPointMCCollection>& handle(*i);
          const StepPointMCCollection& steps(*handle);

          std::cout<<steps.size()<<std::endl;
          for (const auto& step : steps )
            std::cout<<step.volumeId()<<" "<<step.totalEDep()<<" "<<step.position()<<" "
                     <<cal.caloUtil().mu2eToCrystal(step.volumeId(),step.position())<<"   "
                     <<cal.caloUtil().mu2eToDisk(cal.crystal(step.volumeId()).diskID(),step.position())<<"   "
                     <<cal.caloUtil().mu2eToDiskFF(cal.crystal(step.volumeId()).diskID(),step.position())<<std::endl;
      }
  }

}

using mu2e::CaloShowerStepMaker;
DEFINE_ART_MODULE(CaloShowerStepMaker)
