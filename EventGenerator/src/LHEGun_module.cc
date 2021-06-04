#include <iostream>
#include <string>
#include <algorithm>

#include "cetlib_except/exception.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"  
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/LHEReader.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "TH1.h"


namespace mu2e {

  class LHEGun : public art::EDProducer {

     public:
        using RTS = mu2e::RootTreeSampler<IO::StoppedParticleF>;
        struct Config 
        {
            using Name    = fhicl::Name;
            using Comment = fhicl::Comment;            
            fhicl::Table<RTS::Config>  muonStops      { Name("muonStops"),     Comment("Stopped muon config") };
            fhicl::Atom<std::string>   LHEfilename    { Name("LHEfilename"),   Comment("LHE data file name") }; 
            fhicl::Atom<std::string>   genID          { Name("genID"),         Comment("Generation ID name") }; 
            fhicl::Atom<bool>          doHistograms   { Name("doHistograms"),  Comment("Do histograms"),false }; 
            fhicl::Atom<int>           diagLevel      { Name("diagLevel"),     Comment("Diag Level"),0 };
        };

        explicit LHEGun(const art::EDProducer::Table<Config>& config);
        virtual void produce(art::Event& event);


     private:
        void fill(GenParticleCollection& output );

        LHEReader                 reader_;
        CLHEP::HepRandomEngine&   engine_;
        RTS                       stops_;
        GenId                     genId_;
        bool                      doHistograms_;
        int                       diagLevel_;
        double                    mumass2_;
        double                    elmass2_;

        TH1F*                     hEnergy_;
        TH1F*                     hPdgId_;
        TH1F*                     hTime_;
        TH1F*                     hZ_;
  };




  //================================================================
  LHEGun::LHEGun(const art::EDProducer::Table<Config>& config) :
     EDProducer{config},
     reader_            (config().LHEfilename()),            
     engine_            (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
     stops_             (engine_, config().muonStops()),
     genId_             (GenId::findByName(config().genID())),
     doHistograms_      (config().doHistograms()),
     diagLevel_         (config().diagLevel())
  {
      produces<GenParticleCollection>();
      mumass2_ = std::pow(GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::type(13)).ref().mass().value(),2);
      elmass2_ = std::pow(GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::type(11)).ref().mass().value(),2);

      if (doHistograms_)
      {
	 art::ServiceHandle<art::TFileService> tfs;
	 hEnergy_ = tfs->make<TH1F>("hEnergy", "Energy"      , 2400,   0.0,  120);
	 hPdgId_  = tfs->make<TH1F>("hPdgId" , "PDG ID"      ,  500,  -250, 250);
	 hTime_   = tfs->make<TH1F>("hTime"  , "Time"        ,  400,   0.0, 2000.);
	 hZ_      = tfs->make<TH1F>("hZ"     , "Z"           ,  500,  5400, 6400);
      }
  }








  //================================================================
  void LHEGun::produce(art::Event& event)
  {
     std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
     fill(*output);
     event.put(std::move(output));
   }




  //================================================================
  void LHEGun::fill(GenParticleCollection& output )
  {

     static constexpr double GeVToMeV(1000.0);

     const auto& stop = stops_.fire();
     const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
     if (doHistograms_) hTime_->Fill(stop.t);

     reader_.readNextEvent();

     //must abort if no more events
     if (!reader_.isValid()) return;

     for (const auto& particle : reader_.event().particles())
     {
         if (particle.istup() !=1 ) continue;

         int pdgId = particle.idup();
         CLHEP::Hep3Vector p3(particle.pup(0)*GeVToMeV,particle.pup(1)*GeVToMeV,particle.pup(2)*GeVToMeV);

         double energy(particle.pup(4));
         if (std::abs(pdgId) == 11) energy = sqrt(p3.mag2()+elmass2_);  //Madgraph sets m_e = m_mu = 0
         if (std::abs(pdgId) == 13) energy = sqrt(p3.mag2()+mumass2_);
         CLHEP::HepLorentzVector fourmom(p3,energy);
         if (diagLevel_>1) std::cout<<"LHEGun::fill adding paarticle p=("<<p3.x()<<","<<p3.y()<<","<<p3.z()
                                    <<") E="<<energy<<std::endl;

	 output.emplace_back(PDGCode::type(pdgId),genId_,pos,fourmom,stop.t);

	 if (doHistograms_)
	 {
	    hPdgId_->Fill(pdgId);
	    hEnergy_->Fill(energy);
	    hZ_->Fill(pos.z());
	 }
     }
   }

}

DEFINE_ART_MODULE(mu2e::LHEGun);
