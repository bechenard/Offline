#include "Offline/CaloMC/inc/CaloNoiseSimGenerator.hh"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Offline/SeedService/inc/SeedService.hh"

#include "TFile.h"
#include "TH2.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TDirectory.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>




namespace mu2e {

   CaloNoiseSimGenerator::CaloNoiseSimGenerator(const Config& config, CLHEP::HepRandomEngine& engine) :
     waveform_      (config.noiseWFSize(),0.0),
     pedestal_      (0.0),
     digiSampling_  (config.digiSampling()),
     noiseRinDark_  (config.rinNphotPerNs() + config.darkNphotPerNs()),
     noiseElec_     (config.elecNphotPerNs()),
     minPeakADC_    (config.minPeakADC()),
     pePerMeV_      (config.pePerMeV()),
     MeVToADC_      (config.MeVToADC()),
     randPoisson_   (engine),
     randGauss_     (engine),
     randFlat_      (engine),
     pulseShape_    (config.pulseFileName(),config.pulseHistName(),digiSampling_),
     diagLevel_     (config.diagLevel())
   {}


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::initialize()
   {
       pulseShape_.buildShapes();
       generateWF();
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::refresh() {generateWF();}

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::generateWF()
   {
       float scaleFactor(MeVToADC_/pePerMeV_);

       std::fill(waveform_.begin(),waveform_.end(),0);

       const auto&      pulse         = pulseShape_.digitizedPulse(0.0);
       const unsigned   pulseSize     = pulse.size();
       const unsigned   bufferSize    = int(0.75*pulseSize);
       const unsigned   noiseSize     = waveform_.size();
       const double     totalTime     = (noiseSize+bufferSize)*digiSampling_;
       const int        noiseLevelPE  = int(totalTime*noiseRinDark_);

       //Generate the radiation induced noise (RIN)
       const int nPh = randPoisson_(noiseLevelPE);
       for (int i=0;i<nPh;++i)
       {
           double t0 = randFlat_.fire(0.0,totalTime);
           const auto& wf = pulseShape_.digitizedPulse(t0);

           int i0 = int(t0/digiSampling_) - bufferSize;
           int l0 = (i0<0) ? -i0 : 0;
           int l1 = std::min(pulseSize,noiseSize-i0);
           for (int l=l0;l<l1;++l) waveform_[i0+l] += wf[l]*scaleFactor;
       }

       //add electronics noise
       double noiseADC = noiseElec_*digiSampling_*scaleFactor;
       for (auto& val : waveform_) val += randGauss_.fire(0.0,noiseADC);

       //estimate pedestal for this waveform - set it to theoretical value for the time being
       pedestal_ = std::trunc(noiseRinDark_*digiSampling_*std::accumulate(pulse.begin(),pulse.end(),0.0)*scaleFactor);
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::addSampleNoise(std::vector<double>& wfVector, unsigned istart, unsigned ilength)
   {
       if (ilength >=waveform_.size())
          throw cet::exception("CATEGORY")<<"[CaloNoiseSimGenerator] noise length request too long";

       unsigned irandom = unsigned(randFlat_.fire(0.,waveform_.size()-ilength));
       for (unsigned i=0;i<ilength;++i) wfVector[istart+i] += waveform_[irandom+i];
   }


   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::plotNoise(const std::string& name)
   {
       std::vector<double> x(waveform_.size()),y;
       std::iota(x.begin(),x.end(),0);
       TGraph gr(x.size(),x.data(),waveform_.data());
       gr.SetTitle("Original waveform");

       TH1F h1("h1","Projection waveform",100,-50,50);
       for (const auto& val: waveform_) h1.Fill(val-pedestal_);

       TCanvas c1("c1","c1");
       c1.Divide(2,2);
       c1.cd(1);
       gr.Draw("AL");
       c1.cd(2);
       h1.Draw();
       c1.SaveAs(name.c_str());
   }

   //------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSimGenerator::dumpNoise(const std::string& fname)
   {
      TFile outfile(fname.c_str(), "RECREATE");

      TH1F h("histo_0","histo_0", waveform_.size(), 0, waveform_.size());
      for (size_t i = 0; i < waveform_.size(); ++i) h.SetBinContent(i+1, waveform_[i]);

      h.Write();
      outfile.Close();
      std::cout<<"CaloNoiseSimGenerator written waveform in "<<fname<<"\n";
   }

}


