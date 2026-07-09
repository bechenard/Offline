#include "cetlib_except/exception.h"
#include "Offline/CaloMC/inc/CaloNoiseSampler.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TFile.h"
#include "TH1F.h"
#include "TKey.h"

#include <memory>
#include <vector>
#include <iostream>


namespace mu2e {


   CaloNoiseSampler::CaloNoiseSampler(const std::string& fileName) :
      fileName_   {fileName},
      prefix_     {"histo_"},
      histoBaseID_{0},
      noiseMap_   {}
   {}


   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSampler::fillCache(int histoBaseID)
   {
      // Cache is for all baseID, clear it
      noiseMap_.clear();
      histoBaseID_ = histoBaseID;

      // Refill the cache with all histos sharing the same baseID
      ConfigFileLookupPolicy resolveFullPath;
      std::string fullFileName = resolveFullPath(fileName_);

      TFile file(fullFileName.c_str());
      if (!file.IsOpen()) throw cet::exception("NOISEREADER")<<"Filename"<<fullFileName.c_str()<<" does not exist\n";

      TIter next(file.GetListOfKeys());
      TKey* key;

      while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj->InheritsFrom(TH1F::Class()))
          continue;

        TH1F* histo = static_cast<TH1F*>(obj);
        std::string name = histo->GetName();

        // Parse the integer after the prefix
        std::string suffix = name.substr(prefix_.size());
        std::istringstream iss(suffix);

        int hid;
        if (!(iss >> hid) || !iss.eof())
          throw cet::exception("NOISEREADER")<<"Hitsogram "<<name.c_str()<<" is invalid\n";

        // Only take histos in the right batch to reduce memory footprint
        // the batch is defined as histogrm id/base

        if (hid/base != histoBaseID/base) continue;

        const float* array = histo->GetArray();
        noiseMap_[hid].assign(array + 1, array + histo->GetNbinsX() + 1);
      }
   }


   //----------------------------------------------------------------------------------------------------------------------
   const std::vector<float>& CaloNoiseSampler::getNoiseWF(int histoID)
   {
     int baseID = histoID/base;

     if (baseID != histoBaseID_) fillCache(baseID);

     auto iter = noiseMap_.find(histoID);
     if (iter == noiseMap_.end())
        throw cet::exception("NOISEREADER")<<"histoID "<<histoID<<" is invalid\n";
     return iter->second;
   }


   //----------------------------------------------------------------------------------------------------------------------
   void CaloNoiseSampler::printCache()
   {
     std::cout<<"CaloNoiseSampler cache\n";
     for (const auto& kv : noiseMap_) std::cout<<"Histo id "<<kv.first<<"  noise waveform length "<<kv.second.size()<<"\n";
   }
}
