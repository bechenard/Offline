#ifndef CaloNoiseSampler_HH
#define CaloNoiseSampler_HH
//
// Cache and provide noise waveforms for readouts
//
#include <map>
#include <vector>


namespace mu2e {

  class CaloNoiseSampler
  {
        CaloNoiseSampler(const std::string& fileName);

        const std::vector<float>& getNoiseWF(int histoID);
        void                printCache();


     private:
        void fillCache(int histoBaseID);

        std::string                      fileName_;
        std::string                      prefix_;
        int                              histoBaseID_;
        std::map<int,std::vector<float>> noiseMap_;

        static constexpr int base = 10000;
   };

}
#endif
