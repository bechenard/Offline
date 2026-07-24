//
// Contains data for the calorimeter class. The non-critical fields are saved into maps, and a few performance critical fields
// accessed throughout the code are cached for efficiency.
//
// See wiki for list of valid keys
//
// Original author B. Echenard
//
#ifndef CalorimeterGeom_CaloG4Info_hh
#define CalorimeterGeom_CaloG4Info_hh

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include <map>
#include <string>


namespace mu2e {

  template <typename T> class CaloG4InfoData
  {
    public:
       CaloG4InfoData() : data_() {};

       void set(const std::string& key, const T& value) {data_[key] = value;}

       const T& get(const std::string& key) const
       {
          auto iter = data_.find(key);
          if (iter == data_.end()) throw cet::exception("CaloG4Info") << " unknown element "<<key<<"\n";
          return iter->second;
       };

    private:
       std::map<std::string,T> data_;
  };



  class CaloG4Info {
    public:
      CaloG4Info() : dataBool_(), dataInt_(), dataDouble_(), dataVInt_(),
                     dataVDouble_(), dataString_(),dataH3V_()
      {}

      void set(const std::string& key, bool value)                       {dataBool_.set(key,value);}
      void set(const std::string& key, int value)                        {dataInt_.set(key,value);}
      void set(const std::string& key, double value)                     {dataDouble_.set(key,value);}
      void set(const std::string& key, const std::vector<int>& value)    {dataVInt_.set(key,value);}
      void set(const std::string& key, const std::vector<double>& value) {dataVDouble_.set(key,value);}
      void set(const std::string& key, const std::string& value)         {dataString_.set(key,value);}
      void set(const std::string& key, const CLHEP::Hep3Vector& value)   {dataH3V_.set(key,value);}

      bool                       getBool   (const std::string& key) const {return dataBool_.get(key);}
      int                        getInt    (const std::string& key) const {return dataInt_.get(key);}
      double                     getDouble (const std::string& key) const {return dataDouble_.get(key);}
      const std::vector<int>&    getVInt   (const std::string& key) const {return dataVInt_.get(key);}
      const std::vector<double>& getVDouble(const std::string& key) const {return dataVDouble_.get(key);}
      const std::string&         getString (const std::string& key) const {return dataString_.get(key);}
      const CLHEP::Hep3Vector&   getHepVec (const std::string& key) const {return dataH3V_.get(key);}



    private:
      CaloG4InfoData<bool>                dataBool_;
      CaloG4InfoData<int>                 dataInt_;
      CaloG4InfoData<double>              dataDouble_;
      CaloG4InfoData<std::vector<int>>    dataVInt_;
      CaloG4InfoData<std::vector<double>> dataVDouble_;
      CaloG4InfoData<std::string>         dataString_;
      CaloG4InfoData<CLHEP::Hep3Vector>   dataH3V_;
  };
}

#endif
