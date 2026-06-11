#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh
//
// Hold information about a crystal
// ID, Neighbors, and position are given in the "Mu2e" global reference frame
// localID and localPosition are given in the local disk refernce frame
//
// Original author B Echenard
//
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

  class Crystal {
    public:
      static constexpr int invalidID_ = CaloConst::_invalid;

      Crystal(int localID, int diskID, const CLHEP::Hep3Vector& localPosition,
              const CLHEP::Hep3Vector& size) :
         ID_(-1),
         localID_(localID),
         diskID_(diskID),
         size_(size),
         localPosition_(localPosition),
         position_(),
         neighbors_(),
         nextNeighbors_()
      {}

      int                       ID           () const {return ID_;}
      int                       localID      () const {return localID_;}
      int                       diskID       () const {return diskID_;}
      const CLHEP::Hep3Vector&  size         () const {return size_;}
      const CLHEP::Hep3Vector&  localPosition() const {return localPosition_;}
      const CLHEP::Hep3Vector&  position     () const {return position_;}
      const std::vector<int>&   neighbors    () const {return neighbors_;}
      const std::vector<int>&   nextNeighbors() const {return nextNeighbors_;}

      void setID           (int ID)                       {ID_ = ID;}
      void setlocalID      (int localID)                  {localID_ = localID;}
      void setdiskID       (int diskID)                   {diskID_ = diskID;}
      void setLocalPosition(const CLHEP::Hep3Vector& pos) {localPosition_ = pos;}
      void setPosition     (const CLHEP::Hep3Vector& pos) {position_ = pos;}
      void setNeighbors    (const std::vector<int>& list) {neighbors_ = list;}
      void setNextNeighbors(const std::vector<int>& list) {nextNeighbors_ = list;}

    private:
      int               ID_;
      int               localID_;
      int               diskID_;
      CLHEP::Hep3Vector size_;
      CLHEP::Hep3Vector localPosition_;
      CLHEP::Hep3Vector position_;
      std::vector<int>  neighbors_;
      std::vector<int>  nextNeighbors_;
  };
}
#endif
