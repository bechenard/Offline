#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/CaloUtil.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <iostream>
#include <algorithm>


namespace mu2e {

  DiskCalorimeter::DiskCalorimeter() :
    ProditionsEntity(Calorimeter::cxname),
    disks_(),
    crystals_(),
    G4Info_(),
    util_(*this)
  {
    rebuildCrystalPtrs();
  }

  DiskCalorimeter::DiskCalorimeter(const DiskCalorimeter& rhs) :
    ProditionsEntity(rhs),
    disks_(rhs.disks_),
    G4Info_(rhs.G4Info_),
    util_(*this)
  {
    rebuildCrystalPtrs();
    util_.trackerCenter(rhs.util_.trackerCenter());
  }

  DiskCalorimeter::DiskCalorimeter(DiskCalorimeter&& rhs) noexcept :
    ProditionsEntity(std::move(rhs)),
    disks_(std::move(rhs.disks_)),
    G4Info_(std::move(rhs.G4Info_)),
    util_(*this)
  {
    rebuildCrystalPtrs();
    util_.trackerCenter(std::move(rhs.util_.trackerCenter()));
  }

  void DiskCalorimeter::rebuildCrystalPtrs()
  {
    crystals_.clear();
    for (const auto& disk : disks_) {
      for (size_t i = 0; i < disk.nCrystals(); ++i) {
        crystals_.push_back(&disk.crystal(i));
      }
    }
  }


  bool DiskCalorimeter::isInsideAnyDisk(const CLHEP::Hep3Vector& pos) const
  {
    for (const auto& disk : disks_){
      CLHEP::Hep3Vector posInDisk = util_.mu2eToDisk(disk.id(),pos);
      if (disk.isInsideDisk(posInDisk)) return true;
    }
    return false;
  }


  bool DiskCalorimeter::isInsideAnyCrystal(const CLHEP::Hep3Vector& pos) const
  {
    for (const auto& disk : disks_){
      CLHEP::Hep3Vector posInDisk = util_.mu2eToDiskFF(disk.id(),pos);
      if (disk.isInsideCrystal(posInDisk)) return true;
    }
    return false;
  }


  bool DiskCalorimeter::isInsideSameDisk(const CLHEP::Hep3Vector& front,
                                         const CLHEP::Hep3Vector& back) const
  {
    for (const auto& disk : disks_) {
      CLHEP::Hep3Vector frontInDisk = util_.mu2eToDiskFF(disk.id(),front);
      CLHEP::Hep3Vector backInDisk  = util_.mu2eToDiskFF(disk.id(),back);
      if (disk.isInsideDisk(frontInDisk) && disk.isInsideDisk(backInDisk)) return true;
    }
    return false;
  }

  void DiskCalorimeter::print(std::ostream &os) const
  {
     os<<"Disk calorimeter "<<std::endl;
     os<<"Number of disks :"<< disks_.size()<<std::endl;
     for (size_t idisk=0;idisk<disks_.size();++idisk) disk(idisk).print(os);
  }
}
