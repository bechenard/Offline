//
// Geometry and identifier info about the Calorimeter.
//
// Original author B. Echenard
//
#include "Offline/CalorimeterGeom/inc/CaloUtil.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  CaloUtil::CaloUtil(const DiskCalorimeter& owner) :
     owner_(owner),
     trackerCenter_()
  {}

  //-----------------------------------------------------------------------------
  void CaloUtil::trackerCenter(const CLHEP::Hep3Vector& vec) {trackerCenter_ = vec;}
  const CLHEP::Hep3Vector& CaloUtil::trackerCenter() const {return trackerCenter_;}

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::mu2eToCrystal(int crystalId, const CLHEP::Hep3Vector& pos) const
  {
    const auto& crystal = owner_.crystal(crystalId);
    const Disk& thisDisk = owner_.disk(crystal.diskID());
    CLHEP::Hep3Vector crysLocalPos = crystal.localPosition();
    return thisDisk.diskInfo().rotation()*(pos-thisDisk.diskInfo().origin())-crysLocalPos;
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::mu2eToDisk(int diskId, const CLHEP::Hep3Vector& pos) const
  {
    const Disk& thisDisk = owner_.disk(diskId);
    return (thisDisk.diskInfo().rotation())*(pos-thisDisk.diskInfo().origin());
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::mu2eToDiskFF(int diskId, const CLHEP::Hep3Vector& pos) const
  {
    const Disk& thisDisk = owner_.disk(diskId);
    return (thisDisk.diskInfo().rotation())*(pos-thisDisk.diskInfo().origin()) -
            thisDisk.diskInfo().originToCrystalOrigin();
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::mu2eToTracker(const CLHEP::Hep3Vector& pos) const
  {
    return pos - trackerCenter_;
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::crystalToMu2e(int crystalId, const CLHEP::Hep3Vector& pos) const
  {
    const auto& crystal = owner_.crystal(crystalId);
    const Disk& thisDisk = owner_.disk(crystal.diskID());
    CLHEP::Hep3Vector crysLocalPos = crystal.localPosition();
    return thisDisk.diskInfo().inverseRotation()*(pos+crysLocalPos) + thisDisk.diskInfo().origin();
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::diskToMu2e(int diskId, const CLHEP::Hep3Vector& pos) const
  {
    const Disk& thisDisk = owner_.disk(diskId);
    return thisDisk.diskInfo().inverseRotation()*pos + thisDisk.diskInfo().origin();
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::diskFFToMu2e(int diskId, const CLHEP::Hep3Vector& pos) const
  {
    const Disk& thisDisk = owner_.disk(diskId);
    return thisDisk.diskInfo().inverseRotation()*(pos + thisDisk.diskInfo().originToCrystalOrigin())
           + thisDisk.diskInfo().origin();
  }

  //-----------------------------------------------------------------------------
  CLHEP::Hep3Vector CaloUtil::trackerToMu2e(const CLHEP::Hep3Vector& pos) const
  {
    return pos + trackerCenter_;
  }

}
