//
// Create a disk and fills it with crystals.
//
// Original author B Echenard
//
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/SquareMapper.hh"
#include "Offline/CalorimeterGeom/inc/SquareShiftMapper.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <algorithm>
#include <map>
#include <fstream>
#include <ranges>

namespace mu2e {

  Disk::Disk(unsigned id, double rCrystalIn, double rCrystalOut,
             double nominalCellSize, double nominalCellLength, unsigned offset,
             const CLHEP::Hep3Vector& crystalOrigInDisk,
             const std::string& crystalFileName) :
    id_(id),
    crystals_(),
    diskInfo_(),
    radiusInCrystal_(rCrystalIn),
    radiusOutCrystal_(rCrystalOut),
    nominalCellSize_(nominalCellSize),
    crystalOffset_(offset),
    mapToCrystal_(),
    crystalToMap_()
  {
     crystalMap_ = std::shared_ptr<CrystalMapper>(new SquareShiftMapper());
     fillCrystals(crystalOrigInDisk, nominalCellLength, crystalFileName );
  }


  //-----------------------------------------------------------------------------
  // take the crystals from their measured location and correlate them with ideal map location
  void Disk::fillCrystals(const CLHEP::Hep3Vector& crystalOrigInDisk,
                          double nominalCellLength,
                          const std::string& fileName)
  {
    //Start by reading the actual location of the crystals from the file
    ConfigFileLookupPolicy configFile;
    std::string fullFileName = configFile(fileName);

    std::ifstream crysFile;
    crysFile.open(fullFileName);
    if (!crysFile.is_open()) {
      throw cet::exception("DISK_OPEN_FAILED")<<"failed to open file " << fullFileName << "\n";
    }

    std::vector<float> realPosX,realPosY;
    float  x,y;
    unsigned did,cidx;
    while (crysFile >> did >> cidx >> x >> y){
      if (did != id_) continue;
      if (realPosX.size() != cidx) throw cet::exception("DISK_READ_FAILED")<<" invalid crystal index " << cidx << "\n";
      realPosX.push_back(x);
      realPosY.push_back(y);
    }
    crysFile.close();


    // Match the crystals to the ideal location map to facilitate navigation
    unsigned nCrystal(0);
    unsigned nRingsMax   = unsigned(1.5*radiusOutCrystal_/nominalCellSize_);
    unsigned nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);
    mapToCrystal_.assign(nCrystalMap, Crystal::invalidID_);

    for (unsigned mapIdx=0;mapIdx<nCrystalMap;++mapIdx)
    {
      CLHEP::Hep2Vector posIdeal  = nominalCellSize_*crystalMap_->xyFromIndex(mapIdx);

      size_t idx=0;
      for (;idx<realPosX.size();++idx){
        if (std::abs(realPosX[idx]-posIdeal.x())/nominalCellSize_ < 0.5 &&
            std::abs(realPosY[idx]-posIdeal.y())/nominalCellSize_ < 0.5) break;
      }
      if (idx == realPosX.size() ) continue;

      checkPosition(realPosX[idx]+crystalOrigInDisk.x(),realPosY[idx]+crystalOrigInDisk.y(),
                    nominalCellSize_, nominalCellSize_, idx);

      CLHEP::Hep3Vector size(nominalCellSize_,nominalCellSize_,nominalCellLength);
      CLHEP::Hep3Vector posFF(realPosX[idx],realPosY[idx],0);
      CLHEP::Hep3Vector pos = posFF + crystalOrigInDisk;

      mapToCrystal_[mapIdx] = nCrystal;
      crystalToMap_.push_back(mapIdx);
      crystals_.push_back(Crystal(nCrystal, id_, pos, size));

      ++nCrystal;
    }
  }


  //-----------------------------------------------------------------------------
  void Disk::checkPosition(float xPos, float yPos, float xsize, float ysize, unsigned crystalID) const
  {
    for (const auto& other : crystals_) {
      float x0 = other.localPosition().x();
      float y0 = other.localPosition().y();
      float dx = std::abs(x0-xPos);
      float dy = std::abs(y0-yPos);
      float sx = 0.4999*(other.size().x() + xsize);
      float sy = 0.4999*(other.size().y() + ysize);

      if (dx < sx && dy < sy) throw cet::exception("DISK_POS_FAILED")<<
          "crystals" << crystalID << " and " << other.localID()<<" overlap\n";
    }
  }


  //-----------------------------------------------------------------------------
  std::vector<unsigned> Disk::idxFromRow(int thisRow) const
  {
    std::vector<unsigned> cryList;
    for (size_t i=0;i<crystals_.size();++i){
      int irow = crystalMap_->rowFromIndex(crystalToMap_[i]);
      if (irow == thisRow) cryList.push_back(i);
    }
    return cryList;
  }


  //-----------------------------------------------------------------------------
  // Position needs to be in the Disk frame
  bool Disk::isInsideDisk(const CLHEP::Hep3Vector& posInDisk) const
  {
    if (posInDisk.z() < -diskInfo_.size().z()/2.0 || posInDisk.z() > diskInfo_.size().z()/2.0) return false;
    float radius = posInDisk.perp();
    if (radius < diskInfo_.size().x() || radius > diskInfo_.size().y()) return false;
    return true;
  }



  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  bool Disk::isInsideCrystal(const CLHEP::Hep3Vector& pos) const
  {
    const double tolerance(1e-6);
    if (pos.z() < -tolerance) return false;
    if (pos.z() > diskInfo_.crystalZLength()+tolerance) return false;
    if (idxFromPosition(pos) == Crystal::invalidID_) return false;
    return true;
  }

  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  unsigned Disk::idxFromPosition(const CLHEP::Hep3Vector& pos) const
  {
    // First, filter out obvious misses
    const double tolerance(1e-6);
    float perp = pos.perp();
    if (perp < diskInfo_.size().x() ||
        perp > diskInfo_.size().y() ||
        pos.z() < -tolerance        ||
        pos.z() > diskInfo_.crystalZLength()+tolerance) return Crystal::invalidID_;

    // now look at crystal closest to ideal position - should work most of the time
    auto mapIdx = crystalMap_->indexFromXY(pos.x()/nominalCellSize_,pos.y()/nominalCellSize_);
    if (isCrystalIdxValid(mapIdx)) {
      if (isInsideCrystal(mapToCrystal_[mapIdx],pos)) return mapToCrystal_[mapIdx];

      // if no match, look at neighbors (likely there)
      const int level(1);
      const auto neighbors(crystalMap_->neighbors(mapIdx,level));
      for (const auto& idx : neighbors) {
         if (isCrystalIdxValid(idx) && isInsideCrystal(mapToCrystal_[idx],pos)) return mapToCrystal_[idx];
      }
    }

    // last chance, look at all crystals
    for (size_t icry=0;icry<crystals_.size();++icry){
      if (isInsideCrystal(icry,pos)) return icry;
    }

    return Crystal::invalidID_;
  }

  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  bool Disk::isInsideCrystal(unsigned icry, const CLHEP::Hep3Vector& pos) const
  {
    if (icry >= crystals_.size()) return false;

    const auto& cpos  = crystals_[icry].localPosition();
    const auto& csize = crystals_[icry].size();
    if (std::abs(cpos.x()-pos.x()) > 0.5*csize.x() ) return false;
    if (std::abs(cpos.y()-pos.y()) > 0.5*csize.y() ) return false;
    if (pos.z() < -1e-6 || pos.z() > csize.z()     ) return false;

    return true;
  }

  //-----------------------------------------------------------------------------
  //find the indices of the crystal neighbors for a given level (level = number of rings away)
  // local in, local out
  std::vector<unsigned> Disk::neighbors(unsigned localId, unsigned level) const
  {
    std::vector<unsigned> list;
    const auto mapNbrs = crystalMap_->neighbors(crystalToMap_.at(localId), level);
    for (const auto& mapIdx : mapNbrs)
      if (isCrystalIdxValid(mapIdx)) list.push_back(mapToCrystal_.at(mapIdx));
    return list;
  }

  //-----------------------------------------------------------------------------
  bool Disk::isCrystalIdxValid(unsigned i) const {
    return i < mapToCrystal_.size() && mapToCrystal_[i] != Crystal::invalidID_;
  }

  //-----------------------------------------------------------------------------
  //Move and rotate the disk. The displacement and rotation are relatives
  void Disk::moveDisk(const CLHEP::Hep3Vector& shift, const CLHEP::HepRotation& rot)
  {
    diskInfo_.setPose(diskInfo_.origin() + shift, diskInfo_.rotation() * rot);
    for (auto& c : crystals_) c.setPosition(diskInfo_.toGlobal(c.localPosition()));
  }

  //-----------------------------------------------------------------------------
  void Disk::moveCrystal(unsigned localId, const CLHEP::Hep3Vector& shift)
  {
    auto& crystal = crystals_.at(localId);
    crystal.setLocalPosition(crystal.localPosition() + shift);
    crystal.setPosition(diskInfo_.toGlobal(crystal.localPosition()));
  }

  //-----------------------------------------------------------------------------
  void Disk::print(std::ostream &os) const
  {
     os<<"Disk                  "<<id_<<std::endl;
     os<<"Number of crystals    "<<crystals_.size()<<std::endl;
     os<<"Crystal offset        "<<crystalOffset_<<std::endl;
     diskInfo_.print(os);
  }

}

