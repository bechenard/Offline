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

  Disk::Disk(int id, double rCrystalIn, double rCrystalOut, double nominalCellSize,
             double nominalCellLength, int offset, const CLHEP::Hep3Vector& crystalOrigInDisk,
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
    int did,cidx;
    while (crysFile >> did >> cidx >> x >> y){
      if (did != id_) continue;
      if (std::ssize(realPosX) != cidx) throw cet::exception("DISK_READ_FAILED")<<" invalid crystal index " << cidx << "\n";
      realPosX.push_back(x);
      realPosY.push_back(y);
    }
    crysFile.close();


    // Match the crystals to the ideal location map to facilitate navigation
    int nCrystal(0);
    int nRingsMax   = int(1.5*radiusOutCrystal_/nominalCellSize_);
    int nCrystalMap = crystalMap_->nCrystalMax(nRingsMax);
    mapToCrystal_.insert(mapToCrystal_.begin(), nCrystalMap, Crystal::invalidID_);

    for (int mapIdx=0;mapIdx<nCrystalMap;++mapIdx)
    {
      CLHEP::Hep2Vector posIdeal  = nominalCellSize_*crystalMap_->xyFromIndex(mapIdx);

      size_t idx=0;
      for (;idx<realPosX.size();++idx){
        if (std::abs(realPosX[idx]-posIdeal.x())/nominalCellSize_ < 0.5 &&
            std::abs(realPosY[idx]-posIdeal.y())/nominalCellSize_ < 0.5) break;
      }
      if (idx == realPosX.size() ) continue;

      checkPosition(realPosX[idx],realPosY[idx],nominalCellSize_,nominalCellSize_);

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
  void Disk::checkPosition(float xPos, float yPos, float xsize, float ysize) const
  {
    for (const auto& other : crystals_) {
      float x0 = other.localPosition().x();
      float y0 = other.localPosition().y();
      float dx = std::abs(x0-xPos);
      float dy = std::abs(y0-yPos);
      float sx = 0.4999*(other.size().x() + xsize);
      float sy = 0.4999*(other.size().y() + ysize);

      if (dx < sx && dy < sy) throw cet::exception("DISK_POS_FAILED")<<
          "crystals" << crystals_.size()+1 << " and " << other.localID()<<" overlap\n";
    }
  }


  //-----------------------------------------------------------------------------
  std::vector<int> Disk::idxFromRow(int thisRow) const
  {
    std::vector<int> cryList;
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
    if (posInDisk.z() < -diskInfo_.size().z() || posInDisk.z() > diskInfo_.size().z()) return false;
    float radius = posInDisk.perp();
    if (radius < diskInfo_.size().x() || radius > diskInfo_.size().y()) return false;
    return true;
  }



  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  bool Disk::isInsideCrystal(const CLHEP::Hep3Vector& posinDiskFF) const
  {
    const double tolerance(1e-6);
    if (posinDiskFF.z() < -tolerance) return false;
    if (posinDiskFF.z() > diskInfo_.crystalZLength()+tolerance) return false;
    if (idxFromPosition(posinDiskFF) == Crystal::invalidID_) return false;
    return true;
  }

  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  int Disk::idxFromPosition(const CLHEP::Hep3Vector& pos) const
  {
    // First, filter out obvious misses
    float perp = pos.perp();
    if (perp < diskInfo_.size().x() ||
        perp > diskInfo_.size().y() ||
        std::abs(pos.z())>diskInfo_.size().z()) return Crystal::invalidID_;

    // now look at crystal closest to ideal position - should work most of the time
    int mapIdx = crystalMap_->indexFromXY(pos.x()/nominalCellSize_,pos.y()/nominalCellSize_);
    if (isCrystalIdxValid(mapIdx) && isInsideCrystal(mapToCrystal_[mapIdx],pos)) return mapToCrystal_[mapIdx];

    // if no match, look at neighbors (likely there)
    const int level(1);
    const auto neighbors(crystalMap_->neighbors(mapIdx,level));
    for (const auto& idx : neighbors) {
       if (isCrystalIdxValid(idx) && isInsideCrystal(mapToCrystal_[idx],pos)) return mapToCrystal_[idx];
    }

    // last chance, look at all crystals
    for (size_t icry=0;icry<crystals_.size();++icry){
      if (isInsideCrystal(icry,pos)) return icry;
    }

    return Crystal::invalidID_;
  }

  //-----------------------------------------------------------------------------
  // Position needs to be in the DiskFF frame
  bool Disk::isInsideCrystal(int icry, const CLHEP::Hep3Vector& pos) const
  {
    if (icry >= std::ssize(crystals_)) return false;

    const auto& cpos  = crystals_[icry].localPosition();
    const auto& csize = crystals_[icry].size();
    if (std::abs(cpos.x()-pos.x()) > 0.5*csize.x() ) return false;
    if (std::abs(cpos.y()-pos.y()) > 0.5*csize.y() ) return false;
    if (pos.z() < -1e-6 || pos.z() > csize.z()     ) return false;

    return true;
  }






  //-----------------------------------------------------------------------------
  //find the global indexes of the crystal neighbors for a given level (level = number of rings away)
  std::vector<int> Disk::neighbors(int crystalId, int level) const
  {
    std::vector<int> list = findLocalNeighbors(crystalId,level);
    for (auto& val : list) val += crystalOffset_;
    return list;
  }


  //-----------------------------------------------------------------------------
  //find the local indexes of the crystal neighbors for a given level (level = number of rings away)
  std::vector<int> Disk::findLocalNeighbors(int crystalId, int level) const
  {
    std::vector<int> list;
    std::vector<int> temp(crystalMap_->neighbors(crystalToMap_.at(crystalId),level));

    for (const auto& mapIdx : temp) {
      if (isCrystalIdxValid(mapIdx)) list.push_back(mapToCrystal_.at(mapIdx));
    }

    return list;
  }

  //-----------------------------------------------------------------------------
  bool Disk::isCrystalIdxValid(int i) const {
    return i < std::ssize(mapToCrystal_) && mapToCrystal_[i] != Crystal::invalidID_;
  }


  //-----------------------------------------------------------------------------
  //Move and rotate the disk. The displacement and rotation are relatives
  void Disk::moveDisk(const CLHEP::Hep3Vector& shift, const CLHEP::HepRotation& rot)
  {
    diskInfo_.origin(diskInfo_.origin()+shift);
    diskInfo_.originLocal(diskInfo_.originLocal()+shift);
    diskInfo_.frontFaceCenter(diskInfo_.frontFaceCenter()+shift);
    diskInfo_.backFaceCenter(diskInfo_.backFaceCenter()+shift);
    diskInfo_.rotation(diskInfo_.rotation()*rot);

    //do not forget to recalculate the global crystal position
    for (auto& crystal : crystals_){
       auto globalPosition = diskInfo_.origin() + diskInfo_.inverseRotation()*(crystal.localPosition());
       crystal.setPosition(globalPosition);
    }
  }

  //-----------------------------------------------------------------------------
  void Disk::moveCrystal(int id, const CLHEP::Hep3Vector& shift)
  {
    auto& crystal = crystals_.at(id-crystalOffset_);
    crystal.setLocalPosition(crystal.localPosition()+shift);
    auto globalPosition = diskInfo_.origin() + diskInfo_.inverseRotation()*(crystal.localPosition());
    crystal.setPosition(globalPosition);
  }


  //-----------------------------------------------------------------------------
  void Disk::print(std::ostream &os) const
  {
     os<<"Disk                  "<<id_<<std::endl;
     os<<"Number of crystals    "<<crystals_.size()<<std::endl;
     os<<"Crystal offset        "<<crystalOffset_<<std::endl;
     os<<"Radius In             "<<diskInfo_.innerEnvelopeR()<<std::endl;
     os<<"Radius Out            "<<diskInfo_.outerEnvelopeR()<<std::endl;
     os<<"origin                "<<diskInfo_.originLocal()<<std::endl;
     os<<"origin Mu2e           "<<diskInfo_.origin()<<std::endl;
     os<<"size                  "<<diskInfo_.size()<<std::endl;
     os<<"rotation              "<<diskInfo_.rotation()<<std::endl;
     os<<"originToCrystalOrigin "<<diskInfo_.originToCrystalOrigin()<<std::endl;
     os<<"z Front               "<<diskInfo_.frontFaceCenter().z()<<std::endl;
     os<<"z Back                "<<diskInfo_.backFaceCenter().z()<<std::endl;
     os<<"r In tracker          "<<diskInfo_.innerEnvelopeR()<<std::endl;
     os<<"r Out tracker         "<<diskInfo_.outerEnvelopeR()<<std::endl;
  }

}

