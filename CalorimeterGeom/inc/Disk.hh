#ifndef CalorimeterGeom_Disk_hh
#define CalorimeterGeom_Disk_hh
//
// Hold information about a disk in the calorimter.
//
// The crystal numbering scheme start at the center of the disk
// to facilitate navigation
//
// Original author B Echenard
//
#include "Offline/CalorimeterGeom/inc/DiskInfo.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <vector>
#include <memory>

namespace mu2e {

  class Disk {
    public:
        Disk(int id, double rCrystalIn, double rCrystalOut, double nominalCellSize,
             double nominalCellLength, int offset, const CLHEP::Hep3Vector& diskOriginToCrystalOrigin,
             const std::string& crystalFileName);

        int                             id()           const {return id_;}

        size_t                          nCrystals()    const {return crystals_.size();}
        const Crystal&                  crystal(int i) const {return crystals_.at(i);}
              Crystal&                  crystal(int i)       {return crystals_.at(i);}
        int                             crystalOffset()const {return crystalOffset_;}

        const DiskInfo&                 diskInfo()     const {return diskInfo_;}
              DiskInfo&                 diskInfo()           {return diskInfo_;}

        std::vector<int>                neighbors      (int crystalId, int level=1)   const;
        std::vector<int>                idxFromRow     (int thisRow)                  const;
        int                             idxFromPosition(const CLHEP::Hep3Vector& pos) const;

        bool                            isInsideDisk   (const CLHEP::Hep3Vector& pos) const;
        bool                            isInsideCrystal(const CLHEP::Hep3Vector& pos) const;

        void                            moveCrystal    (int id, const CLHEP::Hep3Vector& disp);
        void                            moveDisk       (const CLHEP::Hep3Vector& disp,
                                                        const CLHEP::HepRotation& rotation);

        void                            print          (std::ostream& os = std::cout) const;



    private:
        void                            fillCrystals      (const CLHEP::Hep3Vector&, double nominalCellLength,
                                                           const std::string& filename);
        void                            checkPosition     (float xPos,float yPos,float xside,float ysize) const;
        bool                            isInsideCrystal   (int icry, const CLHEP::Hep3Vector& pos)        const;
        bool                            isCrystalIdxValid (int i)                                         const;
        std::vector<int>                findLocalNeighbors(int crystalId, int level=1)                    const;

        int                             id_;
        std::vector<Crystal>            crystals_;
        DiskInfo                        diskInfo_;
        double                          radiusInCrystal_;
        double                          radiusOutCrystal_;
        double                          nominalCellSize_;
        int                             crystalOffset_;
        std::shared_ptr<CrystalMapper>  crystalMap_;
        std::vector<int>                mapToCrystal_;
        std::vector<int>                crystalToMap_;
  };
}
#endif
