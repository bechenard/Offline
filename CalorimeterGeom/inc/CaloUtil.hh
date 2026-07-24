//
// Contains geometry utilities: coordinates transformations and checks if you are inside disks
// FF is the abrevation of front face
//
// Original author B. Echenard
//

#ifndef CalorimeterGeom_CaloUtil_hh
#define CalorimeterGeom_CaloUtil_hh

#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"

#include <vector>
#include <memory>


namespace mu2e {

    class DiskCalorimeter;

    class CaloUtil {

       public:
         CaloUtil(const DiskCalorimeter& owner);

         void trackerCenter (const CLHEP::Hep3Vector& vec);
         const CLHEP::Hep3Vector& trackerCenter()  const;

         CLHEP::Hep3Vector mu2eToCrystal(unsigned crystalId, const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector mu2eToDisk   (unsigned diskId,    const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector mu2eToDiskFF (unsigned diskId,    const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector mu2eToTracker(const CLHEP::Hep3Vector& pos)                     const;

         CLHEP::Hep3Vector crystalToMu2e(unsigned crystalId, const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector diskToMu2e   (unsigned diskId,    const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector diskFFToMu2e (unsigned diskId,    const CLHEP::Hep3Vector& pos) const;
         CLHEP::Hep3Vector trackerToMu2e(const CLHEP::Hep3Vector& pos)                     const;


       private:
         const DiskCalorimeter& owner_;
         CLHEP::Hep3Vector      trackerCenter_;
     };
}

#endif
