#ifndef CalorimeterGeom_SquareMapper_hh
#define CalorimeterGeom_SquareMapper_hh
//
// Square position map generator:
//   tesselate a plane with squares starting from the center of the plane
//
//  original author : Bertrand Echenard (Caltech)
//
// Use basis vector, l and k, defined as
// l = right
// k = up
//
/*

 -----------------------
 |       |       |      |
 | -1 1  | 0 1   | 1 1  |
 |       |       |      |
 -----------------------
 |       |       |      |  l,k coordinates
 | -1 0  | 0 0   | 1 0  |
 |       |       |      |
 -----------------------
 |       |       |      |
 | -1 -1 | 0 -1  | 1 -1 |
 |       |       |      |
 -----------------------

  steps :  (1,0), (0,-1), (-1,0), (0,1) (clockwise from top left corner)
  segment: top=0, right=1, bottom=2,left=3

*/
// Tesselation algorithm: tessalate in "rings" from the center
//   for each ring, start at -l,+l (top left corner),
//   then go n time each step to create the ring
//
// Neighbors add (+-1,0) or (0,+-1)
// next ring of neighbours, add +(-2,2) and go around the ring,...
//
#include "Offline/CalorimeterGeom/inc/CrystalMapper.hh"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>


namespace mu2e {

  class SquLK {
    public:
      SquLK()             : l_(0),k_(0) {}
      SquLK(int l, int k) : l_(l),k_(k) {}

      void add(const SquLK &x) {l_+=x.l_;k_+=x.k_;}

      int l_;
      int k_;
  };


  class SquareMapper : public CrystalMapper {
    public:
      SquareMapper();

      int                    nCrystalMax    (int maxRing)                const override;
      CLHEP::Hep2Vector      xyFromIndex    (int thisIndex)              const override;
      int                    indexFromXY    (double x, double y)         const override;
      int                    indexFromRowCol(int nRow, int nCol)         const override;
      int                    rowFromIndex   (int thisIndex)              const override;
      int                    colFromIndex   (int thisIndex)              const override;
      int                    numNeighbors   (int level)                  const override;
      std::vector<int>       neighbors      (int thisIndex, int level=1) const override;
      const std::vector<double>& apexX() const override {return apexX_;}
      const std::vector<double>& apexY() const override {return apexY_;}

    private:
      SquLK    lk(int index)          const;
      int      index(const SquLK& lk) const;
      int      ring(const  SquLK& lk) const;

      std::vector<SquLK>   step_;
      std::vector<double>  apexX_;
      std::vector<double>  apexY_;
  };
}


#endif

