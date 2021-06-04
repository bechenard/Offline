#ifndef GeneralUtilities_LHEParticle_hh
#define GeneralUtilities_LHEParticle_hh

#include <string>
#include <array>
#include <vector>

namespace mu2e {

   class LHEParticle {

      public:
	 LHEParticle() {};
	 explicit LHEParticle(const std::string& data);
	 ~LHEParticle() = default;

	 int    idup()        const { return idup_; }
	 int    istup()       const { return istup_; }
	 int    mothup(int i) const { return mothup_.at(i); }
	 int    icolup(int i) const { return icolup_.at(i); }
	 double pup(int i)    const { return pup_.at(i); }
	 double vtimup()      const { return vtimup_; }
	 double spinup()      const { return spinup_; }

	 void   print(std::ostream& stream) const;


      private:
	 int                   idup_   = 0;
	 int                   istup_  = 0; 
	 std::array<int,2>     mothup_ = {{-1,-1}};
	 std::array<int,2>     icolup_ = {{-1,-1}};
	 std::array<double,5>  pup_    = {{0.,0.,0.,0.,0.}} ;
	 double                vtimup_ = 0;
	 int                   spinup_ = 0;
   };

}

#endif
