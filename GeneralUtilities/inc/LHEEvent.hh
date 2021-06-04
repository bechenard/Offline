#ifndef GeneralUtilities_LHEEvent_hh
#define GeneralUtilities_LHEEvent_hh

#include "GeneralUtilities/inc/LHEParticle.hh"
#include <vector>
#include <array>


namespace mu2e {

   class LHEEvent {

     public:
	 LHEEvent()  = default;
	 ~LHEEvent() = default;

	 void newEvent(const std::string& data);
	 void addParticle(const std::string& line);
	 void setVertex(const std::string& line);
	 void setVertex(double x, double y, double z);
	 void reset();

	 int    nup()                                 const { return nup_; }
	 int    idprup()                              const { return idprup_; }
	 double wgtup()                               const { return xwgtup_; }
	 double scalup()                              const { return scalup_; }
	 double aqedup()                              const { return aqedup_; }
	 double aqcdup()                              const { return aqcdup_; }
	 double vertexTime()                          const { return vtxt_; }
	 const std::array<double,3>&     vertex()     const { return vtx_; }
	 const std::vector<LHEParticle>& particles()  const { return particles_; }

     private:
	 int                      nup_       = 0;
	 int                      idprup_    = 0;
	 double                   xwgtup_    = 0;
	 double                   scalup_    = 0;
	 double                   aqedup_    = 0;
	 double                   aqcdup_    = 0;
	 double                   vtxt_      = 0;
	 std::array<double,3>     vtx_       = {{0,0,0}};
	 std::vector<LHEParticle> particles_ = {};
   };
}

#endif

