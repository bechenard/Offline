#ifndef GeneralUtilities_LHEReader_hh
#define GeneralUtilities_LHEReader_hh

#include "GeneralUtilities/inc/LHEEvent.hh"
#include <fstream>


namespace mu2e {

   class LHEReader {

      public:
        LHEReader(const std::string& filename);
        virtual ~LHEReader();

        void  readNextEvent();
        const LHEEvent& event() const {return event_;}
        bool  isValid()         const {return valid_;}

      private:
	std::ifstream ifile_;
	LHEEvent      event_;
        bool          valid_;

   };

}
#endif
