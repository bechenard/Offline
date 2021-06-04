#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/LHEParticle.hh"
#include <iostream>
#include <sstream>

namespace mu2e {

  LHEParticle::LHEParticle(const std::string& line)
  {
     std::stringstream ss(line);
     std::string word;
     std::vector<std::string> tokens;
     while (ss >> word) tokens.push_back(word);

     if (tokens.size() != 13) throw cet::exception("CATEGORY")<<"Wrong number of tokens in LHE particle record.";

     idup_       = atof(tokens[0].c_str());
     istup_      = atoi(tokens[1].c_str());
     mothup_[0]  = atoi(tokens[2].c_str());
     mothup_[1]  = atoi(tokens[3].c_str());
     icolup_[0]  = atoi(tokens[4].c_str());
     icolup_[1]  = atoi(tokens[5].c_str());
     pup_[0]     = atof(tokens[6].c_str());
     pup_[1]     = atof(tokens[7].c_str());
     pup_[2]     = atof(tokens[8].c_str());
     pup_[3]     = atof(tokens[9].c_str());
     pup_[4]     = atof(tokens[10].c_str());
     vtimup_     = atof(tokens[11].c_str());
     spinup_     = atof(tokens[12].c_str());
  }


  void LHEParticle::print(std::ostream& stream) const 
  {
    stream << "LHEParticle { "
           << "idup: " << idup() << ", istup: " << istup()
           << ", mothup[0]: " << mothup(0) << ", mothup[1]: " << mothup(1)
           << ", icolup[0]: " << icolup(0) << ", icolup[1]: " << icolup(1)
           << ", pup[0]: " << pup(0) << ", pup[1]: " << pup(1)
           << ", pup[2]: " << pup(2) << ", pup[3]: " << pup(3)
           << ", pup[4]: " << pup(4) << ", vtimup: " << vtimup()
           << ", spinup: " << spinup() << " }";
  }

}
