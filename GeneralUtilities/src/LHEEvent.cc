#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/LHEEvent.hh"
#include <iostream>
#include <sstream>

namespace mu2e {

  void LHEEvent::reset()
  {
     nup_ = idprup_ = 0;
     xwgtup_ = scalup_ = aqedup_ = aqcdup_ = vtxt_ = 0.0;
     vtx_    = std::array<double,3>{0.,0.,0.};
     particles_.clear();
  }


  void LHEEvent::newEvent(const std::string& line)
  {
     reset();

     std::stringstream ss(line);
     std::string word;
     std::vector<std::string> tokens;
     while (ss >> word) tokens.push_back(word);

     if (tokens.size() != 6) throw cet::exception("CATEGORY")<<"Wrong number of tokens in LHE event information record.";

     nup_    = atoi(tokens[0].c_str());
     idprup_ = atoi(tokens[1].c_str());
     xwgtup_ = atof(tokens[2].c_str());
     scalup_ = atof(tokens[3].c_str());
     aqedup_ = atof(tokens[4].c_str());
     aqcdup_ = atof(tokens[5].c_str());
  }


  void LHEEvent::addParticle(const std::string& line)
  {
     particles_.push_back(LHEParticle(line));
  }


  void LHEEvent::setVertex(double x, double y, double z) 
  {
     vtx_= std::array<double,3>{x,y,z};
  }

  void LHEEvent::setVertex(const std::string& line)
  {
      std::stringstream ss(line);
      std::string word;
      std::vector<std::string> tokens;
      while (ss >> word) tokens.push_back(word);

      if (tokens.size() != 4 && tokens.size() != 5) {
        throw cet::exception("CATEGORY")<<"Wrong number of tokens in LHE event vertex information";
      }

      vtx_[0] = atof(tokens[1].c_str());
      vtx_[1] = atof(tokens[2].c_str());
      vtx_[2] = atof(tokens[3].c_str());
      if (tokens.size() > 4) vtxt_ = atof(tokens[4].c_str()); 
  }

}
