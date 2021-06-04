#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/LHEReader.hh"
#include <iostream>


namespace mu2e {

  LHEReader::LHEReader(const std::string& filename) : valid_(true)
  {
     ifile_.open(filename.c_str(), std::ifstream::in);
     if (!ifile_.is_open()) throw cet::exception("CATEGORY")<< "No LHE File "<<filename;
     std::cout << "Opened LHE file " << filename << std::endl;
  }

  LHEReader::~LHEReader()
  { 
     ifile_.close();
  }


  void LHEReader::readNextEvent()
  {
     event_.reset();

     std::string line;
     bool foundEventElement(false);
     while (getline(ifile_, line))
     {
	 if (line == "<event>") {
            foundEventElement = true;
            break;
	 }
     }

     if (!foundEventElement) {valid_=false; return;}

     getline(ifile_, line);
     event_.newEvent(line);

     while (getline(ifile_, line))
     {
	 if (line == "</event>") break;

	 if (line.find("#") == std::string::npos) event_.addParticle(line); 
	 else if (line.find("#vertex") != std::string::npos) event_.setVertex(line);
     }
  }

}
