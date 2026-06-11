#ifndef CaloConditions_AlignedCalConfig_hh
#define CaloConditions_AlignedCalConfig_hh
//
// Initialize CaloDAQMap from fcl
//
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include <string>

namespace mu2e {

  struct AlignedCalConditions {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    fhicl::Atom<std::string> filenameDisk {
      Name("filenameDisk"), Comment("Filename for disk alignment file")};
    fhicl::Atom<std::string> filenameCrystal{
      Name("filenameCrystal"), Comment("Filename for crystal alignment file")};
    fhicl::Atom<int> verbose{
      Name("verbose"), Comment("verbosity")};
    fhicl::Atom<bool> useDb{
      Name("useDb"), Comment("use database or fcl")};
  };

}

#endif
