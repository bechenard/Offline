//
// Instance name information for TrkExt data product 
//
//  $Id: TrkExtInstanceName.hh,v 1.1 2012/10/23 16:26:40 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/10/23 16:26:40 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtInstanceName_HH
#define TrkExtInstanceName_HH

#include <string>
#include <vector>
#include "TrkBase/TrkParticle.hh"

namespace mu2e {

  namespace TrkExtParticleMass {
    // Mass of particles. Need to change to read from condition service. However the accuracy is no good for some particles (electron)
    const double MASS_ELECTRON = 0.510998910;
    const double MASS_MUON = 105.658367;
    const double MASS_PION = 139.57018;
    const double MASS_PROTON = 938.272013;
  };

  class TrkExtInstanceNameEntry {

  public:
    TrkExtInstanceNameEntry (TrkParticle::type _hepid, bool _updown, std::string _fitterName="");
    ~TrkExtInstanceNameEntry() {}
    std::string name;
    bool updown;
    int charge;
    double mass2;
    int hepid;
    std::string fitterName;
    unsigned int ntrk;
  };


  class TrkExtInstanceName {

  public:
    TrkExtInstanceName();
    ~TrkExtInstanceName() { }
    void clear (void) { _entries.clear(); }
    unsigned int size(void) const {return _entries.size();}
    TrkExtInstanceNameEntry & get(unsigned int i) { return _entries[i]; }
    bool updown (unsigned int i) const { return _entries[i].updown; }
    int charge (unsigned int i) const { return _entries[i].charge; }
    double mass2 (unsigned int i) const { return _entries[i].mass2; }
    int hepid (unsigned int i) const { return _entries[i].hepid; }
    std::string name (unsigned int i) const { return _entries[i].name; }
    std::string fitterName (unsigned int i) const { return _entries[i].fitterName; }
    void addTrack(unsigned int i) { _entries[i].ntrk += 1; }
    unsigned int ntrk (unsigned int i) const { return _entries[i].ntrk; }

//    void constructAll() ;
    void construct (TrkParticle::type _hepid, bool _up_down, std::string _fitterName) ;

  private:
    std::vector<TrkExtInstanceNameEntry> _entries;
  };



} // end namespace mu2e


#endif
