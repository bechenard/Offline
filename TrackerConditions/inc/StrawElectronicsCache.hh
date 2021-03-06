#ifndef TrackerConditions_StrawElectronicsCache_hh
#define TrackerConditions_StrawElectronicsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkDelayPanel.hh"
#include "DbTables/inc/TrkPreampRStraw.hh"
#include "DbTables/inc/TrkPreampStraw.hh"
#include "DbTables/inc/TrkThresholdRStraw.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

#include "TrackerConditions/inc/StrawElectronicsMaker.hh"


namespace mu2e {
  class StrawElectronicsCache : public ProditionsCache {
  public: 
    StrawElectronicsCache(StrawElectronicsConfig const& config):
      ProditionsCache(StrawElectronics::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      _eventTiming_p = std::make_unique<ProditionsHandle<EventTiming> >();
      if(_useDb) {
	_tdp_p  = std::make_unique<DbHandle<TrkDelayPanel>>();
	_tprs_p = std::make_unique<DbHandle<TrkPreampRStraw>>();
	_tps_p  = std::make_unique<DbHandle<TrkPreampStraw>>();
	_ttrs_p = std::make_unique<DbHandle<TrkThresholdRStraw>>();
      }
    }
    
    set_t makeSet(art::EventID const& eid) {
      auto et = _eventTiming_p->get(eid);
      ProditionsEntity::set_t cids = et.getCids();
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdp_p->get(eid);
	_tprs_p->get(eid);
	_tps_p->get(eid);
	_ttrs_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tdp_p->cid());
	cids.insert(_tprs_p->cid());
	cids.insert(_tps_p->cid());
	cids.insert(_ttrs_p->cid());
      }
      return cids;
    }
    
    DbIoV makeIov(art::EventID const& eid) {
      _eventTiming_p->get(eid);
      DbIoV iov = _eventTiming_p->iov();
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdp_p->get(eid);
	_tprs_p->get(eid);
	_tps_p->get(eid);
	_ttrs_p->get(eid);
	// restrict the valid range ot the overlap
	iov.overlap(_tdp_p->iov());
	iov.overlap(_tprs_p->iov());
	iov.overlap(_tps_p->iov());
	iov.overlap(_ttrs_p->iov());
      }
      return iov;
    }
    
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      auto et = _eventTiming_p->getPtr(eid);
      if(_useDb) {
	return _maker.fromDb( _tdp_p->getPtr(eid),
			      _tprs_p->getPtr(eid), 
			      _tps_p->getPtr(eid),
			      _ttrs_p->getPtr(eid), 
                              et);
      } else {
	return _maker.fromFcl(et);
      }
    }
    
  private:
    bool _useDb;
    StrawElectronicsMaker _maker;


    // these handles are not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<EventTiming> > _eventTiming_p;
    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<TrkDelayPanel>> _tdp_p;
    std::unique_ptr<DbHandle<TrkPreampRStraw>> _tprs_p;
    std::unique_ptr<DbHandle<TrkPreampStraw>> _tps_p;
    std::unique_ptr<DbHandle<TrkThresholdRStraw>> _ttrs_p;

  };
};

#endif
