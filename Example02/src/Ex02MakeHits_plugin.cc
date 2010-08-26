/*----------------------------------------------------------------------

A producer module that makes a collection of overly simplified "hits"
and adds them to the event.

$Id: Ex02MakeHits_plugin.cc,v 1.11 2010/08/26 19:58:17 kutschke Exp $
$Author: kutschke $
$Date: 2010/08/26 19:58:17 $
   
Original author Rob Kutschke


----------------------------------------------------------------------*/

// C++ includes.
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Other external includes.
#include <boost/shared_ptr.hpp>

// Mu2e includes.
#include "ToyDP/inc/ToyHitCollection.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CTrackerGeom/inc/CTracker.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "GeneralUtilities/inc/pow.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

using namespace std;

namespace mu2e {

  static double const czmax = 0.5;

  class Ex02MakeHits : public edm::EDProducer {
    
  public:
    explicit Ex02MakeHits(edm::ParameterSet const& pSet) : 
      
      // Information from the parameter set.
      minPulseHeight_(pSet.getParameter<double>("minPulseHeight")),

      // Random number distributions share one engine.
      _randGauss( createEngine( get_seed_value(pSet)) ),
      _angularDistribution( _randGauss.engine(), -czmax, czmax )
    {
      produces<ToyHitCollection>();

    }
    virtual ~Ex02MakeHits() { }

    virtual void produce(edm::Event& e, edm::EventSetup const& c);
    
  private:

    // Parameters from run time configuration.
    double minPulseHeight_;

    CLHEP::RandGaussQ _randGauss;
    RandomUnitSphere  _angularDistribution;

  };

  void
  Ex02MakeHits::produce(edm::Event& event, edm::EventSetup const&) {

    // Geometry information about the CTracker.
    GeomHandle<CTracker> ctracker;

    // Extract the event number from the event.
    int eventNumber = event.id().event();

    // Make the collection to hold the hits that we will make.
    auto_ptr<ToyHitCollection> p(new ToyHitCollection);

    // Generate a unit vector, random over the unit sphere ( over a 
    // restricted range of cos(theta).  This will be used as the 
    // direction of the generated particle.
    CLHEP::Hep3Vector vunit(_angularDistribution.fire());

    // Radius of each layer of the tracker.
    // Access by reference so there is no copying.
    vector<double> const& radius = ctracker->r();

    // Loop over the layers of the CTracker.
    for ( vector<double>::size_type i=0; i<radius.size(); ++i){

      // Path length, s, from origin to intersection with this layer.
      double sinTheta = sqrt( 1. - square(vunit.cosTheta()) );
      double s = radius[i]/sinTheta;

      // Point of intersection between the track and the layer.
      CLHEP::Hep3Vector intersectionPoint = s*vunit;

      // Generate a gaussian distributed pulse height.
      double ph = _randGauss.fire(10.,2.);
      
      // Keep only hits with large enough pulse height.
      if ( ph > minPulseHeight_ ){
        p->push_back(ToyHit(intersectionPoint,ph));
      }
    }

    // Put the generated hits into the event.
    //event.put(p);

    // At this point p is no longer usable.  If you try to
    // use p it will generate a run time error.

    // There is an alternative syntax for put.
    // Comment out the call to put and uncomment the following lines:
    //
    // Return value is a "read only pointer" to the data product.
    
    edm::OrphanHandle<ToyHitCollection> q = event.put(p);
    edm::LogInfo("Hits") << "Number of hits: " 
                         << q->size();
  }

} // end namespace mu2e

using mu2e::Ex02MakeHits;
DEFINE_FWK_MODULE(Ex02MakeHits);
