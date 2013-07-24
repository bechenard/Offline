//
// Read the table with data about DIO spectrum, and
// merge the spectrum with the analytic expression
// in the endpoint region taken from Czarnecki spectrum
// Czarneckki et al 10.1103/PhysRevD.84.013006
//
// $Id: CzarneckiSpectrum.cc,v 1.10 2013/07/24 18:48:24 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/07/24 18:48:24 $
//

// Mu2e includes
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/PhysicsParams.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"

// Framework includes
#include "cetlib/pow.h"

// C++ includes
#include <fstream>
#include <iostream>
#include <vector>

namespace mu2e {

  CzarneckiSpectrum::CzarneckiSpectrum() :
    _table ( loadTable<2>( "ConditionsService/data/czarnecki_"+
                           GlobalConstantsHandle<PhysicsParams>()->getStoppingTarget()+".tbl" ) )
  {}

  double CzarneckiSpectrum::getWeight(double E) {

    const unsigned iRow = _table.findLowerBoundRow( E );
    double weight       = _table( iRow ) ;
    if ( iRow == _table.getNrows()-1 || iRow == 0 ) return weight;

    auto const & row        = _table.row( iRow   );
    auto const & row_before = _table.row( iRow-1 );
    auto const & row_after  = _table.row( iRow+1 );

    weight = interpolate( E, row_after, row, row_before );
    if ( weight < 0 ) weight = interpolateE5 ( E, row );

    return weight;
  }
  

  double CzarneckiSpectrum::interpolate(const double E, 
                                        const TableRow<2>& row_after,
                                        const TableRow<2>& row,
                                        const TableRow<2>& row_before ) {

    const double e1(  row_after[0] );  const double p1(  row_after[1] );
    const double e2(        row[0] );  const double p2(        row[1] );
    const double e3( row_before[0] );  const double p3( row_before[1] );

    const double discr = e1*e1*e2 + e1*e3*e3 + e2*e2*e3 - e3*e3*e2 - e1*e1*e3 - e1*e2*e2;

    const double A = (p1*e2 + p3*e1 + p2*e3 - p3*e2 - p1*e3 - p2*e1) / discr;

    const double B = (e1*e1*p2 + e3*e3*p1 + e2*e2*p3 - e3*e3*p2 - e1*e1*p3 - e2*e2*p1) / discr;

    const double C = (e1*e1*e2*p3 + e3*e3*e1*p2 + e2*e2*e3*p1 -
                      e3*e3*e2*p1 - e1*e1*e3*p2 - e2*e2*e1*p3) / discr;

    return (A*E*E + B*E + C);

  }
  
  double CzarneckiSpectrum::interpolateE5(double E, TableRow<2> val ) {
    
    GlobalConstantsHandle<PhysicsParams> phy;

    const double energy = val.at(0);
    const double weight = val.at(1);

    const double b = weight/cet::pow<5>( phy->getMuonEnergy()-energy-cet::square(energy)/(2*phy->getAtomicMass()));

    if ( phy->getEndpointEnergy() - E < 0 ) return 0.;
    else return b*cet::pow<5>( phy->getMuonEnergy() - E - cet::square(E)/(2*phy->getAtomicMass()) );
  }

}

