//
// Free function to construct the stopping targets.
//
//
// Original author Peter Shanahan
//
// Notes:


// C++ includes
#include <iostream>
#include <string>

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes
#include "Offline/Mu2eG4/inc/constructStoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/Mu2eG4/inc/StrawSD.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestCons.hh"
#include "Offline/Mu2eG4/inc/checkForOverlaps.hh"
#include "Offline/MECOStyleProtonAbsorberGeom/inc/MECOStyleProtonAbsorber.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Colour.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4LogicalVolumeStore.hh"

using namespace std;

namespace mu2e {

    VolumeInfo constructStoppingTarget( VolumeInfo   const& parent,
                                      SimpleConfig const& config ){

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "stoppingTarget", "stoppingTarget");

    const bool stoppingTargetIsVisible = geomOptions->isVisible("stoppingTarget");
    const bool stoppingTargetIsSolid   = geomOptions->isSolid("stoppingTarget");
    const bool forceAuxEdgeVisible     = geomOptions->forceAuxEdgeVisible("stoppingTarget");
    const bool doSurfaceCheck          = geomOptions->doSurfaceCheck("stoppingTarget");
    const bool placePV                 = geomOptions->placePV("stoppingTarget");

    bool const inGaragePosition = config.getBool("inGaragePosition",false);
    bool const OPA_IPA_ST_Extracted = (inGaragePosition) ? config.getBool("garage.extractOPA_IPA_ST") : false;
    double zOffGarage = (inGaragePosition && OPA_IPA_ST_Extracted) ? config.getDouble("garage.zOffset") : 0.;
    CLHEP::Hep3Vector relPosFake(0.,0., zOffGarage); //for offsetting target in garage position

    int verbosity(config.getInt("stoppingTarget.verbosity",0));

    if ( verbosity > 1 ) std::cout << "In constructStoppingTarget" << std::endl;
    // Master geometry for the Target assembly
    GeomHandle<StoppingTarget> target;

    Mu2eG4Helper    & _helper = *(art::ServiceHandle<Mu2eG4Helper>());
    AntiLeakRegistry & reg = _helper.antiLeakRegistry();


    TubsParams targetMotherParams(0., target->cylinderRadius(), target->cylinderLength()/2.);

    VolumeInfo targetInfo;
    std::string targetMotherName = "StoppingTargetMother";

    targetInfo = nestTubs(targetMotherName,
                            targetMotherParams,
                            findMaterialOrThrow(target->fillMaterial()),
                            0,
                            target->centerInMu2e() - parent.centerInMu2e() + relPosFake,
                            parent,
                            0,
                            false/*visible*/,
                            G4Colour::Black(),
                            false/*solid*/,
                            forceAuxEdgeVisible,
                            placePV,
                            doSurfaceCheck
                            );
    // now create the individual target foils

    G4VPhysicalVolume* pv;

    for (int itf=0; itf<target->nFoils(); ++itf) {

        TargetFoil foil=target->foil(itf);

        VolumeInfo foilInfo;
        G4Material* foilMaterial = findMaterialOrThrow(foil.material());

        std::ostringstream os;
        os << std::setfill('0') << std::setw(2) << itf;
        foilInfo.name = "Foil_" + os.str();

        if ( verbosity > 0 )  std::cout << __func__ << " " << foilInfo.name << std::endl;

        foilInfo.solid = new G4Tubs(foilInfo.name
                                    ,foil.rIn()
                                    ,foil.rOut()
                                    ,foil.halfThickness()
                                    ,0.
                                    ,CLHEP::twopi
                                    );

        foilInfo.logical = new G4LogicalVolume( foilInfo.solid
                                                , foilMaterial
                                                , foilInfo.name
                                                );


        // rotation matrix...
        G4RotationMatrix* rot = 0; //... will have to wait

        G4ThreeVector foilOffset(foil.centerInMu2e() - targetInfo.centerInMu2e() + relPosFake);
        if ( verbosity > 1 ) std::cout << "foil "
                                  << itf
                                  << " centerInMu2e="
                                  << foil.centerInMu2e()
                                  << ", offset="<< foilOffset<< std::endl;

        // G4 manages the lifetime of this object.
        pv = new G4PVPlacement( rot,
                                foilOffset,
                                foilInfo.logical,
                                "Target"+foilInfo.name,
                                targetInfo.logical,
                                0,
                                itf,
                                false);

        doSurfaceCheck && checkForOverlaps( pv, config, verbosity>0);

        if (!stoppingTargetIsVisible) {
          foilInfo.logical->SetVisAttributes(G4VisAttributes::GetInvisible());
        } else {
          G4VisAttributes* visAtt = reg.add(G4VisAttributes(true, G4Colour::Magenta()));
          visAtt->SetForceAuxEdgeVisible(config.getBool("g4.forceAuxEdgeVisible",false));
          visAtt->SetForceSolid(stoppingTargetIsSolid);
          foilInfo.logical->SetVisAttributes(visAtt);
        }
      }// target foils

    return targetInfo;
  }

} // end namespace mu2e
