//
// Object to perform helix fit to straw hits
//
// $Id: CalHelixFinderAlg.hh,v 1.8 2014/05/18 13:56:50 murat Exp $
// $Author: murat $
// $Date: 2014/05/18 13:56:50 $
//
#ifndef CalHelixFinderAlg_HH
#define CalHelixFinderAlg_HH

// data
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
//#include "TString.h"

#include "CalPatRec/inc/LsqSums2.hh"
#include "CalPatRec/inc/LsqSums4.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "CalPatRec/inc/CalHelixPoint.hh"
#include "CalPatRec/inc/CalHelixFinderData.hh"

class TH1F;

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  class Calorimeter;
  class TTracker;
//-----------------------------------------------------------------------------
// output struct
//-----------------------------------------------------------------------------
  class CalHelixFinderAlg {
  public:
    enum { kMaxNHits = 10000 } ;
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    const TTracker*            _tracker;
    const Calorimeter*         _calorimeter;

    const CalTimePeak*         fTimePeak;
    const TimeCluster*         fTimeCluster; //needed for debugging
    
    double                     fCaloTime;
    double                     fCaloX;   
    double                     fCaloY;   
    double                     fCaloZ;   
    
    std::vector<CalHelixPoint> _xyzp;        // normally includes only hits from the time peak
//-----------------------------------------------------------------------------
// for diagnostics purposes save several states of _xyzp (only if _diag > 0)
//
// [0]: copy of initial _xyzp state
// [1]: after filterDist
// [2]: doPatternRecognition: after rescueHitsBeforeSeed
// [3]: doPatternRecognition: after rescueHits  (inside doPatternRecognition)
// [4]: doPatternRecognition: after final findZ (inside doPatternRecognition)
// [5]: doPatternRecognition: after doPatternRecognition
//
// each XYZPHack point has an 'outlier' flag telling if it belongs to the track
//-----------------------------------------------------------------------------
    struct SaveResults_t {
      std::vector<CalHelixPoint> _xyzp;    // points with used flags
      CalHelixFinderData         _helix;   // helix, found at this stage
    };

    SaveResults_t        _results[6];      // diagnostic buffers

    int                  fSeedIndex;
    int                  fCandIndex;   // index of the hit candiate 
    int                  fLastIndex;
    int                  fUseDefaultDfDz;

    int                  _diag;
    int                  _debug;
    int                  _debug2;
    int                  _smartTag;     //flag used to test addiotional layer of rejection after the search for the "best triplet"
    StrawHitFlag         _hsel;         // good hit selection
    StrawHitFlag         _bkgsel;       // background hit selection
    double               _maxElectronHitEnergy;
    int                  _minnhit;      // minimum # of hits to work with parameters for AGE center determination
    double               _maxdz;        // stereo selection parameters
                                        // 2014-03-10 Gianipez and P. Murat: limit
                                        // the dfdz value in the pattern-recognition stage
    double               _mpDfDz;
    int                  _minNSt;       // minimum number of active stations found in the ::findDfDZ(...) function
    double               _dzOverHelPitchCut; //cut on the ratio between the Dz and the predicted helix-pitch used in ::findDfDz(...)
    double               _maxDfDz;
    double               _minDfDz;
    double               _sigmaPhi;     // hit phi resolution (wrt the trajectory axis, assume R=25-30cm)
    double               _weightXY;     // scale factor for makeing the xy-chi2 with a mean close to 1
    double               _weightZPhi;
    double               _weight3D;
    double               _ew;           // error along the wire (mm)
    double               _maxXDPhi;     // max normalized hit residual in phi (findRZ)

    double               _hdfdz;        // estimated d(phi)/dz value
    double               _sdfdz;        // estimated d(phi)/dz error
    double               _hphi0;
					// 201-03-31 Gianipez added for changing the value of the
					// squared distance requed bewtween a straw hit and its predicted 
					// position used in the patter recognition procedure
    double               _distPatRec;

    double               _rhomin;
    double               _rhomax;        // crude cuts on tranvservse radius for stereo hits
    double               _mindist;       // minimum distance between points used in circle initialization
    double               _maxdist;       // maximum distance in hits
    double               _pmin, _pmax;   // range of total momentum
    double               _tdmin, _tdmax; // range of abs(tan(dip)
    double               _rcmin,_rcmax;  // maximum transverse radius of circle
    double               _sfactor;       // stereo hit error factor
    bool                 _xyweights;
    bool                 _zweights;      // weight points by estimated errors
    bool                 _filter;        // filter hits
    bool                 _plotall;       // plot also failed fits
    bool                 _usetarget;     // constrain to target when initializing
    mutable double       _bz;            // cached value of Field Z component at the tracker origin
//-----------------------------------------------------------------------------
// cached value of radius and pitch sign: these depend on the particle type
// and direction
//-----------------------------------------------------------------------------
    double    _rmin, _rmax, _smin, _smax, _dfdzsign;
//-----------------------------------------------------------------------------//
// store the paramters value of the most reliable track candidate
//-----------------------------------------------------------------------------//
    double    _x0, _y0, _phi0, _radius, _dfdz;
    int       _goodPointsTrkCandidate;
    int       _minPointsTrkCandidate;
    double    _chi2TrkCandidate;
    double    _maxChi2TrkCandidate;
    int       _markCandidateHits;        // apparently, always set to 0
                                         // thresholds for the worst hit chi2, total XY and ZPhi fit chi2's
    double    _hitChi2Max;
    double    _chi2xyMax;
    double    _chi2zphiMax;
    double    _chi2hel3DMax;

    // indices, distance from prediction and distance along z axis from the seeding hit
    // of the hits found in the pattern recognition

    int       _indicesTrkCandidate[kMaxNHits];
    double    _distTrkCandidate   [kMaxNHits];
    double    _dzTrkCandidate     [kMaxNHits];

    double    _phiCorrected       [kMaxNHits];
    int       _phiCorrectedDefined;


    double    _dfdzErr;                 // error on dfdz by ::findDfDz

    TH1F*     _hDist;
    double    _chi2nFindZ;
    double    _eventToLook;
    TH1F*     _hDfDzRes;
    TH1F*     _hPhi0Res;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
                                        // parameter set should be passed in on construction

    explicit CalHelixFinderAlg(fhicl::ParameterSet const&);
    virtual ~CalHelixFinderAlg();
                                        // cached bfield accessor
    double bz() const;

    void   calculateDfDz    (double phi0, double phi1, double z0,  double z1, double &dfdz);
    void   calculateDphiDz_2(const int* HitIndex, int NHits, double X0, double Y0, double& DphiDz);

    //projects the straw hit error along the radial direction of the circle-helix
    double calculateWeight     (const CLHEP::Hep3Vector& HitPos, 
				const CLHEP::Hep3Vector& StrawDir, 
				const CLHEP::Hep3Vector& HelCenter, 
				double                   Radius);

    double calculatePhiWeight  (const CLHEP::Hep3Vector& HitPos   , 
				const CLHEP::Hep3Vector& StrawDir , 
				const CLHEP::Hep3Vector& HelCenter, 
				double                   Radius   , 
				int                      Print    , 
				const char*              Banner=NULL);

    //calculates the residual along the radial direction of the helix-circle
    double calculateRadialDist (const CLHEP::Hep3Vector& HitPos, 
				const CLHEP::Hep3Vector& HelCenter, 
				double                   Radius);

    void   calculateTrackParameters(const CLHEP::Hep3Vector& p1, 
				    const CLHEP::Hep3Vector& p2,
                                    const CLHEP::Hep3Vector& p3,
				    CLHEP::Hep3Vector&       Center, 
				    double&                  Radius,
                                    double&                  Phi0, 
				    double&                  TanLambda);

				    // CalHelixFinderData&      mytrk,
                                    // bool                     cleanPattern=false);

    static double deltaPhi   (double phi1, double phi2);

   // returns the index of the hit which provides the highest contribute to the chi2
    void   cleanUpWeightedCircleFit(::LsqSums4&         TrkSxy,
				    int                SeedIndex,
				    int*               IdVec,
				    // CLHEP::Hep3Vector& HelCenter,
				    // double&            Radius,
				    double*            Weights,
				    int&               Iworst);

    bool   doLinearFitPhiZ     (CalHelixFinderData& Helix, 
				int                 SeedIndex, 
				int*                indexVec,
				int                 UseInteligentWeight=0, 
				int                 DoCleanUp =1);

   //perfoms the weighted circle fit, update the helix parameters (HelicCenter, Radius) and
    // fills the vector Weights which holds the calculated weights of the hits
    void   doWeightedCircleFit (::LsqSums4&        TrkSxy, 
				int                SeedIndex,
				int*               IdVec,
                                CLHEP::Hep3Vector& HelCenter, 
				double&            Radius, 
				double*            Weights,
                                int                Print=0, 
				const char*        Banner=NULL);

    void doPatternRecognition(CalHelixFinderData& mytrk);
  
    void defineHelixParams(CalHelixFinderData& Helix) const;

    // TH1F* hDist() {return _hDist;}

    int   isHitUsed(int index);

    void fillHitLayer                 (CalHelixFinderData& Helix);
    void fillXYZP                     (CalHelixFinderData& Helix);
    void filterDist                   ();
    void filterUsingPatternRecognition(CalHelixFinderData& Helix);
    bool findHelix                    (CalHelixFinderData& Helix);
    bool findHelix                    (CalHelixFinderData& Helix, const CalTimePeak* TimePeak);
    bool findHelix                    (CalHelixFinderData& Helix, const TimeCluster* TimePeak );
    int  findDfDz                     (CalHelixFinderData& Helix, int SeedIndex, int *indexVec, int  Diag_flag=0);
    void findTrack                    (int                 SeedIndex,
			               CalHelixFinderData& Helix,
			               bool                UseDefaultDfDz = false,
			               int                 UseMPVdfdz     = 0);
//-----------------------------------------------------------------------------
// setters
//-----------------------------------------------------------------------------
    void  setTracker    (const TTracker*    Tracker) { _tracker     = Tracker; }
    void  setCalorimeter(const Calorimeter* Cal    ) { _calorimeter = Cal    ; }
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    void   plotXY               (int ISet);

    void   plotZPhi             (int ISet);
    void   printInfo            (CalHelixFinderData& Helix);
    void   printXYZP            (CalHelixFinderData& Helix);

    int    refineHelixParameters(CalHelixFinderData& Helix,
				 int seedIndex, int *indexVec,
				 const char*  Banner=NULL,
				 int Print=0);

                                        // 12-10-2013 Gianipez: new pattern recognition functions
    void   rescueHitsBeforeSeed (CalHelixFinderData&  Helix);

    void   rescueHits           (CalHelixFinderData&  Helix, int   seedIndex       ,
				 int *indexVec             , int UsePhiResiduals = 0);

    void   resolve2PiAmbiguity  (const CLHEP::Hep3Vector& Center, double DfDz, double Phi0);

    void   resetTrackParamters  ();
//-----------------------------------------------------------------------------
// save intermediate results in diagnostics mode
//-----------------------------------------------------------------------------
    void   saveResults                    (std::vector<CalHelixPoint>& Xyzp, 
					   CalHelixFinderData&         Helix, 
					   int                         Index);

    void   searchWorstHitWeightedCircleFit(int                SeedIndex,
                                           int*               IdVec,
                                           const CLHEP::Hep3Vector& HelCenter,
                                           double&            Radius,
                                           double*            Weights,
                                           int&               Iworst ,
                                           double&            HitChi2Worst);

  };
}
#endif