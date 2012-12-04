#ifndef ILCTRACKPOINTARRAY_H
#define ILCTRACKPOINTARRAY_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                          Class IlcTrackPoint                             //
//   This class represent a single track space-point.                       //
//   It is used to access the points array defined in IlcTrackPointArray.   //
//   Note that the space point coordinates are given in the global frame.   //
//                                                                          //
//   cvetan.cheshkov@cern.ch 3/11/2005                                      //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class IlcTrackPoint : public TObject {

 public:
  IlcTrackPoint();
  IlcTrackPoint(Float_t x, Float_t y, Float_t z, const Float_t *cov, UShort_t volid);
  IlcTrackPoint(const Float_t *xyz, const Float_t *cov, UShort_t volid);
  IlcTrackPoint(const IlcTrackPoint &p);
  IlcTrackPoint& operator= (const IlcTrackPoint& p);
  virtual ~IlcTrackPoint() {}

  // Multiplication with TGeoMatrix and distance between points (chi2) to be implemented

  void     SetXYZ(Float_t x, Float_t y, Float_t z, const Float_t *cov = 0);
  void     SetXYZ(const Float_t *xyz, const Float_t *cov = 0);
  void     SetVolumeID(UShort_t volid) { fVolumeID = volid; }

  Float_t  GetX() const { return fX; }
  Float_t  GetY() const { return fY; }
  Float_t  GetZ() const { return fZ; }
  void     GetXYZ(Float_t *xyz, Float_t *cov = 0) const;
  const Float_t *GetCov() const { return &fCov[0]; }
  UShort_t GetVolumeID() const { return fVolumeID; }

  Float_t  GetResidual(const IlcTrackPoint &p, Bool_t weighted = kFALSE) const;

  Float_t  GetAngle() const;
  IlcTrackPoint& Rotate(Float_t alpha) const;
  IlcTrackPoint& MasterToLocal() const;

  void     Print(Option_t *) const;

 private:

  Float_t  fX;        // X coordinate
  Float_t  fY;        // Y coordinate
  Float_t  fZ;        // Z coordinate
  Float_t  fCov[6];   // Cov matrix
  UShort_t fVolumeID; // Volume ID

  ClassDef(IlcTrackPoint,1)
};

//////////////////////////////////////////////////////////////////////////////
//                          Class IlcTrackPointArray                        //
//   This class contains the ESD track space-points which are used during   //
//   the alignment procedures. Each space-point consist of 3 coordinates    //
//   (and their errors) and the index of the sub-detector which contains    //
//   the space-point.                                                       //
//   cvetan.cheshkov@cern.ch 3/11/2005                                      //
//////////////////////////////////////////////////////////////////////////////

class IlcTrackPointArray : public TObject {

 public:

  IlcTrackPointArray();
  IlcTrackPointArray(Int_t npoints);
  IlcTrackPointArray(const IlcTrackPointArray &array);
  IlcTrackPointArray& operator= (const IlcTrackPointArray& array);
  virtual ~IlcTrackPointArray();

  //  Bool_t    AddPoint(Int_t i, IlcCluster *cl, UShort_t volid);
  Bool_t    AddPoint(Int_t i, const IlcTrackPoint *p);

  Int_t     GetNPoints() const { return fNPoints; }
  Int_t     GetCovSize() const { return fSize; }
  Bool_t    GetPoint(IlcTrackPoint &p, Int_t i) const;
  // Getters for fast access to the coordinate arrays
  const Float_t*  GetX() const { return &fX[0]; }
  const Float_t*  GetY() const { return &fY[0]; }
  const Float_t*  GetZ() const { return &fZ[0]; }
  const Float_t*  GetCov() const { return &fCov[0]; }
  const UShort_t* GetVolumeID() const { return &fVolumeID[0]; }

  Bool_t    HasVolumeID(UShort_t volid) const;

 private:

  Int_t     fNPoints;    // Number of space points
  Float_t   *fX;         //[fNPoints] Array with space point X coordinates
  Float_t   *fY;         //[fNPoints] Array with space point Y coordinates
  Float_t   *fZ;         //[fNPoints] Array with space point Z coordinates
  Int_t     fSize;       // Size of array with cov matrices = 6*N of points
  Float_t   *fCov;       //[fSize] Array with space point coordinates cov matrix
  UShort_t  *fVolumeID;  //[fNPoints] Array of space point volume IDs

  ClassDef(IlcTrackPointArray,1)
};

#endif
