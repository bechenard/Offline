#ifndef ILCCLUSTER_H
#define ILCCLUSTER_H
/* Copyright(c) 2005-2006, ILC Project Experiment, All rights reserved.   *
 * See cxx source for full Copyright notice                               */

/* $Id: IlcCluster.h,v 1.1 2012/12/04 00:51:27 tassiell Exp $ */

//-------------------------------------------------------------------------
//                          Class IlcCluster
//                Generic class for clusters
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------

#include <TObject.h>

//_____________________________________________________________________________
class IlcCluster : public TObject {
public:
  IlcCluster();
  virtual ~IlcCluster() {;}
  IlcCluster(Int_t *lab, Float_t *hit);
  void SetLabel(Int_t lab, Int_t i) {fTracks[i]=lab;}

  virtual void SetX(Float_t /*x*/){}
  virtual void SetY(Float_t y)      {fY=y;}
  virtual void SetZ(Float_t z)      {fZ=z;}
  virtual void SetDetector(Int_t /*detector*/) {}
  void SetSigmaY2(Float_t sy2)      {fSigmaY2=sy2;}
  void SetSigmaZ2(Float_t sz2)      {fSigmaZ2=sz2;}

  Int_t GetLabel(Int_t i) const {return fTracks[i];}
  virtual Float_t GetX()          const {return 0;}
  Float_t GetY()          const {return fY;}
  Float_t GetZ()          const {return fZ;}
  virtual Int_t GetDetector()  const {return 0;}
  Float_t GetSigmaY2()    const {return fSigmaY2;}
  Float_t GetSigmaZ2()    const {return fSigmaZ2;}

  //PH  virtual void Use() = 0;
  //PH The pure virtual function has been temporarily replaced by 
  //PH virtual function with empty body. This correction somehow helps
  //PH to write/read TClonesArray with IlcVXDclusterV2 objects, but obviously
  //PH hides some more tricky problems (to be investigated)
  virtual void Use(Int_t = 0) {;}

protected:
  Int_t     fTracks[3];//labels of overlapped tracks
  Float_t   fY ;       //Y of cluster
  Float_t   fZ ;       //Z of cluster
  Float_t   fSigmaY2;  //Sigma Y square of cluster
  Float_t   fSigmaZ2;  //Sigma Z square of cluster
  
  ClassDef(IlcCluster,2)  // Time Projection Chamber clusters
};

#endif


