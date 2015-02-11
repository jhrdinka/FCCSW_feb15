//
//  BoundaryDiscSurface.cxx
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#include "RecoGeometry/BoundaryDiscSurface.h"

/*Reco::BoundaryDiscSurface::BoundaryDiscSurface(Volume* nextVolume, Volume* previousVolume) :
BoundarySurface(nextVolume, previousVolume)
{}

Reco::BoundaryDiscSurface::BoundaryDiscSurface(VolumeArray* nextVolumes, VolumeArray* previousVolumes) :
BoundarySurface(nextVolumes, previousVolumes)
{}
*/

Reco::BoundaryDiscSurface::BoundaryDiscSurface(TGeoNode* node, TGeoConeSeg* tube) :
BoundarySurface(),
DiscSurface(node,tube)
{}

Reco::BoundaryDiscSurface::BoundaryDiscSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf) :
BoundarySurface(),
DiscSurface(tube, transf)
{}

Reco::BoundaryDiscSurface::BoundaryDiscSurface(std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax) :
BoundarySurface(),
DiscSurface(transf,Rmin,Rmax)
{}

Reco::BoundaryDiscSurface::~BoundaryDiscSurface()
{}

//const Volume* Reco::BoundaryDiscSurface::getNextVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const
//{
    
//}
/*
void Reco::BoundaryDiscSurface::setNextVolume(Volume* nextVolume) const
{
    Reco::BoundarySurface::setNextVolume(nextVolume);
}

void Reco::BoundaryDiscSurface::setPreviousVolume(Volume* previousVolume) const
{
    Reco::BoundarySurface::setPreviousVolume(previousVolume);
}

void Reco::BoundaryDiscSurface::setNextVolumes(VolumeArray* nextVolumes) const
{
    Reco::BoundarySurface::setNextVolumes(nextVolumes);
}

void Reco::BoundaryDiscSurface::setPreviousVolumes(VolumeArray* previousVolumes) const
{
    Reco::BoundarySurface::setPreviousVolumes(previousVolumes);
}
*/
const Reco::DiscSurface* Reco::BoundaryDiscSurface::surfaceRepresentation() const
{
    return (this);
}