//
//  BoundaryCylinderSurface.cxx
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#include "RecoGeometry/BoundaryCylinderSurface.h"

Reco::BoundaryCylinderSurface::BoundaryCylinderSurface(TGeoNode* node, TGeoConeSeg* tube) :
BoundarySurface(),
CylinderSurface(node,tube)
{}

Reco::BoundaryCylinderSurface::BoundaryCylinderSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf) :
BoundarySurface(),
CylinderSurface(tube,transf)
{}

Reco::BoundaryCylinderSurface::BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D> transf, double radius, double halfZ) :
BoundarySurface(),
CylinderSurface(transf,radius,halfZ)
{}

Reco::BoundaryCylinderSurface::~BoundaryCylinderSurface()
{}

Reco::BoundaryCylinderSurface& Reco::BoundaryCylinderSurface::operator=(const BoundaryCylinderSurface& boundarycylindersurface)
{
    if (this!=&boundarycylindersurface) {
        Reco::BoundarySurface::operator=(boundarycylindersurface);
        Reco::Surface::operator=(boundarycylindersurface);
    }
    return (*this);
}
/*
void Reco::BoundaryCylinderSurface::setNextVolume(Volume* nextVolume) const
{
    Reco::BoundarySurface::setNextVolume(nextVolume);
}

void Reco::BoundaryCylinderSurface::setPreviousVolume(Volume* previousVolume) const
{
    Reco::BoundarySurface::setPreviousVolume(previousVolume);
}

void Reco::BoundaryCylinderSurface::setNextVolumes(VolumeArray* nextVolumes) const
{
    Reco::BoundarySurface::setNextVolumes(nextVolumes);
}

void Reco::BoundaryCylinderSurface::setPreviousVolumes(VolumeArray* previousVolumes) const
{
    Reco::BoundarySurface::setPreviousVolumes(previousVolumes);
}
*/
const Reco::CylinderSurface* Reco::BoundaryCylinderSurface::surfaceRepresentation() const
{
    return (this);
} 