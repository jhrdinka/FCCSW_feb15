//
//  ContainerCylinderVolume.cxx
//  
//
//  Created by Julia Hrdinka on 27/01/15.
//
//

#include "RecoGeometry/ContainerCylinderVolume.h"

Reco::ContainerCylinderVolume::ContainerCylinderVolume(VolumeArray* volumes, std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax, double halfZ) :
Reco::ContainerVolume(volumes),
Reco::CylinderVolume(transf, Rmin, Rmax, halfZ)
{}

Reco::ContainerCylinderVolume::ContainerCylinderVolume(VolumeArray* volumes, TGeoNode* node, TGeoConeSeg* tube) :
Reco::ContainerVolume(volumes),
Reco::CylinderVolume(node,tube)
{}

Reco::ContainerCylinderVolume::ContainerCylinderVolume(const ContainerCylinderVolume& containercylindervolume) :
Reco::Volume(containercylindervolume),
Reco::ContainerVolume(containercylindervolume),
Reco::CylinderVolume(containercylindervolume)
{}

Reco::ContainerCylinderVolume* Reco::ContainerCylinderVolume::clone() const
{
    return (new Reco::ContainerCylinderVolume(*this));
}

Reco::ContainerCylinderVolume::~ContainerCylinderVolume()
{}

bool Reco::ContainerCylinderVolume::isInside(const Alg::Point3D& glopos, double tol) const
{
    return (CylinderVolume::isInside(glopos,tol));
}

bool Reco::ContainerCylinderVolume::setBoundarySurface(size_t n, std::shared_ptr<const BoundarySurface> boundarysurface) const
{
    return (Reco::Volume::setBoundarySurface(n,boundarysurface));
}