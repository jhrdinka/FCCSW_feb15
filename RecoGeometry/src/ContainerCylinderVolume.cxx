//
//  ContainerCylinderVolume.cxx
//  
//
//  Created by Julia Hrdinka on 27/01/15.
//
//

#include "RecoGeometry/ContainerCylinderVolume.h"

Reco::ContainerCylinderVolume::ContainerCylinderVolume(VolumeArray* volumes, std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax, double halfZ) :
ContainerVolume(volumes),
CylinderVolume(transf, Rmin, Rmax, halfZ)
{}

Reco::ContainerCylinderVolume::ContainerCylinderVolume(VolumeArray* volumes, TGeoNode* node, TGeoConeSeg* tube) :
ContainerVolume(volumes),
CylinderVolume(node,tube)
{}

Reco::ContainerCylinderVolume::~ContainerCylinderVolume()
{}

bool Reco::ContainerCylinderVolume::isInside(const Alg::Point3D& glopos, double tol) const
{
    return (CylinderVolume::isInside(glopos,tol));
}