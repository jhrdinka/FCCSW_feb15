//
//  ContainerVolume.cxx
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#include "RecoGeometry/ContainerVolume.h"


Reco::ContainerVolume::ContainerVolume(VolumeArray* volumes) :
Reco::Volume(),
m_volumes(volumes)
{}

Reco::ContainerVolume::ContainerVolume(const ContainerVolume& containervolume):
Reco::Volume(containervolume),
m_volumes(0)
{
    m_volumes = containervolume.m_volumes;
}

Reco::ContainerVolume::~ContainerVolume()
{
    m_volumes = 0;
}

bool Reco::ContainerVolume::setBoundarySurface(size_t n, std::shared_ptr<const BoundarySurface> boundarysurface) const
{
    return (Reco::Volume::setBoundarySurface(n,boundarysurface));
}

const Reco::Volume* Reco::ContainerVolume::getVolume(const Alg::Point3D& glopos) const
{
    return (m_volumes->object(glopos));
}