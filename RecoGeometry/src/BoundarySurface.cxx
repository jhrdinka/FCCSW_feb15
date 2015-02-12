//
//  BoundarySurface.h
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#include "RecoGeometry/BoundarySurface.h"

Reco::BoundarySurface::BoundarySurface() :
m_nextVolume(0),
m_previousVolume(0),
m_nextVolumes(0),
m_previousVolumes(0)
{}
/*
Reco::BoundarySurface::BoundarySurface(Volume* nextVolume, Volume* previousVolume) :
m_nextVolume(nextVolume),
m_previousVolume(previousVolume),
m_nextVolumes(0),
m_previousVolumes(0)
{}

Reco::BoundarySurface::BoundarySurface(VolumeArray* nextVolumes, VolumeArray* previousVolumes) :
m_nextVolume(0),
m_previousVolume(0),
m_nextVolumes(nextVolumes),
m_previousVolumes(previousVolumes)
{}
*/
Reco::BoundarySurface::~BoundarySurface()
{
    delete m_nextVolumes;
}

Reco::BoundarySurface& Reco::BoundarySurface::operator=(const BoundarySurface& boundarysurface)
{
    if(this!=&boundarysurface)
    {
        m_nextVolume     = boundarysurface.m_nextVolume;
        m_previousVolume = boundarysurface.m_previousVolume;
        delete m_nextVolumes;    m_nextVolumes = 0;
        delete m_previousVolumes; m_previousVolumes = 0;
    }
    return (*this);
    
}

void Reco::BoundarySurface::setNextVolume(std::shared_ptr<const Reco::Volume> nextVolume) const
{
    m_nextVolume = nextVolume;
}

void Reco::BoundarySurface::setPreviousVolume(std::shared_ptr<const Reco::Volume> previousVolume) const
{
    m_previousVolume = previousVolume;
}

void Reco::BoundarySurface::setNextVolumes(const VolumeArray* nextVolumes) const
{
    m_nextVolumes = nextVolumes;
}

void Reco::BoundarySurface::setPreviousVolumes(const VolumeArray* previousVolumes) const
{
    m_previousVolumes = previousVolumes;
}

const Reco::Volume* Reco::BoundarySurface::getNextVolume() const
{
    return (m_nextVolume.get());
}

const Reco::Volume* Reco::BoundarySurface::getPreviousVolume() const
{
    return (m_previousVolume.get());
}

const Reco::BoundarySurface::VolumeArray* Reco::BoundarySurface::getNextVolumes() const
{
    return (m_nextVolumes);
}

const Reco::BoundarySurface::VolumeArray* Reco::BoundarySurface::getPreviousVolumes() const
{
    return (m_previousVolumes);
}


