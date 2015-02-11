//
//  ContainerVolume.cxx
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#include "RecoGeometry/ContainerVolume.h"


Reco::ContainerVolume::ContainerVolume(VolumeArray* volumes) :
Volume(),
m_volumes(volumes)
{}

Reco::ContainerVolume::~ContainerVolume()
{
    delete m_volumes;
}
