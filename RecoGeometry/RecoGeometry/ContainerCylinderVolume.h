//
//  ContainerCylinderVolume.h
//  
//
//  Created by Julia Hrdinka on 27/01/15.
//
//

#ifndef RECO_CONTAINERCYLINDERVOLUME_H
#define RECO_CONTAINERCYLINDERVOLUME_H

#include "RecoGeometry/ContainerVolume.h"
#include "RecoGeometry/CylinderVolume.h"

namespace Reco {
    
    class ContainerCylinderVolume : public ContainerVolume, public CylinderVolume {
        
    public:
        
        enum binZType {
            nEnd      = 0,
            nBarrel   = 1,
            pBarrel   = 2,
            pEnd      = 3,
            zLength   = 4
        };
        
        enum binRType {
            inner   = 0,
            middle  = 1,
            outer   = 2,
            rLength = 3
        };
        
        ContainerCylinderVolume(VolumeArray* volumes, std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax, double halfZ);
        ContainerCylinderVolume(VolumeArray* volumes, TGeoNode* node, TGeoConeSeg* tube);
        virtual bool isInside(const Alg::Point3D& glopos, double tol=0.) const;
        ~ContainerCylinderVolume();
    
    };
}

#endif //RECO_CONTAINERCYLINDERVOLUME_H
