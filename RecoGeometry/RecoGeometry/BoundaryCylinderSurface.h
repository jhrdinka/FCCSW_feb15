//
//  BoundaryCylinderSurface.h
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#ifndef RECO_BOUNDARYCYLINDERSURFACE_H
#define RECO_BOUNDARYCYLINDERSURFACE_H

#include "RecoGeometry/CylinderSurface.h"
#include "RecoGeometry/BoundarySurface.h"

namespace Reco {

    class BoundaryCylinderSurface : public BoundarySurface, public CylinderSurface {
    
    public:
        
        BoundaryCylinderSurface(TGeoNode* node, TGeoConeSeg* tube);
        BoundaryCylinderSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf);
        BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D> transf, double radius, double halfZ);
        ~BoundaryCylinderSurface();
        
        BoundaryCylinderSurface& operator=(const BoundaryCylinderSurface& boundarycylindersurface);
        
 //       virtual void setNextVolume(Volume* nextVolume) const override;
 //       virtual void setPreviousVolume(Volume* previousVolume) const override;
 //       virtual void setNextVolumes(VolumeArray* nextVolumes) const override;
 //       virtual void setPreviousVolumes(VolumeArray* previousVolumes) const override;
//        virtual const Volume* getNextVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const;
//        virtual const Volume* getPreviousVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const;
        virtual const CylinderSurface* surfaceRepresentation() const;
        
    };
}


#endif //RECO_BOUNDARYCYLINDERSURFACE_H
