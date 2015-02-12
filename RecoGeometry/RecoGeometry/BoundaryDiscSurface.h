//
//  BoundaryDiscSurface.h
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#ifndef RECO_BOUNDARYDISCSURFACE_H
#define RECO_BOUNDARYDISCSURFACE_H

#include "RecoGeometry/DiscSurface.h"
#include "RecoGeometry/BoundarySurface.h"

namespace Reco {
    
    class BoundaryDiscSurface : public BoundarySurface, public DiscSurface {
        
    public:
        
        BoundaryDiscSurface(TGeoNode* node, TGeoConeSeg* tube);
        BoundaryDiscSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf);
        BoundaryDiscSurface(std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax);
   //     BoundaryDiscSurface(Volume* nextVolume, Volume* previousVolume);
   //     BoundaryDiscSurface(VolumeArray* nextVolumes, VolumeArray* previousVolumes);
        ~BoundaryDiscSurface();
        
        BoundaryDiscSurface& operator=(const BoundaryDiscSurface& boundarydiscsurface);
        
//        virtual void setNextVolume(Volume* nextVolume) const override;
//        virtual void setPreviousVolume(Volume* previousVolume) const override;
//        virtual void setNextVolumes(VolumeArray* nextVolumes) const override;
//        virtual void setPreviousVolumes(VolumeArray* previousVolumes) const override;
//        virtual const Volume* getNextVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const;
//        virtual const Volume* getPreviousVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const;
        virtual const DiscSurface* surfaceRepresentation() const;
    
    };
}
#endif //RECO_BOUNDARYDISCSURFACE_H
