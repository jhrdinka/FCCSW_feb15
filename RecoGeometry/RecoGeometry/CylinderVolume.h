//
//  CylinderVolume.h
//  
//
//  Created by Julia Hrdinka on 09/01/15.
//
//

#ifndef RECO_CYLINDERVOLUME_H
#define RECO_CYLINDERVOLUME_H

#include "RecoGeometry/Volume.h"
#include "TGeoCone.h"

namespace Reco {
    
    class CylinderVolume : public virtual Volume {
        
    public:
        
        enum BoundarySurfaceType {
            negDisc       = 0,
            posDisc       = 1,
            outerCylinder = 2,
            innerCylinder = 3,
            length        = 4
        };
        
        enum Coordinates {
            Rmin    = 0,
            Rmax    = 1,
            halfZ   = 2
        };
        
        CylinderVolume(std::shared_ptr<const Alg::Transform3D> transf, LayerArray* layers, double Rmin, double Rmax, double halfZ);
   //     CylinderVolume(double Rmin, double Rmax, double halfZ);
        CylinderVolume(std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax, double halfZ);
   //     CylinderVolume(LayerArray* layers, double Rmin, double Rmax, double halfZ);
        CylinderVolume(LayerArray* layers, TGeoNode* node, TGeoConeSeg* tube);
        CylinderVolume(TGeoNode* node, TGeoConeSeg* tube);
        virtual ~CylinderVolume();
        
        double getRmin() const;
        double getRmax() const;
        double getHalfZ() const;
        const BoundarySurface* getPosDisc() const;
        const BoundarySurface* getNegDisc() const;
        const BoundarySurface* getOuterCylinder() const;
        const BoundarySurface* getInnerCylinder() const;
        
        
        virtual bool isInside(const Alg::Point3D& glopos, double tol=0.) const override;
 //       virtual std::vector<Layer*> materialLayersOrdered(const Alg::Point3D& position,const Alg::Vector3D& momentum, double charge) const override;
        
    private:
        double          m_Rmin;
        double          m_Rmax;
        double          m_halfZ;
  //      SurfaceVector   m_boundarySurfaces;
    
    };
}

#endif //RECO_CYLINDERVOLUME_H
