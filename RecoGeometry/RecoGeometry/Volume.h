//
//  Volume.h
//  
//
//  Created by Julia Hrdinka on 09/01/15.
//
//

#ifndef RECO_VOLUME_H
#define RECO_VOLUME_H

#include "Algebra/AlgPrimitives.h"
#include "TrkGeometryUtils/BinnedArray.h"
#include "RecoGeometry/BoundarySurface.h"
#include "RecoGeometry/Surface.h"
#include "RecoGeometry/Layer.h"
#include "TGeoNode.h"

namespace Reco {
    
    class Volume {
        
   //     class Layer;
        
    public:
        
        enum VolumeType {
            nEndCap   = 0,
            barrel    = 1,
            pEndCap   = 2,
            container = 3,
            none      = 4
        };
        
        typedef Trk::BinnedArray<Layer> LayerArray;
//        typedef std::vector<std::pair<std::shared_ptr<const Layer>,Alg::Point3D>> LayerArray;
        typedef std::vector<std::shared_ptr<const BoundarySurface>> SurfaceVector;
        
        virtual ~Volume();
        const Alg::Point3D& center() const;
        const Alg::Transform3D& transform() const;
  //      /afs/cern.ch/user/j/jhrdinka/FCC/FCCSW/RecoGeometry/RecoGeometry/Volume.hconst LayerArray& layers() const;
        virtual bool isInside(const Alg::Point3D& glopos, double tol=0.) const = 0;
        std::vector<const Layer*> materialLayersOrdered(const Alg::Point3D& glopos,const Alg::Vector3D& mom, double charge) const;
 //       SurfaceVector& boundarySurfaces();
        int NumberOfSurfaces() const;
        const BoundarySurface* getBoundarySurface (int n) const;
        //ueberlegen wie schoener moeglich
        void setType(VolumeType volumeType) const;
        VolumeType type() const;
        double getCoordinate(int n) const;
        
    
    protected:
        Volume();
        Volume(std::shared_ptr<const Alg::Transform3D> transf);
        Volume(LayerArray* layers);
        Volume(std::shared_ptr<const Alg::Transform3D> transf, LayerArray* layers);
        Volume(LayerArray* layers, TGeoNode* node);
        Volume(TGeoNode* node);
        
        mutable Alg::Point3D*                         m_center;
        std::shared_ptr<const Alg::Transform3D>       m_transform;
  //      const LayerArray*     m_layers;
        LayerArray*             m_layers;
        
        static Alg::Point3D     s_origin;
        static Alg::Transform3D s_idTransform;
        SurfaceVector           m_boundarySurfaces;
        
        mutable VolumeType      m_volumeType;
        std::vector<double>     m_coordinates;

    };
}


#endif //RECO_VOLUME_H
