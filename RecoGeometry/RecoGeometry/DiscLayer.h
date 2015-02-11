//
//  DiscLayer.h
//  
//
//  Created by Julia Hrdinka on 08/12/14.
//
//

#ifndef RECO_DISCLAYER_H
#define RECO_DISCLAYER_H

#include "RecoGeometry/Layer.h"
#include "RecoGeometry/DiscSurface.h"
#include "TrkGeometryUtils/BinnedArray.h"

namespace Reco  {
    
    class DiscLayer : public DiscSurface, public Layer {
    
    public:
        
        DiscLayer(TGeoNode* node, TGeoConeSeg* disc, Trk::BinnedArray<Surface>* sf);
        DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* disc, Trk::BinnedArray<Surface>* sf);
        DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double dz, Trk::BinnedArray<Surface>* sf);
        ~DiscLayer();
        
        virtual const Surface* getModule(const Alg::Point3D& glopos) const override;
        //locpos on layer
        virtual const Surface* getModule(const Alg::Point2D& locpos) const override;
   //     virtual bool insideLayer(const Alg::Point3D& glopos, double tol=0.) const override;
        double getHalfZ() const;
        virtual void setNextLayer(std::shared_ptr<const Layer> layer) const override;
        virtual void setPreviousLayer(std::shared_ptr<const Layer> layer) const override;
        virtual const Layer* getNextLayer() const override;
        virtual const Layer* getPreviousLayer() const override;
        virtual const Layer* getNextLayer(const Alg::Vector3D& dir) const override;
        virtual const Layer* getPreviousLayer(const Alg::Vector3D& dir) const override;
        virtual bool onLayer(const Alg::Point3D& glopos) const override;
        virtual LayerType type() const override;
        virtual const DiscSurface* surfaceRepresentation() const override;
        
    private:
        
        DiscLayer();
        const double m_dz;
    };
}

#endif //RECO_DISCLAYER_H
