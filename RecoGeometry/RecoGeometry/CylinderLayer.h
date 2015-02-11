//
//  CylinderLayer.h
//  
//
//  Created by Julia Hrdinka on 26/11/14.
//
//

#ifndef RECO_CYLINDERLAYER_H
#define RECO_CYLINDERLAYER_H

#include "RecoGeometry/Layer.h"
#include "RecoGeometry/CylinderSurface.h"
#include "TrkGeometryUtils/BinnedArray.h"

namespace Reco {
    
    class CylinderLayer : public Layer, public CylinderSurface {
    
    public:
        
        CylinderLayer(TGeoNode* node, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf);
        CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf);
        CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double halfZ, Trk::BinnedArray<Surface>* sf);
        ~CylinderLayer();
        
        virtual const Surface* getModule(const Alg::Point3D& glopos) const override;
        //locpos on layer
        virtual const Surface* getModule(const Alg::Point2D& locpos) const override;
//        virtual bool insideLayer(const Alg::Point3D& glopos, double tol=0.) const override;
        double getRmin() const;
        double getRmax() const;
        double getDR()   const;
        virtual void setNextLayer(std::shared_ptr<const Layer> layer) const override;
        virtual void setPreviousLayer(std::shared_ptr<const Layer> layer) const override;
        virtual const Layer* getNextLayer() const override;
        virtual const Layer* getPreviousLayer() const override;
        virtual const Layer* getNextLayer(const Alg::Vector3D& dir) const override;
        virtual const Layer* getPreviousLayer(const Alg::Vector3D& dir) const override;
        virtual bool onLayer(const Alg::Point3D& glopos) const override;
        virtual LayerType type() const override;
        virtual const CylinderSurface* surfaceRepresentation() const override;
        
    private:
        
        CylinderLayer();
        const double m_Rmin;
        const double m_Rmax;
    };
    
}


#endif //RECO_CYLINDERLAYER_H
