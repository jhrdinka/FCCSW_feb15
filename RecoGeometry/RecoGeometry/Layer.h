//
//  Layer.h
//  
//
//  Created by Julia Hrdinka on 26/11/14.
//
//

#ifndef RECO_LAYER_H
#define RECO_LAYER_H

#include "TrkGeometryUtils/BinnedArray.h"

namespace Reco {
    
    class Surface;
    
    typedef Trk::BinnedArray<Surface> SurfaceArray;
    typedef std::vector<std::shared_ptr<const Surface>> SurfaceVector;
    
    class Layer {
    
    public:
        
        enum LayerType {
            navigation  = 0,
            cylinder    = 1,
            disc        = 2,
        };

 //       Layer (SurfaceArray* sf, Reco::Layer* nextLayer, Reco::Layer* previousLayer);
        Layer (SurfaceArray* sf);
        Layer (const Layer& layer);
        virtual ~Layer();

        const SurfaceArray* getSurfaces() const;
        virtual const Surface* getModule(const Alg::Point3D& glopos) const;
        //locpos on layer
        virtual const Surface* getModule(const Alg::Point2D& locpos) const;
 //       virtual bool insideLayer(const Alg::Point3D& glopos, double tol=0.) const = 0;
        virtual void setNextLayer(std::shared_ptr<const Layer> layer) const;
        virtual void setPreviousLayer(std::shared_ptr<const Layer> layer) const;
        virtual const Layer* getNextLayer() const;
        virtual const Layer* getPreviousLayer() const;
        virtual const Layer* getNextLayer(const Alg::Vector3D& dir) const;
        virtual const Layer* getPreviousLayer(const Alg::Vector3D& dir) const;
        virtual bool onLayer(const Alg::Point3D& glopos) const = 0;
        virtual LayerType type() const = 0;
        virtual const Surface* surfaceRepresentation() const = 0;
        
        //Methoden mit nextLayer(direction) und previousLayer(direction) spaeter
        
    protected:
        Layer ();
        
        const SurfaceArray*                     m_surfaces;
  //      SurfaceVector*                          m_boundaries;
        mutable std::shared_ptr<const Layer>    m_nextLayer;
        mutable std::shared_ptr<const Layer>    m_previousLayer;
    };
}

#endif //RECO_LAYER_H
