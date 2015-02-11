//
//  NavigationLayer.h
//  
//
//  Created by Julia Hrdinka on 14/01/15.
//
//

#ifndef RECO_NAVIGATIONLAYER_H
#define RECO_NAVIGATIONLAYER_H

#include "RecoGeometry/Layer.h"

namespace Reco {
    
    class NavigationLayer : public Layer {
    
    public:
        
        NavigationLayer();
        NavigationLayer(std::shared_ptr<const Layer> next,std::shared_ptr<const Layer> previous);
        ~NavigationLayer();
        virtual void setNextLayer(std::shared_ptr<const Layer> layer) const override;
        virtual void setPreviousLayer(std::shared_ptr<const Layer> layer) const override;
        virtual const Layer* getNextLayer() const override;
        virtual const Layer* getPreviousLayer() const override;
        virtual const Layer* getNextLayer(const Alg::Vector3D& dir) const override;
        virtual const Layer* getPreviousLayer(const Alg::Vector3D& dir) const override;
        //only because of inheritance
        virtual LayerType type() const override;
        //gibt false zurueck
        virtual bool onLayer(const Alg::Point3D& glopos) const override;
        virtual const Surface* surfaceRepresentation() const override;
    };
}

#endif //RECO_NAVIGATIONLAYER_H
