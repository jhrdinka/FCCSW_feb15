//
//  MaterialMap.h
//  
//
//  Created by Julia Hrdinka on 05/01/15.
//
//

#ifndef RECO_MATERIALMAP_H
#define RECO_MATERIALMAP_H

#include "TrkGeometryUtils/BinUtility.h"
#include "RecoGeometry/Material.h"

namespace Reco {
    
    class MaterialMap {
    
    public:
        MaterialMap (Trk::BinUtility* binutility, std::map<std::pair<int,int>, Material>* material);
        ~MaterialMap();
        size_t bin(const Alg::Point2D& lposition, size_t ba=0) const;
        const Material& material(Alg::Point2D& locpos) const;
        
    private:
        const Trk::BinUtility* m_binutility;
        const std::map<std::pair<int,int>, Material>* m_material;
    };
}
#endif //RECO_MATERIALMAP_H
