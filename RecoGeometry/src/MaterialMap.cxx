//
//  MaterialMap.cxx
//  
//
//  Created by Julia Hrdinka on 05/01/15.
//
//

#include "RecoGeometry/MaterialMap.h"

Reco::MaterialMap::MaterialMap(Trk::BinUtility* binutility, std::map<std::pair<int,int>,Reco::Material>* material) :
    m_binutility (binutility),
    m_material (material)
{}

Reco::MaterialMap::~MaterialMap()
{
    m_binutility  = 0;
    m_material    = 0;
}

size_t Reco::MaterialMap::bin(const Alg::Point2D& lposition, size_t ba) const
{
    return (m_binutility->bin(lposition, ba));
}

const Reco::Material& Reco::MaterialMap::material(Alg::Point2D& locpos) const
{
    int binx = bin(locpos,0);
    int biny = bin(locpos,1);
    //schauen, ob eh passt
    return (m_material->at(std::make_pair(binx,biny)));
}