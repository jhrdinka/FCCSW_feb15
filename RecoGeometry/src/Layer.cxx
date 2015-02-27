//
//  Layer.cxx
//  
//
//  Created by Julia Hrdinka on 27/11/14.
//
//

#include "RecoGeometry/Layer.h"


Reco::Layer::Layer() :
m_surfaces(0),
m_nextLayer(0),
m_previousLayer(0)
{}

Reco::Layer::Layer(const Reco::Layer& layer) :
m_surfaces(layer.m_surfaces),
m_nextLayer(layer.m_nextLayer),
m_previousLayer(layer.m_previousLayer)
{}

Reco::Layer::Layer(Reco::SurfaceArray* sf) :
m_surfaces(sf),
m_nextLayer(0),
m_previousLayer(0)
{}

Reco::Layer::~Layer()
{
    delete m_surfaces;
}


const Reco::SurfaceArray* Reco::Layer::getSurfaces() const
{
    return (m_surfaces);
}

void Reco::Layer::setNextLayer(std::shared_ptr<const Layer> layer) const
{
    m_nextLayer = layer;
}

void Reco::Layer::setPreviousLayer(std::shared_ptr<const Layer> layer) const
{
    m_previousLayer = layer;
}

const Reco::Layer* Reco::Layer::getNextLayer() const
{
    return (m_nextLayer.get());
}

const Reco::Layer* Reco::Layer::getPreviousLayer() const
{
    return (m_previousLayer.get());
}

const Reco::Layer* Reco::Layer::getNextLayer(const Alg::Vector3D&) const
{
    return getNextLayer();
}

const Reco::Layer* Reco::Layer::getPreviousLayer(const Alg::Vector3D&) const
{
    return getPreviousLayer();
}