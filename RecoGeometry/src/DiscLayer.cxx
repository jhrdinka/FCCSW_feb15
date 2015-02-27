//
//  DiscLayer.cxx
//  
//
//  Created by Julia Hrdinka on 08/12/14.
//
//

#include "RecoGeometry/DiscLayer.h"

Reco::DiscLayer::DiscLayer() :
DiscSurface(),
Layer (),
m_dz (0.)
{}

Reco::DiscLayer::DiscLayer(TGeoNode* node, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
DiscSurface(node, tube),
Layer(sf),
m_dz(tube->GetDz())
{}

Reco::DiscLayer::DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
DiscSurface(tube,transf),
Layer(sf),
m_dz(tube->GetDZ())
{}

Reco::DiscLayer::DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double dz, Trk::BinnedArray<Surface>* sf) :
DiscSurface(transf, rmin, rmax),
Layer(sf),
m_dz(dz)
{}

Reco::DiscLayer::~DiscLayer()
{}

const Reco::Surface* Reco::DiscLayer::getModule(const Alg::Point3D& glopos) const
{
    if (this->onLayer(glopos)) {
        Alg::Point2D locpos(0.,0.);
        const Alg::Vector3D mom(0.,0.,0.);
        if (m_surfaces->object(locpos) && m_surfaces->object(locpos)->globalToLocal(glopos,mom,locpos)) return m_surfaces->object(locpos);
    }
    return (0);
}

const Reco::Surface* Reco::DiscLayer::getModule(const Alg::Point2D& locpos) const
{
    Alg::Point3D glopos(0.,0.,0.);
    const Alg::Vector3D mom(0.,0.,0.);
    this->localToGlobal(locpos,mom,glopos);
    return (Reco::DiscLayer::getModule(glopos));
}

double Reco::DiscLayer::getHalfZ() const
{
    return(m_dz);
}

void Reco::DiscLayer::setNextLayer(std::shared_ptr<const Layer> layer) const
{
    Reco::Layer::setNextLayer(layer);
}

void Reco::DiscLayer::setPreviousLayer(std::shared_ptr<const Layer> layer) const
{
    Reco::Layer::setPreviousLayer(layer);
}

const Reco::Layer* Reco::DiscLayer::getNextLayer() const
{
    return Reco::Layer::getNextLayer();
}

const Reco::Layer* Reco::DiscLayer::getPreviousLayer() const
{
    return Reco::Layer::getPreviousLayer();
}

const Reco::Layer* Reco::DiscLayer::getNextLayer(const Alg::Vector3D& dir) const
{
    if (dir.Z()*normal().Z()>0) return getNextLayer();
    if (dir.Z()*normal().Z()<0) return getPreviousLayer();
    else return 0;
    //noch ueberlegen wie das genau gehen soll eventuell mit punkt
    //was wenn
}

const Reco::Layer* Reco::DiscLayer::getPreviousLayer(const Alg::Vector3D& dir) const
{
    if (dir.Z()*normal().Z()>0) return getPreviousLayer();
    if (dir.Z()*normal().Z()<0) return getNextLayer();
    else return 0;
}

bool Reco::DiscLayer::onLayer(const Alg::Point3D& glopos) const
{
    Alg::Point3D locpos (transform().Inverse()*glopos);
    double r    = sqrt(locpos.X()*locpos.X()+locpos.Y()*locpos.Y());
    return ((r>=getRmin()) && (r<=getRmax()) && (fabs(locpos.Z()) <= m_dz));
}

Reco::Layer::LayerType Reco::DiscLayer::type() const
{
    return (Reco::Layer::disc);
}

const Reco::DiscSurface* Reco::DiscLayer::surfaceRepresentation() const
{
    return (this);
}



