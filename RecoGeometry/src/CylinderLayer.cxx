//
//  CylinderLayer.cxx
//  
//
//  Created by Julia Hrdinka on 26/11/14.
//
//

#include "RecoGeometry/CylinderLayer.h"

Reco::CylinderLayer::CylinderLayer() :
Layer (),
CylinderSurface(),
m_Rmin(0.),
m_Rmax(0.)
{}

Reco::CylinderLayer::CylinderLayer(TGeoNode* node, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
Layer(sf),
CylinderSurface(node, tube),
m_Rmin(tube->GetRmin1()),
m_Rmax(tube->GetRmax1())
{}

Reco::CylinderLayer::CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
Layer(sf),
CylinderSurface(tube, transf),
m_Rmin(tube->GetRmin1()),
m_Rmax(tube->GetRmax1())
{}

Reco::CylinderLayer::CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double halfZ, Trk::BinnedArray<Surface>* sf) :
Layer(sf),
CylinderSurface(transf, 0.5*(rmin+rmax), halfZ),
m_Rmin(rmin),
m_Rmax(rmax)
{}

Reco::CylinderLayer::~CylinderLayer()
{}

const Reco::Surface* Reco::CylinderLayer::getModule(const Alg::Point3D& glopos) const
{
    if (this->onLayer(glopos)) {
    Alg::Point2D locpos(0.,0.);
    const Alg::Vector3D mom(0.,0.,0.);
    if (m_surfaces->object(locpos) && m_surfaces->object(locpos)->globalToLocal(glopos,mom,locpos)) return m_surfaces->object(locpos);
    }
    return (0);
}

const Reco::Surface* Reco::CylinderLayer::getModule(const Alg::Point2D& locpos) const
{
    Alg::Point3D glopos(0.,0.,0.);
    const Alg::Vector3D mom(0.,0.,0.);
    this->localToGlobal(locpos,mom,glopos);
    return (Reco::CylinderLayer::getModule(glopos));
}

double Reco::CylinderLayer::getRmin() const
{
    return (m_Rmin);
}

double Reco::CylinderLayer::getRmax() const
{
    return (m_Rmax);
}
double Reco::CylinderLayer::getDR() const
{
    return (fabs(m_Rmax-m_Rmin));
}

void Reco::CylinderLayer::setNextLayer(std::shared_ptr<const Layer> layer) const
{
    Reco::Layer::setNextLayer(layer);
}

void Reco::CylinderLayer::setPreviousLayer(std::shared_ptr<const Layer> layer) const
{
    Reco::Layer::setPreviousLayer(layer);
}

const Reco::Layer* Reco::CylinderLayer::getNextLayer() const
{
    return Reco::Layer::getNextLayer();
}

const Reco::Layer* Reco::CylinderLayer::getPreviousLayer() const
{
    return Reco::Layer::getPreviousLayer();
}

const Reco::Layer* Reco::CylinderLayer::getNextLayer(const Alg::Vector3D&) const
{
    return getNextLayer(); //noch ueberlegen wie das genau gehen soll eventuell mit punkt
}
const Reco::Layer* Reco::CylinderLayer::getPreviousLayer(const Alg::Vector3D&) const
{
    return getPreviousLayer();
}
bool Reco::CylinderLayer::onLayer(const Alg::Point3D& glopos) const
{
    Alg::Point3D locpos (transform().Inverse()*glopos);
    double r    = sqrt(locpos.X()*locpos.X()+locpos.Y()*locpos.Y());
    return ((r>=m_Rmin) && (r<=m_Rmax) && (fabs(locpos.Z()) <= getHalfZ()));
}

Reco::Layer::LayerType Reco::CylinderLayer::type() const
{
    return (Reco::Layer::cylinder); 
}

const Reco::CylinderSurface* Reco::CylinderLayer::surfaceRepresentation() const
{
    return (this);
}





