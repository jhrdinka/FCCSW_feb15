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
{
    //inner Cylinder
 /*   if (m_Rmin>0.) {
        std::shared_ptr<Reco::CylinderSurface> innerCylinder (new Reco::CylinderSurface(new Alg::Transform3D(transform()),m_Rmin,getHalfZ()));
        m_boundaries->push_back(innerCylinder);
    }
    //outer Cylinder
    std::shared_ptr<Reco::CylinderSurface> outerCylinder (new Reco::CylinderSurface(new Alg::Transform3D(transform()),m_Rmax,getHalfZ()));
    m_boundaries->push_back(outerCylinder);*/
}

Reco::CylinderLayer::CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
Layer(sf),
CylinderSurface(tube, transf),
m_Rmin(tube->GetRmin1()),
m_Rmax(tube->GetRmax1())
{
/*    //inner Cylinder
    if (m_Rmin>0.) {
        std::shared_ptr<Reco::CylinderSurface> innerCylinder (new Reco::CylinderSurface(transf,m_Rmin,getHalfZ()));
        m_boundaries->push_back(innerCylinder);
    }
    //outer Cylinder
    std::shared_ptr<Reco::CylinderSurface> outerCylinder (new Reco::CylinderSurface(transf,m_Rmax,getHalfZ()));
    m_boundaries->push_back(outerCylinder); */
}

Reco::CylinderLayer::CylinderLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double halfZ, Trk::BinnedArray<Surface>* sf) :
Layer(sf),
CylinderSurface(transf, 0.5*(rmin+rmax), halfZ),
m_Rmin(rmin),
m_Rmax(rmax)
{
/*    //inner Cylinder
    if (m_Rmin>0.) {
        std::shared_ptr<Reco::CylinderSurface> innerCylinder (new Reco::CylinderSurface(transf,m_Rmin,halfZ));
        m_boundaries->push_back(innerCylinder);
    }
    //outer Cylinder
    std::shared_ptr<Reco::CylinderSurface> outerCylinder (new Reco::CylinderSurface(transf,m_Rmax,halfZ));
    m_boundaries->push_back(outerCylinder); */
}

Reco::CylinderLayer::~CylinderLayer()
{}

const Reco::Surface* Reco::CylinderLayer::getModule(const Alg::Point3D& glopos) const
{
    Alg::Point2D* locpos = new Alg::Point2D();
    const Alg::Vector3D* mom = new Alg::Vector3D();
    if (this->globalToLocal(glopos, *mom, *locpos)) return (getModule(*locpos));
    return (0);
}

const Reco::Surface* Reco::CylinderLayer::getModule(const Alg::Point2D& locpos) const
{
    return (m_surfaces->object(locpos));
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
    return getNextLayer(); //noch ueberlegen wie das gehen soll
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

/*bool Reco::CylinderLayer::insideLayer(const Alg::Point3D& glopos, double tol) const
{
    double radius = sqrt(glopos.X()*glopos.X()+glopos.Y()*glopos.Y());
    return ((radius < (m_Rmax+tol)) && (radius > (m_Rmin-tol)) && (fabs(glopos.Z()) < (this->getHalfZ()+tol)));
}*/
const Reco::CylinderSurface* Reco::CylinderLayer::surfaceRepresentation() const
{
    return (this);
}





