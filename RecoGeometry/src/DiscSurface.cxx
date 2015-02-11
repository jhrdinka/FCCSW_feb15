//
//  DiscSurface.cxx
//  
//
//  Created by Julia Hrdinka on 02/12/14.
//
//

#include "RecoGeometry/DiscSurface.h"
#include "Algebra/RealQuadraticEquation.h"

Reco::DiscSurface::DiscSurface() :
Reco::Surface(),
m_Rmin(0.),
m_Rmax(0.)
{}

Reco::DiscSurface::DiscSurface(std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax) :
Reco::Surface(transf),
m_Rmin(Rmin),
m_Rmax(Rmax)
{}

Reco::DiscSurface::DiscSurface(TGeoNode* node, TGeoConeSeg* tube) :
Reco::Surface(node)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    //vl noch ueberpruefung mit z machen??

}

Reco::DiscSurface::DiscSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface(transf)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    //vl noch ueberpruefung mit z machen??
    
}

Reco::DiscSurface::DiscSurface(TGeoNode* node, TGeoConeSeg* tube, Reco::MaterialMap* materialmap) :
Reco::Surface(node, materialmap)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    //vl noch ueberpruefung mit z machen??
}

Reco::DiscSurface::DiscSurface(TGeoConeSeg* tube, Reco::MaterialMap* materialmap, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface(materialmap, transf)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    //vl noch ueberpruefung mit z machen??
}

Reco::DiscSurface::~DiscSurface()
{}

double Reco::DiscSurface::getRmin() const
{
    return (m_Rmin);
}

double Reco::DiscSurface::getRmax() const
{
    return (m_Rmax);
}

const Alg::Vector3D& Reco::DiscSurface::normal() const
{
    return (Reco::Surface::normal());
}

const Alg::Vector3D* Reco::DiscSurface::normal(const Alg::Point2D&) const
{
    return (new Alg::Vector3D(normal()));
}

bool Reco::DiscSurface::isInside(const Alg::Point2D& locpos, double tol1, double tol2) const
{
    return ((fabs(locpos.Y()) <= (M_PI + tol2)) && (locpos.X() > (m_Rmin - tol1)) && (locpos.X() < (m_Rmax + tol1)));
}

void Reco::DiscSurface::localToGlobal(const Alg::Point2D& locpos, const Alg::Vector3D&, Alg::Point3D& glopos) const
{
    double x        = locpos.X()*cos(locpos.Y());
    double y        = locpos.X()*sin(locpos.Y());
    Alg::Point3D loc3D(x,y,0.);
    glopos          = transform()*loc3D;
}

bool Reco::DiscSurface::globalToLocal(const Alg::Point3D& glopos, const Alg::Vector3D&, Alg::Point2D& locpos) const
{
    Alg::Point3D loc3D(transform().Inverse()*glopos);
    double r        = sqrt(loc3D.X()*loc3D.X()+loc3D.Y()*loc3D.Y());
    double phi      = atan(loc3D.Y()/loc3D.X());
    locpos.SetCoordinates(r,phi);
    
    return (isInside(locpos,s_onSurfaceTolerance,s_onSurfaceTolerance) && loc3D.Z()*loc3D.Z()<s_onSurfaceTolerance*s_onSurfaceTolerance);
}