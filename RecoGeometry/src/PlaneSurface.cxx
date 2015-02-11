//
//  PlaneSurface.cxx
//  
//
//  Created by Julia Hrdinka on 25/11/14.
//
//

#include "RecoGeometry/PlaneSurface.h"

Reco::PlaneSurface::PlaneSurface() :
Reco::Surface(),
m_halfX(0.),
m_halfY(0.)
{}

Reco::PlaneSurface::PlaneSurface(std::shared_ptr<const Alg::Transform3D> transf, double halfX, double halfY) :
Reco::Surface (transf),
m_halfX(halfX),
m_halfY(halfY)
{}

Reco::PlaneSurface::PlaneSurface(TGeoNode* node, TGeoBBox* box) :
Reco::Surface (node)
{
    m_halfX = box->GetDX();
    m_halfY = box->GetDY();
}

Reco::PlaneSurface::PlaneSurface(TGeoBBox* box, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface (transf)
{
    m_halfX = box->GetDX();
    m_halfY = box->GetDY();
}

Reco::PlaneSurface::PlaneSurface(TGeoNode* node,  TGeoBBox* box, Reco::MaterialMap* materialmap) :
Reco::Surface (node, materialmap)
{
    m_halfX = box->GetDX();
    m_halfY = box->GetDY();
}

Reco::PlaneSurface::PlaneSurface(TGeoBBox* box, Reco::MaterialMap* materialmap, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface (materialmap, transf)
{
    m_halfX = box->GetDX();
    m_halfY = box->GetDY();
}

Reco::PlaneSurface::~PlaneSurface()
{}

double Reco::PlaneSurface::halflengthX() const
{
    return (m_halfX);
}

double Reco::PlaneSurface::halflengthY() const
{
    return (m_halfY);
}

const Alg::Vector3D& Reco::PlaneSurface::normal() const
{
    return (Reco::Surface::normal());
}

const Alg::Vector3D* Reco::PlaneSurface::normal(const Alg::Point2D&) const
{
    return (new Alg::Vector3D(normal()));
}

bool Reco::PlaneSurface::isInside(const Alg::Point2D& locpos, double tol1, double tol2) const
{
    return ((fabs(locpos.X()) < (m_halfX + tol1)) && (fabs(locpos.Y()) < (m_halfY + tol2)));
}

void Reco::PlaneSurface::localToGlobal(const Alg::Point2D& locpos, const Alg::Vector3D&, Alg::Point3D& glopos) const
{
    Alg::Point3D loc3D (locpos.X(),locpos.Y(),0.);
    glopos = (transform())*loc3D;
}

bool Reco::PlaneSurface::globalToLocal(const Alg::Point3D& glopos, const Alg::Vector3D&, Alg::Point2D& locpos) const
{
    Alg::Point3D loc3D (transform().Inverse()*glopos);
    locpos.SetCoordinates(loc3D.X(),loc3D.Y());
    return(isInside(locpos,s_onSurfaceTolerance,s_onSurfaceTolerance) && loc3D.Z()*loc3D.Z()<s_onSurfaceTolerance*s_onSurfaceTolerance);
}
