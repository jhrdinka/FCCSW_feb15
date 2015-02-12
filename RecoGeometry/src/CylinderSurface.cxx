//
//  CylinderSurface.cxx
//  
//
//  Created by Julia Hrdinka on 25/11/14.
//
//

#include "RecoGeometry/CylinderSurface.h"
#include "Algebra/RealQuadraticEquation.h"

Reco::CylinderSurface::CylinderSurface() :
Reco::Surface(),
m_R(0.),
m_halfZ(0.)
{}

Reco::CylinderSurface::CylinderSurface(std::shared_ptr<const Alg::Transform3D> transf, double radius, double halfZ) :
Reco::Surface(transf),
m_R(radius),
m_halfZ(halfZ)
{}

Reco::CylinderSurface::CylinderSurface(double radius, double halfZ) :
Reco::Surface(),
m_R(radius),
m_halfZ(halfZ)
{}

Reco::CylinderSurface::CylinderSurface(TGeoNode* node, TGeoConeSeg* tube) :
Reco::Surface(node)
{
    m_R = 0.5*(tube->GetRmax1()+tube->GetRmin1());
    //now is cone -> average from rmin1, rmin2, rmax1, rmax2??? brauch ich ned, weil myconeseg nur mit einem radius def is
    //maybe if (rmax1==rmax2...)
    m_halfZ = tube->GetDz();
}

Reco::CylinderSurface::CylinderSurface(TGeoConeSeg* tube, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface(transf)
{
    m_R = 0.5*(tube->GetRmax1()+tube->GetRmin1());
    //now is cone -> average from rmin1, rmin2, rmax1, rmax2??? brauch ich ned, weil myconeseg nur mit einem radius def is
    //maybe if (rmax1==rmax2...)
    m_halfZ = tube->GetDz();
}

Reco::CylinderSurface::CylinderSurface(TGeoNode* node, TGeoConeSeg* tube, Reco::MaterialMap* materialmap) :
Reco::Surface(node,materialmap)
{
    m_R = 0.5*(tube->GetRmax1()+tube->GetRmin1());
    //now is cone -> average from rmin1, rmin2, rmax1, rmax2???
    //maybe if (rmax1==rmax2...)
    m_halfZ = tube->GetDz();
}

Reco::CylinderSurface::CylinderSurface(TGeoConeSeg* tube, Reco::MaterialMap* materialmap, std::shared_ptr<const Alg::Transform3D> transf) :
Reco::Surface(materialmap, transf)
{
    m_R = 0.5*(tube->GetRmax1()+tube->GetRmin1());
    //now is cone -> average from rmin1, rmin2, rmax1, rmax2???
    //maybe if (rmax1==rmax2...)
    m_halfZ = tube->GetDz();
}

Reco::CylinderSurface::~CylinderSurface()
{}

Reco::CylinderSurface& Reco::CylinderSurface::operator=(const CylinderSurface& cylindersurface)
{
    if (this!=&cylindersurface) {
        Reco::Surface::operator=(cylindersurface);
        m_R = cylindersurface.m_R;
        m_halfZ = cylindersurface.m_halfZ;
    }
    return (*this);
}

double Reco::CylinderSurface::getR() const
{
    return (m_R);
}

double Reco::CylinderSurface::getHalfZ() const
{
    return (m_halfZ);
}

const Alg::Vector3D& Reco::CylinderSurface::normal() const
{
    return (Reco::Surface::normal());
}

const Alg::Vector3D* Reco::CylinderSurface::normal(const Alg::Point2D& locpos) const
{
    double phi      = locpos.X()/m_R;
    Alg::Vector3D normal(m_R*cos(phi), m_R*sin(phi), locpos.Y());
    normal          = transform()*normal;
    
    return (new Alg::Vector3D(normal.Unit()));
}

bool Reco::CylinderSurface::isInside(const Alg::Point2D& locpos, double tol1, double tol2) const
{
    return ((fabs(locpos.Y()) < (m_halfZ + tol2)) && (fabs(locpos.X()) < (M_PI*m_R + tol1)));
}

void Reco::CylinderSurface::localToGlobal(const Alg::Point2D& locpos, const Alg::Vector3D&, Alg::Point3D& glopos) const
{
    double phi      = locpos.X()/m_R;
    double x        = m_R*cos(phi);
    double y        = m_R*sin(phi);
    Alg::Point3D loc3D(x,y,locpos.Y());
    glopos          = transform()*loc3D;
}

bool Reco::CylinderSurface::globalToLocal(const Alg::Point3D& glopos, const Alg::Vector3D&, Alg::Point2D& locpos) const
{
    Alg::Point3D loc3D(transform().Inverse()*glopos);
    double r        = sqrt(loc3D.X()*loc3D.X()+loc3D.Y()*loc3D.Y());
    double phi      = atan(loc3D.Y()/loc3D.X());
    locpos.SetCoordinates(r*phi,loc3D.Z());
    
    return (isInside(locpos,s_onSurfaceTolerance,s_onSurfaceTolerance));
}

Trk::Intersection Reco::CylinderSurface::straightLineIntersection(const Alg::Point3D& pos,
 const Alg::Vector3D& dir,
 bool forceDir) const
{
    Alg::Point3D point1    = transform().Inverse()*pos;
    Alg::Vector3D direction = transform().Inverse().Rotation()*dir;
    Alg::Point3D point2    = point1 + direction;
    double R  = getR();
    double t1 = 0.;
    double t2 = 0.;
    
     if (direction.x())
     {
         // get line and circle constants
         double k = (direction.Y())/(direction.X());
         double d = (point2.X()*point1.Y() - point1.X()*point2.Y())/(point2.X()-point1.X());
    
         // and solve the qaudratic equation
         Alg::RealQuadraticEquation pquad(1 + k*k, 2*k*d, d*d - R*R);
         if (pquad.solutions != Alg::none)
         {
             // the solutions in the 3D frame of the cylinder
             t1 = (pquad.first  - point1.X())/direction.X();
             t2 = (pquad.second - point1.X())/direction.X();
         }
         else  // bail out if no solution exists
             return Trk::Intersection(pos, 0., false);
     }
     else
     {
         // x value is the one of point1
         // x^2 + y^2 = R^2
         // y = sqrt(R^2-x^2)
         double x  = point1.X();
         double r2mx2 = R*R-x*x;
         // bail out if no solution
         if (r2mx2 < 0. ) return Trk::Intersection(pos, 0., false);
         double y  = sqrt(r2mx2);
         // assign parameters and solutions
         t1 = y-point1.y();
         t2 = -y-point1.y();
     }
    
    Alg::Point3D sol1raw(point1 + t1 * direction);
    Alg::Point3D sol2raw(point1 + t2 * direction);
    // now reorder and return
    Alg::Point3D solution(0,0,0);
    double path = 0.;
    
    // first check the validity of the direction
    bool isValid = true;
    
    // both solutions are of same sign, take the smaller, but flag as false if not forward
    if (t1*t2 > 0 || !forceDir)
    {
        // asign validity
        isValid = forceDir ? ( t1 > 0. ) : true ;
        // assign the right solution
        if (t1*t1 < t2*t2)
        {
            solution = sol1raw;
            path = t1;
        }
        else
        {
            solution = sol2raw;
            path = t2;
        }
    }
    else
    {
        if (t1 > 0.) {
            solution = sol1raw;
            path = t1;
        }
        else
        {
            solution = sol2raw;
            path = t2;
        }
    }
    // the solution is still in the local 3D frame, direct check
    Alg::Point2D loc;
    if (isValid) isValid = (globalToLocal(pos,dir,loc) && isInside(loc,s_onSurfaceTolerance,s_onSurfaceTolerance));
    
    // now return
    return Trk::Intersection(transform()*solution, path, isValid );
}


