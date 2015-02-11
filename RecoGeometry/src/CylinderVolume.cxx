//
//  CylinderVolume.cxx
//  
//
//  Created by Julia Hrdinka on 09/01/15.
//
//

#include "RecoGeometry/CylinderVolume.h"
#include "RecoGeometry/BoundaryCylinderSurface.h"
#include "RecoGeometry/BoundaryDiscSurface.h"

Reco::CylinderVolume::CylinderVolume(std::shared_ptr<const Alg::Transform3D> transf, LayerArray* layers, double Rmin, double Rmax, double halfZ) :
Reco::Volume(transf, layers),
m_Rmin(Rmin),
m_Rmax(Rmax),
m_halfZ(halfZ)
{
    //Discs
    Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transf->Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transf->Translation().Vect();
    std::shared_ptr<const Alg::Transform3D> rightTransform(new Alg::Transform3D(rotation,center+lz*halfZ));
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    std::shared_ptr<const Alg::Transform3D> leftTransform(new Alg::Transform3D(leftrotation,center+lz*halfZ));
    
    std::shared_ptr<Reco::BoundaryDiscSurface> posDisc (new Reco::BoundaryDiscSurface(rightTransform,Rmin,Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> negDisc (new Reco::BoundaryDiscSurface(leftTransform,Rmin,Rmax));
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D>(transf),Rmax,halfZ));
    
    m_boundarySurfaces.push_back(negDisc);
    m_boundarySurfaces.push_back(posDisc);
    m_boundarySurfaces.push_back(outerCylinder);
    //inner Cylinder
    if (Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D>(transf),Rmin,halfZ));
        m_boundarySurfaces.push_back(innerCylinder);
    }
    //set Coordinates
    m_coordinates[Reco::CylinderVolume::Rmin]   = m_Rmin;
    m_coordinates[Reco::CylinderVolume::Rmax]   = m_Rmax;
    m_coordinates[Reco::CylinderVolume::halfZ]  = m_halfZ;
    
}

/*Reco::CylinderVolume::CylinderVolume(double Rmin, double Rmax, double halfZ) :
Reco::Volume(),
m_Rmin(Rmin),
m_Rmax(Rmax),
m_halfZ(halfZ)
{
    //Discs
    Alg::Vector3D zdir(0.,0.,1.);
    std::shared_ptr<Reco::BoundaryDiscSurface> rightDisc (new Reco::BoundaryDiscSurface(new Alg::Transform3D(zdir*halfZ),Rmin,Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> leftDisc (new Reco::BoundaryDiscSurface(new Alg::Transform3D(Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI),zdir*halfZ),Rmin,Rmax));
    rightDisc.get()->setType(Reco::Volume::posDisc);
    leftDisc.get()->setType(Reco::Volume::negDisc);
    
    m_surfaces->push_back(rightDisc);
    m_surfaces->push_back(leftDisc);
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(Rmax,halfZ));
    outerCylinder.get()->setType(Reco::Volume::outerCylinder);
    m_surfaces->push_back(outerCylinder);
    
    //inner Cylinder
    if (Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(Rmin,halfZ));
        innerCylinder.get()->setType(Reco::Volume::innerCylinder);
        m_surfaces->push_back(innerCylinder);
    }
}
*/

Reco::CylinderVolume::CylinderVolume(std::shared_ptr<const Alg::Transform3D> transf, double Rmin, double Rmax, double halfZ) :
Reco::Volume(transf),
m_Rmin(Rmin),
m_Rmax(Rmax),
m_halfZ(halfZ)
{
    //Discs
    Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transf->Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transf->Translation().Vect();
    std::shared_ptr<const Alg::Transform3D> rightTransform(new Alg::Transform3D(rotation,center+lz*halfZ));
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    std::shared_ptr<const Alg::Transform3D> leftTransform(new Alg::Transform3D(leftrotation,center+lz*halfZ));
    
    std::shared_ptr<Reco::BoundaryDiscSurface> posDisc (new Reco::BoundaryDiscSurface(rightTransform,Rmin,Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> negDisc (new Reco::BoundaryDiscSurface(leftTransform,Rmin,Rmax));
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D>(transf),Rmax,halfZ));
    
    m_boundarySurfaces.push_back(negDisc);
    m_boundarySurfaces.push_back(posDisc);
    m_boundarySurfaces.push_back(outerCylinder);
    
    //inner Cylinder
    if (Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(std::shared_ptr<const Alg::Transform3D>(transf),Rmin,halfZ));
        m_boundarySurfaces.push_back(innerCylinder);
    }
    //set Coordinates
    m_coordinates[Reco::CylinderVolume::Rmin]   = m_Rmin;
    m_coordinates[Reco::CylinderVolume::Rmax]   = m_Rmax;
    m_coordinates[Reco::CylinderVolume::halfZ]  = m_halfZ;
}

/*
Reco::CylinderVolume::CylinderVolume(LayerArray* layers, double Rmin, double Rmax, double halfZ) :
Reco::Volume(layers),
m_Rmin(Rmin),
m_Rmax(Rmax),
m_halfZ(halfZ)
{
    //Discs
    Alg::Vector3D zdir(0.,0.,1.);
    std::shared_ptr<Reco::BoundaryDiscSurface> rightDisc (new Reco::BoundaryDiscSurface(new Alg::Transform3D(zdir*halfZ),Rmin,Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> leftDisc (new Reco::BoundaryDiscSurface(new Alg::Transform3D(Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI),zdir*halfZ),Rmin,Rmax));
    rightDisc.get()->setType(Reco::Volume::posDisc);
    leftDisc.get()->setType(Reco::Volume::negDisc);
    
    m_surfaces->push_back(rightDisc);
    m_surfaces->push_back(leftDisc);
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(Rmax,halfZ));
    outerCylinder.get()->setType(Reco::Volume::outerCylinder);
    m_surfaces->push_back(outerCylinder);
    
    //inner Cylinder
    if (Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(Rmin,halfZ));
        innerCylinder.get()->setType(Reco::Volume::innerCylinder);
        m_surfaces->push_back(innerCylinder);
    }
}
*/
Reco::CylinderVolume::CylinderVolume(LayerArray* layers, TGeoNode* node, TGeoConeSeg* tube) :
Reco::Volume(layers,node)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    m_halfZ= tube->GetDZ();
    
    //Discs
    Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transform().Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transform().Translation().Vect();
    std::shared_ptr<const Alg::Transform3D> rightTransform(new Alg::Transform3D(rotation,center+lz*m_halfZ));
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    std::shared_ptr<const Alg::Transform3D> leftTransform(new Alg::Transform3D(leftrotation,center+lz*m_halfZ));
    
    std::shared_ptr<Reco::BoundaryDiscSurface> posDisc (new Reco::BoundaryDiscSurface(rightTransform,m_Rmin,m_Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> negDisc (new Reco::BoundaryDiscSurface(leftTransform,m_Rmin,m_Rmax));
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(std::make_shared<const Alg::Transform3D>(Alg::Transform3D(transform())),m_Rmax,m_halfZ));
    
    m_boundarySurfaces.push_back(negDisc);
    m_boundarySurfaces.push_back(posDisc);
    m_boundarySurfaces.push_back(outerCylinder);
    //inner Cylinder
    if (m_Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(std::make_shared<const Alg::Transform3D>(Alg::Transform3D(transform())),m_Rmin,m_halfZ));
        m_boundarySurfaces.push_back(innerCylinder);
    }
    //setCoordinates
    m_coordinates[Reco::CylinderVolume::Rmin]   = m_Rmin;
    m_coordinates[Reco::CylinderVolume::Rmax]   = m_Rmax;
    m_coordinates[Reco::CylinderVolume::halfZ]  = m_halfZ;
}

Reco::CylinderVolume::CylinderVolume(TGeoNode* node, TGeoConeSeg* tube) :
Reco::Volume(node)
{
    m_Rmin = tube->GetRmin1();
    m_Rmax = tube->GetRmax1();
    m_halfZ= tube->GetDZ();
    
    //Discs
    Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transform().Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transform().Translation().Vect();
    std::shared_ptr<const Alg::Transform3D> rightTransform(new Alg::Transform3D(rotation,center+lz*m_halfZ));
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    std::shared_ptr<const Alg::Transform3D> leftTransform (new Alg::Transform3D(leftrotation,center+lz*m_halfZ));
    
    std::shared_ptr<Reco::BoundaryDiscSurface> posDisc (new Reco::BoundaryDiscSurface(rightTransform,m_Rmin,m_Rmax));
    std::shared_ptr<Reco::BoundaryDiscSurface> negDisc (new Reco::BoundaryDiscSurface(leftTransform,m_Rmin,m_Rmax));
    
    //outer Cylinder
    std::shared_ptr<Reco::BoundaryCylinderSurface> outerCylinder (new Reco::BoundaryCylinderSurface(std::make_shared<const Alg::Transform3D>(Alg::Transform3D(transform())),m_Rmax,m_halfZ));
    
    m_boundarySurfaces.push_back(negDisc);
    m_boundarySurfaces.push_back(posDisc);
    m_boundarySurfaces.push_back(outerCylinder);
    
    //inner Cylinder
    if (m_Rmin>0.) {
        std::shared_ptr<Reco::BoundaryCylinderSurface> innerCylinder (new Reco::BoundaryCylinderSurface(std::make_shared<const Alg::Transform3D>(Alg::Transform3D(transform())),m_Rmin,m_halfZ));
        m_boundarySurfaces.push_back(innerCylinder);
    }
    //setCoordinates
    m_coordinates[Reco::CylinderVolume::Rmin]   = m_Rmin;
    m_coordinates[Reco::CylinderVolume::Rmax]   = m_Rmax;
    m_coordinates[Reco::CylinderVolume::halfZ]  = m_halfZ;
}

Reco::CylinderVolume::~CylinderVolume()
{}

double Reco::CylinderVolume::getRmin() const
{
    return m_Rmin;
}

double Reco::CylinderVolume::getRmax() const
{
    return m_Rmax;
}

double Reco::CylinderVolume::getHalfZ() const
{
    return m_halfZ;
}

const Reco::BoundarySurface* Reco::CylinderVolume::getPosDisc() const
{
    return (m_boundarySurfaces[Reco::CylinderVolume::posDisc].get());
}

const Reco::BoundarySurface* Reco::CylinderVolume::getNegDisc() const
{
    return (m_boundarySurfaces[Reco::CylinderVolume::negDisc].get());
}

const Reco::BoundarySurface* Reco::CylinderVolume::getOuterCylinder() const
{
    return (m_boundarySurfaces[Reco::CylinderVolume::outerCylinder].get());
}

const Reco::BoundarySurface* Reco::CylinderVolume::getInnerCylinder() const
{
    return (m_boundarySurfaces[Reco::CylinderVolume::innerCylinder].get());
}

bool Reco::CylinderVolume::isInside(const Alg::Point3D& glopos, double tol) const
{
    Alg::Point3D locpos (transform().Inverse()*glopos);
    std::cout << "CylinderVolume: " << transform() << std::endl;
    double r    = sqrt(locpos.X()*locpos.X()+locpos.Y()*locpos.Y());
    return ((r>=(m_Rmin-tol)) && (r<=(m_Rmax+tol)) && (fabs(locpos.Z()) <= (m_halfZ+tol)));
}

//std::vector<Reco::Layer*> materialLayersOrdered(const Alg::Point3D& glopos,const Alg::Vector3D& mom, double charge) const
//{
   // Reco::Layer* layer = m_layers->object(glopos);
    //type funktion machen

   // if (layer && navilayer) {
        
  //  }
    
    
//}













