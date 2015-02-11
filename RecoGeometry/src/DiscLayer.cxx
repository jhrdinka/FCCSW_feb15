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
{
    //Discs
  /*  Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transform().Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transform().Translation().Vect();
    Alg::Transform3D* rightTransform    = new Alg::Transform3D(rotation,center+lz*m_dz);
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Transform3D* leftTransform     = new Alg::Transform3D(leftrotation,center+lz*m_dz);
    
    std::shared_ptr<Reco::DiscSurface> rightDisc (new Reco::DiscSurface(rightTransform,getRmin(),getRmax()));
    std::shared_ptr<Reco::DiscSurface> leftDisc (new Reco::DiscSurface(leftTransform,getRmin(),getRmax()));
    
    m_boundaries->push_back(rightDisc);
    m_boundaries->push_back(leftDisc); */
}

Reco::DiscLayer::DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, TGeoConeSeg* tube, Trk::BinnedArray<Surface>* sf) :
DiscSurface(tube,transf),
Layer(sf),
m_dz(tube->GetDZ())
{
    //Discs
/*    Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transf->Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transf->Translation().Vect();
    Alg::Transform3D* rightTransform    = new Alg::Transform3D(rotation,center+lz*m_dz);
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Transform3D* leftTransform     = new Alg::Transform3D(leftrotation,center+lz*m_dz);
    
    std::shared_ptr<Reco::DiscSurface> rightDisc (new Reco::DiscSurface(rightTransform,getRmin(),getRmax()));
    std::shared_ptr<Reco::DiscSurface> leftDisc (new Reco::DiscSurface(leftTransform,getRmin(),getRmax()));
    
    m_boundaries->push_back(rightDisc);
    m_boundaries->push_back(leftDisc);*/
}

Reco::DiscLayer::DiscLayer(std::shared_ptr<const Alg::Transform3D> transf, double rmin, double rmax, double dz, Trk::BinnedArray<Surface>* sf) :
DiscSurface(transf, rmin, rmax),
Layer(sf),
m_dz(dz)
{
    //Discs
  /*  Alg::Vector3D lx(0.,0.,0.);
    Alg::Vector3D ly(0.,0.,0.);
    Alg::Vector3D lz(0.,0.,0.);
    Alg::RotationMatrix3D rotation      = transf->Rotation();
    rotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Vector3D         center        = transf->Translation().Vect();
    Alg::Transform3D* rightTransform    = new Alg::Transform3D(rotation,center+lz*m_dz);
    Alg::RotationMatrix3D leftrotation  = rotation*Alg::RotationY(M_PI)*Alg::RotationZ(-M_PI);
    leftrotation.GetComponents<Alg::Vector3D>(lx,ly,lz);
    Alg::Transform3D* leftTransform     = new Alg::Transform3D(leftrotation,center+lz*m_dz);
    
    std::shared_ptr<Reco::DiscSurface> rightDisc (new Reco::DiscSurface(rightTransform,rmin,rmax));
    std::shared_ptr<Reco::DiscSurface> leftDisc (new Reco::DiscSurface(leftTransform,rmin,rmax));
    
    m_boundaries->push_back(rightDisc);
    m_boundaries->push_back(leftDisc); */
}

Reco::DiscLayer::~DiscLayer()
{}

const Reco::Surface* Reco::DiscLayer::getModule(const Alg::Point3D& glopos) const
{
    Alg::Point2D* locpos = new Alg::Point2D();
    const Alg::Vector3D* mom = new Alg::Vector3D();
   if (this->globalToLocal(glopos, *mom, *locpos)) return (getModule(*locpos));
    return (0);
}

const Reco::Surface* Reco::DiscLayer::getModule(const Alg::Point2D& locpos) const
{
    return (m_surfaces->object(locpos));
}

/*bool Reco::DiscLayer::insideLayer(const Alg::Point3D& glopos, double tol=0.) const
{
    return ();
}*/

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



