//
//  Volume.cxx
//  
//
//  Created by Julia Hrdinka on 09/01/15.
//
//

#include "RecoGeometry/Volume.h"

Alg::Point3D    Reco::Volume::s_origin;
Alg::Transform3D Reco::Volume::s_idTransform;

Reco::Volume::Volume() :
m_center(0),
m_transform(0),
m_layers(0),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{}

Reco::Volume::Volume(std::shared_ptr<const Alg::Transform3D> transf, LayerArray* layers) :
m_center(0),
m_transform(transf),
m_layers(layers),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{}

Reco::Volume::Volume(std::shared_ptr<const Alg::Transform3D> transf) :
m_center(0),
m_transform(transf),
m_layers(0),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{}

Reco::Volume::Volume(LayerArray* layers) :
m_center(0),
m_transform(0),
m_layers(layers),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{}

Reco::Volume::Volume(LayerArray* layers, TGeoNode* node) :
m_center(0),
m_transform(0),
m_layers(layers),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{
    const Double_t* rotation    = (node->GetMatrix()->GetRotationMatrix());
    const Double_t* translation = (node->GetMatrix()->GetTranslation());
    
    m_transform = std::make_shared<const Alg::Transform3D>(Alg::Transform3D(
                                       rotation[0], rotation[1], rotation[2], translation[0],
                                       rotation[3], rotation[4], rotation[5], translation[1],
                                       rotation[6], rotation[7], rotation[8], translation[2]));
    
}

Reco::Volume::Volume(TGeoNode* node) :
m_center(0),
m_transform(0),
m_layers(0),
m_boundarySurfaces(),
m_volumeType(Reco::Volume::none),
m_coordinates(3,0.)
{
    const Double_t* rotation    = (node->GetMatrix()->GetRotationMatrix());
    const Double_t* translation = (node->GetMatrix()->GetTranslation());
    
    m_transform = m_transform = std::make_shared<const Alg::Transform3D>(Alg::Transform3D(
                                       rotation[0], rotation[1], rotation[2], translation[0],
                                       rotation[3], rotation[4], rotation[5], translation[1],
                                       rotation[6], rotation[7], rotation[8], translation[2]));

}

Reco::Volume::~Volume()
{
    delete m_center;
    delete m_layers;
    
}

const Alg::Point3D& Reco::Volume::center() const
{
    if (!m_transform && !m_center) return (s_origin);
    else if(!m_center) m_center = new Alg::Point3D(m_transform->Translation().Vect());
        return (*m_center);
    
}

const Alg::Transform3D& Reco::Volume::transform() const
{
    if(m_transform) return (*m_transform);
    return (s_idTransform);
}

//const LayerArray& Reco::Volume::layers() const
//{
//    return m_layers;
//}

/*Reco::Volume::SurfaceVector& Reco::Volume::boundarySurfaces()
{
    return (m_boundarySurfaces);
}*/

int Reco::Volume::NumberOfSurfaces() const
{
    return (m_boundarySurfaces.size());
}

const Reco::BoundarySurface* Reco::Volume::getBoundarySurface(int n) const
{
    return (m_boundarySurfaces.at(n).get());
}

void Reco::Volume::setType(Reco::Volume::VolumeType volumeType) const
{
    m_volumeType = volumeType;
}

Reco::Volume::VolumeType Reco::Volume::type() const
{
    return m_volumeType;
}

double Reco::Volume::getCoordinate(int n) const
{
    return (m_coordinates.at(n));
}

std::vector<const Reco::Layer*> Reco::Volume::materialLayersOrdered(const Alg::Point3D& glopos,const Alg::Vector3D& mom, double) const
{   std::vector<const Reco::Layer*> result;
    std::vector<const Reco::Layer*> fulllayers(m_layers->arrayObjects());
    if (isInside(glopos)) {
        std::cout << "###INSIDE###" << std::endl;
        const Trk::BinUtility* binutil = m_layers->binUtility();
        Alg::Point3D newpos = glopos;
        if (binutil->orderDirection(glopos,mom)==1)
        {
            std::cout << "orderdirection ==1 " << std::endl;
            int counter = 0;
            for(auto it = (fulllayers.cbegin()+binutil->bin(glopos)); it!=fulllayers.cend();++it)
            {
                double path = 0;
                //navigationlayer
                if ((*it)->type()==0) {
                    if (counter ==0) {
                        const Reco::Layer* nextlayer = ((*it)->getNextLayer());
                        if (nextlayer) {
                            path = (nextlayer)->surfaceRepresentation()->straightLineIntersection(newpos,mom.Unit()).pathlength;
                            newpos = newpos + mom.Unit()*path;
                            }
                        else return (result);
                    }
                }
                else {
                     if((*it)->onLayer(newpos)) {
                     result.push_back(*it);
                     const Reco::Layer* nextlayer = ((*it)->getNextLayer());
                     if (nextlayer) {
                        path = (nextlayer)->surfaceRepresentation()->straightLineIntersection(newpos,mom.Unit()).pathlength;
                         newpos = newpos + mom.Unit()*path;
                         }
                        else return (result);
                     }
                    else return (result);
                }
                
                ++counter;
            } //for
        } //pos direction
        //noch einbaun wenn normal drauf steht!!!
        else {
            std::reverse(fulllayers.begin(),fulllayers.end());
            
            std::cout << "orderdirection ==-1 " << std::endl;
            int counter = 0;
            int k = fulllayers.size() - (binutil->bin(glopos)+1);
    //        Alg::Point2D locpos;
   //         if(fulllayers.at(k)->surfaceRepresentation()->globalToLocal(glopos,mom,locpos)) {
  //              Alg::Vector3D normal = fulllayers[k].normal(locpos);
   //             if (normal.Dot(mom)==0) {
   //                 return (result);
   //             }
   //         }
            for(auto it = (fulllayers.cbegin()+k); it!=fulllayers.cend();++it)
            {
                
                double path = 0;
                //navigationlayer
                if ((*it)->type()==0) {
                    if (counter ==0) {
                        const Reco::Layer* nextlayer = ((*it)->getPreviousLayer());
                        if (nextlayer) {
                            path = (nextlayer)->surfaceRepresentation()->straightLineIntersection(newpos,mom.Unit()).pathlength;
                            newpos = newpos + mom.Unit()*path;
                        }
                        else return (result);
                    }
                }
                else {
                    if((*it)->onLayer(newpos)) {
                        result.push_back(*it);
                        const Reco::Layer* nextlayer = ((*it)->getPreviousLayer());
                        if (nextlayer) {
                            path = (nextlayer)->surfaceRepresentation()->straightLineIntersection(newpos,mom.Unit()).pathlength;
                            newpos = newpos + mom.Unit()*path;
                        }
                        else return (result);
                    }
                    else return (result);
                }
                
                ++counter;
            } //for

            
        }
    
    } //isInside volume
    else std::cout << "###OUTSIDE###" << std::endl;
       return (result);
}







