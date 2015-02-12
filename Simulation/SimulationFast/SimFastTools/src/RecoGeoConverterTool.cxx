//
//  RecoGeoConverter.cxx
//  
//
//  Created by Julia Hrdinka on 10/11/14.
//
//

#include "SimFastTools/RecoGeoConverterTool.h"
#include "DetExtensions/IExtension.h"
#include "DetExtensions/Extension.h"
#include "DetExtensions/DetCylinderLayer.h"
#include "DetExtensions/DetDiscLayer.h"
#include "DetExtensions/DetModule.h"
#include "DetExtensions/DetCylinderVolume.h"
#include "DetExtensions/DetDiscVolume.h"
#include "TrkGeometryUtils/BinUtility.h"
#include "TrkGeometryUtils/BinnedArray1D.h"
#include "TrkGeometryUtils/BinnedArray2D.h"
#include "RecoGeometry/CylinderLayer.h"
#include "RecoGeometry/DiscLayer.h"
#include "RecoGeometry/CylinderVolume.h"
#include "RecoGeometry/NavigationLayer.h"
#include "RecoGeometry/ContainerVolume.h"
#include "RecoGeometry/ContainerCylinderVolume.h"

#include "TFile.h"

//Gaudi random generator
#include "GaudiKernel/RndmGenerators.h"


DECLARE_COMPONENT(RecoGeoConverterTool)

// Standard Constructor
RecoGeoConverterTool::RecoGeoConverterTool(const std::string& type,
                                   const std::string& name,
                                   const IInterface* parent)
    : AlgTool( type, name, parent ),
    m_log(msgSvc(), name),
    m_lcdd (0)
{
        m_out.open("detelements.txt");
        
        
// Declare additional interface
    declareInterface<IGeoConverterTool>(this);

//Declare Properties of the specific tool
        //=> for JobOptions
    
//    if (service("THistSvc", m_ths).isFailure()) m_log << MSG::ERROR << "Couldn't get THistSvc" << endmsg;
    
}

double RecoGeoConverterTool::flatrand(double min, double max) {
    std::uniform_real_distribution<double> unif(min,max);
    std::random_device rd;
    std::mt19937 gen(rd());
    
    return (unif(gen));
}

StatusCode RecoGeoConverterTool::convert(DD4hep::Geometry::LCDD* lcdd) {
    
    m_lcdd = lcdd;
    DD4hep::Geometry::DetElement detworld = m_lcdd->world();
            translateDetector(detworld);
    m_counter=0;
//    const  DD4hep::Geometry::DetElement::Children& children = detworld.children();
    //schon hier durch children gehen, weil world detelement keine extension hat und man auch keine extension setzen kann
//    for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i != children.end(); ++i) {
      //  scanDetector((*i).second);
//    }
    
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::scanDetector(DD4hep::Geometry::DetElement detelement) {
    
    if (detelement.isValid()){
            translateDetector(detelement);
        const  DD4hep::Geometry::DetElement::Children& children = detelement.children();
        for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i != children.end(); ++i) {
                    scanDetector((*i).second);
        }
    }
    else std::cout << "detelement inValid!!!!!" << std::endl;
    
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::translateDetector(DD4hep::Geometry::DetElement detelement){
    
    //   Rndm::Numbers numbers (&*(m_rndmSvcHandle), Rndm::Flat(0,1));
   // std::cout << numbers.shoot();
    
    std::multimap<int, std::shared_ptr<const Reco::Volume>> volumes;
    std::vector<std::shared_ptr<const Reco::Volume>> previousvol(4,0);
    std::vector<std::shared_ptr<const Reco::Volume>> currentvol(4,0);
    std::vector<float> bZValues(Reco::ContainerCylinderVolume::zLength);
    std::vector<float> bRValues(Reco::ContainerCylinderVolume::rLength);
    std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>* binRvector = new std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>(2);
    std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>* binZvector = new std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>(3);
//    Trk::BinnedArray1D<Reco::Volume>* binnedvolumes;
    //translate rest of detector
    translateVolume(detelement, volumes);
    //putting keys in a vector;
    std::vector<std::pair<int,std::shared_ptr<const Reco::Volume>>> keys_dedup;
    std::unique_copy(begin(volumes),
                end(volumes),
                back_inserter(keys_dedup),
                     [](const std::pair<int,std::shared_ptr<const Reco::Volume>> &entry1,
                        const std::pair<int,std::shared_ptr<const Reco::Volume>> &entry2) {
                    return (entry1.first == entry2.first);}
                     );
    std::cout << "KEYS SIZE: " << keys_dedup.size() << std::endl;
    for (const auto& itkeys : keys_dedup) {
        std::pair <std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator, std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator> ret;
        ret = volumes.equal_range(itkeys.first);
        std::cout << "itkeys" << itkeys.first << std::endl;
 //       std::cout << "last element: " << ret.second.first << std::endl;
        for (std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator itvol=ret.first; itvol!=ret.second; ++itvol) {
            std::cout << "itvol type: " << itvol->second->type() << std::endl;
            std::shared_ptr<const Reco::Volume> volume_ptr(itvol->second);
            currentvol.at(itvol->second->type()) = itvol->second;
        } //for volumevectors fill
        
       if (itkeys.first==0) {
            std::cout << "itkeys==0" << std::endl;
 //           Alg::Point3D center = currentvol.at(Reco::Volume::barrel)->center();
//            binRvector->at(0)  = (make_pair(currentvol.at(Reco::Volume::barrel),center));
            previousvol = currentvol;
        }
        
        else if (itkeys.first==1) {
            std::cout << "itkeys==1" << std::endl;
            //falls ich mehr volumes mache, unterscheidung in volume type
            //set Volumes for BoundarySurfaces between the volumes
            currentvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            std::cout << "Test1" << std::endl;
              //beamtube
            if (!previousvol.empty()) {
                std::cout << "!empty" << std::endl;
            currentvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            }
            std::cout << "Test2" << std::endl;
            //was wenn volumes nicht genau aufeinander treffen??
            Alg::Point3D center1(0.,0.,0.); // = new Alg::Point3D();
            Alg::Point3D center2(0.,0.,0.);
            double halfZ = 0;
            double conthalfZ = 0;
            double contRmax = 0;
            double contRmin = 0;
            //fill bValues
            //neg EndCap
            center1 = currentvol.at(Reco::Volume::nEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nEnd) = (center1.Z()-halfZ);
            binZvector->at(0) = (make_pair(currentvol.at(Reco::Volume::nEndCap),center1));
            std::cout << "Test3" << std::endl;
            //barrel
            center1 = currentvol.at(Reco::Volume::barrel)->center();
            halfZ   = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nBarrel) = (center1.Z()-halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pBarrel) = (center1.Z()+halfZ);
            binZvector->at(1) = (make_pair(currentvol.at(Reco::Volume::barrel),center1));
            std::cout << "Test4" << std::endl;
            //pos EndCap
            center2 = currentvol.at(Reco::Volume::pEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::pEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pEnd) = (center2.Z()+halfZ);
            binZvector->at(2) = (make_pair(currentvol.at(Reco::Volume::pEndCap),center2));
            std::cout << "Test5" << std::endl;
            //length of the ContainerVolume
            conthalfZ = sqrt((center2.Z()-center1.Z())*(center2.Z()-center1.Z()))+halfZ;
            //use Rmin & Rmax of barrel eventuell abmessungen bei Hierarchie hinzufuegen??
            contRmin = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmin);
            contRmax = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmax);
            std::shared_ptr<const Alg::Transform3D> conttrans(new Alg::Transform3D(currentvol.at(Reco::Volume::barrel)->transform()));
            Trk::BinUtility* binutility = new Trk::BinUtility(bZValues,Trk::open,Trk::binZ);
            Trk::BinnedArray1D<Reco::Volume>* binnedvolumes = new Trk::BinnedArray1D<Reco::Volume>(*binZvector,binutility);
            std::shared_ptr<const Reco::ContainerCylinderVolume> contvol(new Reco::ContainerCylinderVolume(binnedvolumes,conttrans,contRmin,contRmax,conthalfZ));
            if (contvol) {
            std::cout << "Test6" << std::endl;
            std::cout << "containervolume" << std::endl;
            m_out << "created containervolume" << std::endl;
            //beampipe
            contvol->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            //rest
            contvol->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedvolumes));
            contvol->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedvolumes));
            contvol->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::nEndCap));
            contvol->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(currentvol.at(Reco::Volume::pEndCap));
            contvol->setType(Reco::Volume::container);
            std::cout << "Test7" << std::endl;
            currentvol.at(Reco::Volume::container) = contvol;
            }
            std::cout << "Test8" << std::endl;

            previousvol = currentvol;
            //add containervolume in the end
            std::cout << "Test9" << std::endl;
        } //if status==1
        else {
            std::cout << "else " << std::endl;
            if(!previousvol.at(Reco::Volume::container)) std::cout << "no container " << std::endl;
                //->getBoundarySurface(Reco::CylinderVolume::outerCylinder);
            //set Volumes for BoundarySurfaces
        //previous container volume
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
        //previous volumes
                //barrel
            std::cout << "2Test1" << std::endl;
                previousvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                previousvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                previousvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                //EndCaps
            std::cout << "2Test2" << std::endl;
                previousvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
                previousvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
            
        //current volumes
            //barrel
            std::cout << "2Test3" << std::endl;
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolumes(previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->getPreviousVolumes());
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
            //EndCaps
             //Grenzen muessen genau ueberein stimmen!
             //make binned Array
            //binutility
            std::cout << "2Test4" << std::endl;
            bRValues.at(Reco::ContainerCylinderVolume::inner) = 0.; //previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmin);
            bRValues.at(Reco::ContainerCylinderVolume::middle) = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmax);
            bRValues.at(Reco::ContainerCylinderVolume::outer) = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmax);
            std::cout << "2Test5" << std::endl;
            Trk::BinUtility* binutil = new Trk::BinUtility(bRValues,Trk::open,Trk::binR);
            Alg::Point3D center1(0.,0.,0.);
            Alg::Point3D center2(0.,0.,0.);
            std::cout << "2Test6" << std::endl;
            center1 = currentvol.at(Reco::Volume::barrel)->center();
            binRvector->at(1)  = (make_pair(currentvol.at(Reco::Volume::barrel),center1));
            std::cout << "2Test6.1" << center1 << std::endl;
            //negativ
            center1 = previousvol.at(Reco::Volume::nEndCap)->center();
            binRvector->at(0) = (make_pair(previousvol.at(Reco::Volume::nEndCap),center1));
            std::cout << "2Test6.2" << center1 << std::endl;
            Trk::BinnedArray1D<Reco::Volume>* binnedRvolumes1 = new Trk::BinnedArray1D<Reco::Volume>(*binRvector,binutil);
            if (!currentvol.at(Reco::Volume::nEndCap)) {
                std::cout << "not" << std::endl;
            }
            //->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolumes(binnedRvolumes1);
            //positiv
            std::cout << "2Test7" << std::endl;
            center1 = previousvol.at(Reco::Volume::pEndCap)->center();
            binRvector->at(0) = (make_pair(previousvol.at(Reco::Volume::pEndCap),center1));
            std::cout << "2Test8" << std::endl;
            Trk::BinnedArray1D<Reco::Volume>* binnedRvolumes2 = new Trk::BinnedArray1D<Reco::Volume>(*binRvector,binutil);
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes2));
            std::cout << "2Test9" << std::endl;
        //make ContainerVolume1 //momentan wird containervolume an containervolume 1 uebergeben und nicht schon aufgeloest in untervolumen -gut? ABER FUER NEext und previousvolumes is schon aufegloest, also passt
            center1 = previousvol.at(Reco::Volume::container)->center();
            binRvector->at(0) = make_pair(previousvol.at(Reco::Volume::container),center1);
            double halfZ = 0;
            double conthalfZ = 0;
            double contRmin = 0;
            double contRmax = 0;
            contRmin = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmin);
            contRmax = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmax);
            std::shared_ptr<const Alg::Transform3D> conttrans(new Alg::Transform3D(previousvol.at(Reco::Volume::container)->transform()));
            conthalfZ = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::halfZ);
            Trk::BinnedArray1D<Reco::Volume>* containerbin1 = new Trk::BinnedArray1D<Reco::Volume>(*binRvector,binutil);
            std::shared_ptr<const Reco::ContainerCylinderVolume> contvol1(new Reco::ContainerCylinderVolume(containerbin1,conttrans,contRmin,contRmax,conthalfZ));
            std::cout << "2Test10" << std::endl;
            //set surfacepointers of ContainerVolume1
            contvol1->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->getNextVolumes());
            contvol1->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            contvol1->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes1));
            contvol1->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes2));
            //beamtube
      //      contvol1->getBoundarySurface(Reco::CylinderVolume::innerCylinder) =

            std::cout << "2Test11" << std::endl;
            //make ContainerVolume2
            //fill bZValues
            //neg EndCap
            center1 = currentvol.at(Reco::Volume::nEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nEnd) = (center1.Z()-halfZ);
            binZvector->at(0) = (make_pair(currentvol.at(Reco::Volume::nEndCap),center1));
            //barrel
            std::cout << "2Test12" << std::endl;
            center1 = contvol1->center();
            halfZ   = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nBarrel) = (center1.Z()-halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pBarrel) = (center1.Z()+halfZ);
            binZvector->at(1) = (make_pair(std::shared_ptr<const Reco::Volume>(contvol1),center1));
            //pos EndCap
            std::cout << "2Test13" << std::endl;
            center2 = currentvol.at(Reco::Volume::pEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::pEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pEnd) = (center2.Z()+halfZ);
            binZvector->at(2) = (make_pair(currentvol.at(Reco::Volume::pEndCap),center2));
            
            Trk::BinUtility* binZutil = new Trk::BinUtility(bZValues, Trk::open, Trk::binZ);
            conthalfZ = sqrt((center2.Z()-center1.Z())*(center2.Z()-center1.Z()))+halfZ;
            Trk::BinnedArray1D<Reco::Volume>* containerbin2 = new Trk::BinnedArray1D<Reco::Volume>(*binZvector,binZutil);
            std::shared_ptr<const Reco::ContainerCylinderVolume> contvol2(new Reco::ContainerCylinderVolume(containerbin2,conttrans,contRmin,contRmax,conthalfZ));
            //make new binnedarray
            std::cout << "2Test14" << std::endl;
            center1 = currentvol.at(Reco::Volume::barrel)->center();
            binZvector->at(1) = (make_pair(std::shared_ptr<const Reco::Volume>(currentvol.at(Reco::Volume::barrel)),center1));
            Trk::BinnedArray1D<Reco::Volume>* containerbin3 = new Trk::BinnedArray1D<Reco::Volume>(*binZvector,binZutil);
            //set surfacepointers of ContainerVolume2
            std::cout << "2Test15" << std::endl;
            const Trk::BinnedArray1D<Reco::Volume>* containerbin4 = previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->getNextVolumes();
            contvol2->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(new Trk::BinnedArray1D<Reco::Volume>(*containerbin4));
            contvol2->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*containerbin3));
            contvol2->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::nEndCap));
            contvol2->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(currentvol.at(Reco::Volume::pEndCap));

            previousvol = currentvol;
            std::cout << "2Test16" << std::endl;
        }
      std::cout << "2Test17" << std::endl;
     //   currentvol.clear();
    } //for keys
    
    



    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::translateVolume(DD4hep::Geometry::DetElement det, std::multimap<int, std::shared_ptr<const Reco::Volume>>& volumes)

{
    const DD4hep::Geometry::DetElement::Children& children = det.children();
    for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i!=children.end(); ++i) {
        
    DD4hep::Geometry::DetElement detelement = (*i).second;
    Det::IExtension* ex = detelement.extension<Det::IExtension>();
    Det::DetCylinderVolume* detcylindervolume = dynamic_cast<Det::DetCylinderVolume*>(ex);
    Det::DetDiscVolume* detdiscvolume = dynamic_cast<Det::DetDiscVolume*>(ex);
  
    if (detcylindervolume) {
        std::cout << "cylindervolume" << std::endl;
        int status = detcylindervolume->status();
        DD4hep::Geometry::PlacedVolume placedvolume = detelement.placement();
        TGeoNode* geonode = placedvolume.ptr();
        if (geonode) {
            TGeoShape* geoshape = geonode->GetVolume()->GetShape();
            const Double_t* rotation    = (geonode->GetMatrix()->GetRotationMatrix());
            const Double_t* translation = (geonode->GetMatrix()->GetTranslation());
            std::shared_ptr<const Alg::Transform3D> transform(new Alg::Transform3D(rotation[0], rotation[1], rotation[2], translation[0],
                                                               rotation[3], rotation[4], rotation[5], translation[1],
                                                               rotation[6], rotation[7], rotation[8], translation[2]));
            
            if(geoshape && geoshape->IsA() == TGeoConeSeg::Class()){
                TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoshape);
                if(tube){
                    LayerVector layers;
                    //translateLayers
                    m_counter = 0;
                    translateLayer(detelement,layers, transform);
                    if(!layers.empty()){
                        double Rmax = tube->GetRmax1();
        //                double halfZ = tube->GetDz();
                        double x = 0;
                        double y = 0;
                        double z = 0;
                        transform->Translation().GetComponents(x,y,z);
                        Alg::Point3D center(x,y,z);
                        std::vector<float> bValues;
                        LayerVector* fulllayers = new LayerVector();
                        binCylinderLayers(layers, *fulllayers, bValues, center, Rmax);
                        if (bValues.empty()||fulllayers->empty()) return::StatusCode::FAILURE;
                        else {
                            std::cout << "size fulllayers: " << fulllayers->size() << std::endl;
                            Trk::BinUtility* binutility = new Trk::BinUtility(bValues,Trk::open,Trk::binR);
                            if(binutility) {
                                Trk::BinnedArray1D<Reco::Layer>* binnedlayers = new Trk::BinnedArray1D<Reco::Layer>(*fulllayers,binutility);
                                std::shared_ptr<const Reco::Volume> cylindervolume(new Reco::CylinderVolume(binnedlayers, geonode, tube));
                                if (cylindervolume) {
   //a                                 m_out << "##### created Volume #####" << std::endl;
                                    std::cout  << "##### created cylindervolume" << std::endl;
                              /*      Alg::Point3D pos(1.,0.,0.);
                                    Alg::Vector3D mom(0.,0.,1.);
                                    std::vector<const Reco::Layer*> materiallayers = cylindervolume->materialLayersOrdered(pos,mom,2.);
                                    std::cout << "CylinderVolume materiallayersize: " << materiallayers.size() << std::endl;*/
                                    
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }

                                    cylindervolume->setType(Reco::Volume::barrel);
                                    std::cout << "status" << status << std::endl;
                                    volumes.emplace(status,cylindervolume);
                                } //if cylindervolume
                            }
                            else return StatusCode::FAILURE;
                        }
                    }//if !layer.empty()
                    else {
                        std::shared_ptr<const Reco::Volume> cylindervolume(new Reco::CylinderVolume(geonode, tube));
                                if (cylindervolume) {
                                    std::cout  << "##### created cylindervolume" << std::endl;
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }
                                    cylindervolume->setType(Reco::Volume::barrel);
                                    std::cout << "status" << status << std::endl;
                                    volumes.emplace(status,cylindervolume);
                                } //if cylindervolume
                            else return StatusCode::FAILURE;
                    } // layer.empty()
                }//if tube
            }//if geoshape
        } //if geonode
        
    } //if detcylindervolume
    
    if (detdiscvolume) {
        int status = detdiscvolume->status();
        DD4hep::Geometry::PlacedVolume placedvolume = detelement.placement();
        TGeoNode* geonode = placedvolume.ptr();
        if (geonode) {
            TGeoShape* geoshape = geonode->GetVolume()->GetShape();
            const Double_t* rotation    = (geonode->GetMatrix()->GetRotationMatrix());
            const Double_t* translation = (geonode->GetMatrix()->GetTranslation());
            std::shared_ptr<const Alg::Transform3D> transform(new Alg::Transform3D(rotation[0], rotation[1], rotation[2], translation[0],
                                                               rotation[3], rotation[4], rotation[5], translation[1],
                                                               rotation[6], rotation[7], rotation[8], translation[2]));
            std::cout << *transform << std::endl;
            if(geoshape && geoshape->IsA() == TGeoConeSeg::Class()){
                TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoshape);
                if(tube){
                    LayerVector layers;          //translateLayers
                    m_counter = 0;
                    translateLayer(detelement,layers, transform);
                    if(!layers.empty()){
                        double halfZ = tube->GetDz();
                        std::cout << "### Disc Rmin: " << tube->GetRmin1() << std::endl;
                        double x = 0;
                        double y = 0;
                        double z = 0;
                        transform->Translation().GetComponents(x,y,z);
                        Alg::Point3D center(x,y,z);
                        std::vector<float> bValues;
                        LayerVector* fulllayers = new LayerVector();
                        binDiscLayers(layers, *fulllayers, bValues, center, halfZ);
                        if (bValues.empty()||fulllayers->empty()) return::StatusCode::FAILURE;
                        else {
                            std::cout << "size fulllayers: " << fulllayers->size() << std::endl;
                            std::cout << "bValues size: " << bValues.size() << std::endl;
                            Trk::BinUtility* binutility = new Trk::BinUtility(bValues,Trk::open,Trk::binZ);
                            std::cout << "binutility size: " << binutility->dimensions() << std::endl;;
                            if(binutility) {
                                Trk::BinnedArray1D<Reco::Layer>* binnedlayers = new Trk::BinnedArray1D<Reco::Layer>(*fulllayers,binutility);
                                std::shared_ptr<const Reco::CylinderVolume> cylindervolume(new Reco::CylinderVolume(binnedlayers, geonode, tube));
                                if (cylindervolume) {
  //a                                  m_out << "##### created Volume #####" << std::endl;
                             /*       Alg::Point3D pos(5.,0.,-40.);
                                    Alg::Vector3D mom(1.,0.,0.);
                                    std::vector<const Reco::Layer*> materiallayers = cylindervolume->materialLayersOrdered(pos,mom,2.);
                                    std::cout << "DiscVolume materiallayersize: " << materiallayers.size() << std::endl; */
 /*
                                    std::cout << "##DiscVolume Rmin: " << cylindervolume->getRmin() << std::endl;
                                    std::shared_ptr<const Reco::CylinderVolume> volume_ptr (cylindervolume);
                                    std::pair<std::shared_ptr<const Reco::CylinderVolume>, Alg::Point3D> volume (volume_ptr,center);
                                    if (center.Z()<0.) {
                                        volumes[Reco::ContainerCylinderVolume::nEndCap] = volume;
                                        bVolValues[Reco::ContainerCylinderVolume::nEnd] = (center.Z()-halfZ);
                                    }
                                    else {
                                        volumes[Reco::ContainerCylinderVolume::pEndCap] = volume;
                                        bVolValues[Reco::ContainerCylinderVolume::pEnd] = (center.Z()+halfZ);
                                    }
  */
                                    std::cout << cylindervolume->transform() << std::endl;
                                    
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }
                                    if (center.Z()<0.) cylindervolume->setType(Reco::Volume::nEndCap);
                                    else cylindervolume->setType(Reco::Volume::pEndCap);
                                    std::cout << "status" << status << std::endl;
                                    volumes.emplace(status,cylindervolume);
                                } //if cylindervolume
                            }
                            else return StatusCode::FAILURE;
                        }
                    }//if !layer.empty()
                }//if tube
            }//if geoshape
        } //if geonode
        
    } //if detdiscvolume

  //  else volumes.clear();
    }
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::translateLayer(DD4hep::Geometry::DetElement det, LayerVector& layers, std::shared_ptr<const Alg::Transform3D>transform)

{
    const  DD4hep::Geometry::DetElement::Children& children = det.children();
    for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i != children.end(); ++i)
    {
        DD4hep::Geometry::DetElement detelement = (*i).second;
        Det::IExtension* ext = detelement.extension<Det::IExtension>();
        Det::DetCylinderLayer* detcylinderlayer = dynamic_cast<Det::DetCylinderLayer*>(ext);
        //CylinderLayer
        if (detelement.isValid() && detcylinderlayer) {
            std::cout << "detcylinderlayer" << std::endl;
        DD4hep::Geometry::PlacedVolume placedvolume = detelement.placement();
        TGeoNode* geonode = placedvolume.ptr();
            if (geonode) {
                TGeoShape* geoshape = geonode->GetVolume()->GetShape();
                const Double_t* rotation    = (geonode->GetMatrix()->GetRotationMatrix());
                const Double_t* translation = (geonode->GetMatrix()->GetTranslation());
                std::shared_ptr<Alg::Transform3D> transf(new Alg::Transform3D(rotation[0], rotation[1], rotation[2], translation[0],
                                                                rotation[3], rotation[4], rotation[5], translation[1],
                                                                rotation[6], rotation[7], rotation[8], translation[2]));
                (*transf) = (*transform)*(*transf);
                if(geoshape && geoshape->IsA() == TGeoConeSeg::Class()){
                    TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoshape);
                    if (tube) {
                        double dz = tube->GetDz();
                        //now create BinnedArray with BinUtility for a Cylinder and Surfaces, for the Layer
                        int binsPhi = detcylinderlayer->modulesPhi();
                        int binsZ   = detcylinderlayer->modulesZ();
                        Trk::BinUtility* currentBinUtility = new Trk::BinUtility(binsPhi,-M_PI,M_PI,Trk::closed,Trk::binPhi);
                        (*currentBinUtility) += Trk::BinUtility(binsZ,-dz,dz,Trk::open,Trk::binZ);
                        std::vector<std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D>> surfaces;
                        //translateModule
                        translateModule(detelement,surfaces,transf);
                        if (surfaces.empty()) {
                            std::cout << "surfaces empty" << std::endl;
                        }
                        if (!surfaces.empty()) {
                            Trk::BinnedArray2D<Reco::Surface>* bin = new Trk::BinnedArray2D<Reco::Surface>(surfaces,currentBinUtility);
                            std::shared_ptr<const Reco::CylinderLayer> cyllayer(new Reco::CylinderLayer(transf,tube,bin));
                            if (cyllayer) {
 //a                               m_out << "### Created CylinderLayer ###" << std::endl;
//a                                m_out << "CylinderLayer RMin: " << cyllayer->getRmin() << std::endl;
                                double R = 0.5*(cyllayer->getRmax()+cyllayer->getRmin());
                                Alg::Point3D center(R,0.,0.);
                                std::cout << "layer type: " << cyllayer->type() << std::endl;
//a                                m_out << "center: "<< center << std::endl;
                                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> layer (cyllayer,center);
                                layers.push_back(layer);
                            } //if cyllayer
                        } //if surfaces filled
                    } //if tube
                            
                } //if shape
                            //           transform = 0;
            } //if node
                ++m_counter;
        } //if detcyinderlayer

        Det::DetDiscLayer* detdisclayer = dynamic_cast<Det::DetDiscLayer*>(ext);
        //DiscLayer
        if (detelement.isValid() && detdisclayer) {
            std::cout << "detdisclayer" << std::endl;
            DD4hep::Geometry::PlacedVolume placedvolume = detelement.placement();
            TGeoNode* geonode = placedvolume.ptr();
            if (geonode) {
                TGeoShape* geoshape = geonode->GetVolume()->GetShape();
                const Double_t* rotation    = (geonode->GetMatrix()->GetRotationMatrix());
                const Double_t* translation = (geonode->GetMatrix()->GetTranslation());
                std::shared_ptr<Alg::Transform3D> transf(new Alg::Transform3D(rotation[0], rotation[1], rotation[2], translation[0],
                                                                rotation[3], rotation[4], rotation[5], translation[1],
                                                                rotation[6], rotation[7], rotation[8], translation[2]));
                (*transf) = (*transform)*(*transf);
                if(geoshape && geoshape->IsA() == TGeoConeSeg::Class()){
                    //weitere if Bedingung zur Unterscheidung zu Disc -> inner r und z klein??
                    TGeoConeSeg* disc = dynamic_cast<TGeoConeSeg*>(geoshape);
                    if (disc) {
                        double rmin = disc->GetRmin1();
                        double rmax = disc->GetRmax1();
                        int binsPhi = detdisclayer->modulesPhi();
                        int binsR   = detdisclayer->modulesR();
                        //now create BinnedArray with BinUtility for a Cylinder and Surfaces, for the Layer
                        Trk::BinUtility* currentBinUtility = new Trk::BinUtility(binsPhi,-M_PI,M_PI,Trk::closed,Trk::binPhi);
                        (*currentBinUtility) += Trk::BinUtility(binsR,rmin,rmax,Trk::open,Trk::binR);
                        std::vector<std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D>> surfaces;
                        //translateModule
                        translateModule(detelement,surfaces,transf);
                        if (!surfaces.empty()) {
                            Trk::BinnedArray2D<Reco::Surface>* bin = new Trk::BinnedArray2D<Reco::Surface>(surfaces,currentBinUtility);
                            std::shared_ptr<const Reco::DiscLayer> disclayer(new Reco::DiscLayer(transf,disc,bin));
                            if (disclayer) {
 //a                               m_out << "### Created DiscLayer ###" << std::endl;
                                Alg::Point3D center = disclayer->center();
//a                                m_out << "surfacecenter: "<< center << *transform<< std::endl;
                                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> layer (disclayer,center);
                                layers.push_back(layer);
                            }//if disclayer
                        }
                    } //if disc
                                        
                } //if shape
        //           transform = 0;
            } //if node
            ++m_counter;
        } //if detdisclayer
    }//for children
    if(layers.empty()) std::cout << "layers empty 2"  << std::endl;
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::translateModule(DD4hep::Geometry::DetElement det, std::vector<std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D>>& surfaces, std::shared_ptr<const Alg::Transform3D> transform)
{
    const  DD4hep::Geometry::DetElement::Children& children = det.children();
    int m = 0;
    for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i != children.end(); ++i)
    {
        m++;
        //           m_out << "Child from Layer: " << m_counter << ", " << (*i).first << std::endl;
        DD4hep::Geometry::PlacedVolume pv = (*i).second.placement();
        Det::IExtension* ext = (*i).second.extension<Det::IExtension>();
        Det::DetModule* detm = dynamic_cast<Det::DetModule*>(ext);
        DD4hep::Geometry::DetElement detelement = (*i).second;
        //MODULE
        if  (detm) {
            //  m_out << "Module" << std::endl;
            //transform *= transform
//a            m_out << "Child from Layer: " << m_counter << ", " << (*i).first << std::endl;
            TGeoNode* node = pv.ptr();
      //      TGeoShape* shape = new TGeoShape();
     //       TGeoBBox* box = new TGeoBBox();
            double halflength = 0.;
            double halfwidth = 0.;
            if (node)
            {
                TGeoShape* shape = node->GetVolume()->GetShape();
                //transform matrix module
                const Double_t* rotation    = (node->GetMatrix()->GetRotationMatrix());
                const Double_t* translation = (node->GetMatrix()->GetTranslation());
                std::shared_ptr<Alg::Transform3D> transf(new Alg::Transform3D(rotation[0], rotation[1], rotation[2], translation[0],
                                                                rotation[3], rotation[4], rotation[5], translation[1],
                                                                rotation[6], rotation[7], rotation[8], translation[2]));

                *transf = (*transform)*(*transf);

                if(shape && shape->IsA() == TGeoBBox::Class())
                {
                    TGeoBBox* box = dynamic_cast<TGeoBBox*>(shape);
                    if (box)
                    {
                        halflength    = box->GetDX();
                        halfwidth     = box->GetDY();
                        int binsX = 100;
                        int binsY = 100;
                        Trk::BinUtility* binutility = new Trk::BinUtility(binsX,-halflength,halflength,Trk::open,Trk::binX);
                        (*binutility) += Trk::BinUtility(binsY,-halfwidth,halfwidth,Trk::open,Trk::binY);
                        //map for the material
                        std::map<std::pair<int,int>, Reco::Material>* map = new std::map<std::pair<int,int>, Reco::Material>();
                        //Histogram
     /*                   m_out << "length: " << halflength << std::endl;
                        m_out << "width: " << halfwidth << std::endl;
                        m_out << "thickness: " << halfthickness << std::endl;
    */
                        //Histograms
     /*                   m_t = new TH2F("Thickness", "Thickness Histogram of a Module", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("Thickness", m_t).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register Thickness Histogram" << endmsg;
                        }
                        m_sens = new TH2F("Sensitive Material", "Percentage of Sensitive Material", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("Sensitive Material", m_sens).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register Sensitive Material Histogram" << endmsg;
                        }
                        m_dens = new TH2F("Density", "Mean Density of the Module Material", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("Density", m_dens).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register Density Histogram" << endmsg;
                        }
                        m_tx0 = new TH2F("T in x0", "Thickness in x0", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("T in x0", m_tx0).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register T in x0 Histogram" << endmsg;
                        }
                        m_tlambda0 = new TH2F("T in lambda0", "Thickness in lambda0", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("T in lambda0", m_tlambda0).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register T in lambda0 Histogram" << endmsg;
                        }
                        m_A = new TH2F("A", "Mass number", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("A", m_A).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register mass number Histogram" << endmsg;
                        }
                        m_Z = new TH2F("Z", "Atomic number", binsX, -halflength, halflength, binsY, -halfwidth, halfwidth);
                        if (m_ths->regHist("Z", m_Z).isFailure()) {
                            m_log << MSG::ERROR << "Couldn't register atomic number Histogram" << endmsg;
                        }
      */
                        double newx = 0.;
                        double newy = 0.;
                        
                        for (int m=0; m<binsX; m++)
                        {
                            newx = binutility->bincenter(m,0);
                            
                            for (int n=0; n<binsY; n++)
                            {
                                //capacity of all module components
                                double capacity = 0;
                                //senitive percentage
                                double sensper = 0;
                                
                                newy = binutility->bincenter(n,1);
                                double z = 1.;
                                double master[3]={newx,newy,z};
                                double local[3]  = {0.,0.,0.};
                                double t_x0 = 0.;
                                double t_lambda0 = 0.;
                                double A = 0.;
                                double Z = 0.;
                                double density = 0;
                                double sumt = 0.;
                                double sumdens = 0;
                                
                                double t = 0.;
                                int comp_num = 0;
                                const  DD4hep::Geometry::DetElement::Children& children = detelement.children();
                                for (DD4hep::Geometry::DetElement::Children::const_iterator j=children.begin(); j != children.end(); ++j)
                                {
                                    //std::cout << "module7" << std::endl;
                                    DD4hep::Geometry::Volume vol   = (*j).second.volume();
                                    DD4hep::Geometry::Material mat = vol.material();
                                    //calculate the whole volume of all components
                                    capacity += vol.ptr()->Capacity();
                                    //calculate the sensitive volume
                                    if(vol.isSensitive()) sensper += vol.ptr()->Capacity();
                                    TGeoNode* childnode = (*j).second.placement().ptr();
                                    if (childnode) {
                                        childnode->MasterToLocal(master,local);
        /*                                if (m==0 && n==0) std::cout << "(" << master[0] << "," << master[1] << "," << master[2] << ") " << " (" << local[0] <<  "," << local[1] << "," << local[2] << ")" << std::endl;
        */
                                        TGeoShape* childshape = childnode->GetVolume()->GetShape();
                                        
                                        if (childshape && childshape->IsA()==TGeoBBox::Class()) {
                                            TGeoBBox* childbox = dynamic_cast<TGeoBBox*>(childshape);
                                            if ((fabs(local[0])<=(childbox->GetDX())) && (fabs(local[1])<=(childbox->GetDY()))) {
                                                //       if (m==0 && n==0)std::cout << " inside " << std::endl;
       //                                         m_out << "inside " << std::endl;
                                                t          = 2*(childbox->GetDZ());
                                                sumt      += t;
      //a                                          if (m==0 && n==0) m_out << childbox->GetDZ() << std::endl;
                                                //         m_out << "##### in component ##### " << "t: " << t << "sumt " << sumt;
         
                                                t_x0      += t/mat.radLength();
                                                t_lambda0 += t/mat.intLength();
                                                A         += mat.density()*mat.A();
                                                Z         += mat.density()*mat.Z();
                                                density   += mat.density()*t;
                                                sumdens   += mat.density();
                                            }
                                        }
                                        
                                    }
                                    
//a                                    if (m==0 && n==0) m_out << " component" << comp_num << " t: " << t << "sumt: " << sumt << std::endl;
                                    ++comp_num;
                                } //for components
                                
                                sensper = sensper/capacity;
                                density = density/sumt;
                                A       = A/sumdens;
                                Z       = Z/sumdens;
                                //fill map with material for each bin
                                map->emplace(std::make_pair(m,n), Reco::Material(A, Z, density, t_x0, t_lambda0,sensper));
    /*
                                m_t->Fill(newx,newy,sumt);
                                m_out << m << " " << n << "            " << newx << "        " << newy << "        " << sumt << std::endl;
                                m_sens->Fill(newx,newy,sensper);
                                m_dens->Fill(newx,newy,density);
                                m_tx0->Fill(newx,newy,t_x0);
                                m_tlambda0->Fill(newx,newy,t_lambda0);
                                m_A->Fill(newx,newy,A);
                                m_Z->Fill(newx,newy,Z);
    */
                            } // for binsY
                        } //for binsX
                        //create MaterialMap
                        Reco::MaterialMap* materialmap = new Reco::MaterialMap(binutility, map);
                        //create Surface
                        std::shared_ptr<const Reco::PlaneSurface> plane(new Reco::PlaneSurface(box,materialmap,transf));
                        if(plane) {
//a                            m_out << "#created plane#" << std::endl;
                            Alg::Point3D center = plane->center();
                            std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D> surface (plane,center);
                            surfaces.push_back(surface);
                        }
                        
                    } //if box
                    
                } //ifshape
                
            } //ifnode
            
            //Histograms
    /*        TFile* file = new TFile("histogram.root","RECREATE");
            if (!file->IsOpen()) m_log << MSG::WARNING << "File couldn't be opened" << endmsg;
            //   m_t->Draw();
            m_t->Write();
            m_sens->Write();
            m_dens->Write();
            m_tx0->Write();
            m_tlambda0->Write();
            m_A->Write();
            m_Z->Write();
            file->Print();
  */
            
        } //if detector module
        
        //WENN GLEICH SENSITIVE MATERIALIEN ONHE MODULE
        /*              if (pv.volume().isSensitive())
         {
         m_out << ", " << "sensitive" << std::endl;
         TGeoNode* node = pv.ptr();
         if (node)
         {
         TGeoShape* shape = node->GetVolume()->GetShape();
         //hier vl fkt fuer die versch surfaces machen und auch dann jeweils speichern
         if(shape && shape->IsA() == TGeoBBox::Class())
         {
         TGeoBBox* box = dynamic_cast<TGeoBBox*>(shape);
         //if (box) std::cout << "box - dynamic cast was successfull" << std::endl;
         if (box) {
         Reco::PlaneSurface* plane = new Reco::PlaneSurface(node, box);
         if(plane) m_out << "Created Plane" << std::endl;
         Alg::Vector3D center = plane->center();
         std::shared_ptr<const Reco::Surface> surf (plane);
         std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D> surface (surf,center);
         surfaces.push_back(surface);
         }
         }
         }
         }*/
    
    } //for layer children
    
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::binCylinderLayers(LayerVector& layers, LayerVector& fulllayers, std::vector<float>& bValues, Alg::Point3D, double Rmax) const
{
    std::cout << "cylinderlayers.size(): " << layers.size() << std::endl;
    //only one layer
    if (layers.size()==1) {
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> currentpair = layers.at(0);
        //cylinder
        std::shared_ptr<const Reco::CylinderLayer> current = std::dynamic_pointer_cast<const Reco::CylinderLayer> (currentpair.first);
        if (current) {
            double currentRmin = current->getRmin();
            double currentRmax = current->getRmax();
            
            if (currentRmin==0. && currentRmax==Rmax) {
 //               std::cout << "Cylinder1Test1" << std::endl;
                bValues.push_back(currentRmin);
                bValues.push_back(currentRmax);
                fulllayers.push_back(currentpair);
            }
            if (currentRmin==0. && currentRmax<Rmax) {
 //               std::cout << "Cylinder1Test2" << std::endl;
                //bins
                bValues.push_back(currentRmin);
                bValues.push_back(currentRmax);
                bValues.push_back(Rmax);
                //Navilayer
                std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                Alg::Point3D navicenter(0.5*(Rmax+currentRmax), 0.,0.);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                //set and fill layers
                navilayer->setPreviousLayer(current);
                fulllayers.push_back(currentpair);
                fulllayers.push_back(navilayer_pair);
            }
            if (currentRmin>0. && currentRmax==Rmax) {
  //              std::cout << "Cylinder1Test3" << std::endl;
                bValues.push_back(0.);
                bValues.push_back(currentRmin);
                bValues.push_back(currentRmax);
                //Navilayer
                std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                Alg::Point3D navicenter(0.5*(currentRmin),0.,0.);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                //set and fill layers
                navilayer->setNextLayer(current);
                fulllayers.push_back(navilayer_pair);
                fulllayers.push_back(currentpair);
            }
            if (currentRmin>0. && currentRmax<Rmax) {
//               std::cout << "Cylinder1Test4" << std::endl;
                bValues.push_back(0.);
                bValues.push_back(currentRmin);
                bValues.push_back(currentRmax);
                bValues.push_back(Rmax);
                //Navilayer
                std::shared_ptr<const Reco::NavigationLayer> navilayer1(new const Reco::NavigationLayer);
                Alg::Point3D navicenter1(0.5*(currentRmin),0.,0.);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer1_pair (navilayer1,navicenter1);
                
                std::shared_ptr<const Reco::NavigationLayer> navilayer2(new const Reco::NavigationLayer);
                Alg::Point3D navicenter2(0.5*(Rmax+currentRmax),0.,0.);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer2_pair (navilayer2,navicenter2);
                //set and fill layers
                navilayer1->setNextLayer(current);
                navilayer2->setPreviousLayer(current);
                fulllayers.push_back(navilayer1_pair);
                fulllayers.push_back(currentpair);
                fulllayers.push_back(navilayer2_pair);
            }
            
        }
    }
    else {
        for (unsigned i=0; i<(layers.size()-1); i++) {
            std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> currentpair = layers.at(i);
            std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> nextpair    = layers.at(i+1);
            std::shared_ptr<const Reco::CylinderLayer> current(std::dynamic_pointer_cast<const Reco::CylinderLayer> (currentpair.first));
            std::shared_ptr<const Reco::CylinderLayer> next(std::dynamic_pointer_cast<const Reco::CylinderLayer> (nextpair.first));
            if (current && next) {
 //               std::cout << "CylinderTest1" << std::endl;
                double currentRmin  = current->getRmin();
                double currentRmax  = current->getRmax();
                double nextRmin     = next->getRmin();
                double nextRmax     = next->getRmax();
            
                //begin
                if (i==0 && currentRmin>0.) {
//                   std::cout << "CylinderTest2" << std::endl;
                    bValues.push_back(0.);
                    //Navilayer
                    std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                    Alg::Point3D navicenter(0.5*(currentRmin),0.,0.);
                    std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                    //set and fill layers
                    navilayer->setNextLayer(current);
                    fulllayers.push_back(navilayer_pair);
                }
                //middle
                if (nextRmin>currentRmax) {
 //                   std::cout << "CylinderTest3" << std::endl;
                    //bins
                    bValues.push_back(currentRmin);
                    bValues.push_back(currentRmax);
                    //Navilayer
                    std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                    Alg::Point3D navicenter(0.5*(nextRmin+currentRmax),0.,0.);
                    std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                    //set and fill layers
                    current->setNextLayer(next);
                    navilayer->setPreviousLayer(current);
                    navilayer->setNextLayer(next);
                    next->setPreviousLayer(current);
                    fulllayers.push_back(currentpair);
                    fulllayers.push_back(navilayer_pair);
                }
                if (nextRmin==currentRmax) {
 //                   std::cout << "CylinderTest4" << std::endl;
                    //bins
                    bValues.push_back(currentRmin);
                    bValues.push_back(currentRmax);
                    //set and fill layers
                    current->setNextLayer(next);
                    next->setPreviousLayer(current);
                    fulllayers.push_back(currentpair);
                }
                //end
                if (i==(layers.size()-2)) {
                    if (nextRmax<Rmax) {
//                        std::cout << "CylinderTest5" << std::endl;
                        //bins
                        bValues.push_back(nextRmin);
                        bValues.push_back(nextRmax);
                        bValues.push_back(Rmax);
                        //Navilayer
                        std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                        Alg::Point3D navicenter(0.5*(Rmax+nextRmax),0.,0.);
                        std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                        //set and fill layers
                        navilayer->setPreviousLayer(next);
                        fulllayers.push_back(nextpair);
                        fulllayers.push_back(navilayer_pair);
                    }
                    if (nextRmax==Rmax) {
//                        std::cout << "CylinderTest6" << std::endl;
                        //bins
                        bValues.push_back(nextRmax-nextRmin);
                        fulllayers.push_back(nextpair);
                    }
                }
            
            }
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::binDiscLayers(LayerVector& layers, LayerVector& fulllayers, std::vector<float>& bValues, Alg::Point3D center, double halfZ) const
{   //for BinendArray
 //   LayerVector fulllayers; //= new LayerVector();
    //counting trough real layers
           //disc
    std::cout << "disclayers.size(): " << layers.size() << std::endl;
    //only one layer
    if (layers.size()==1) {
        std::cout << "layer size 1" << std::endl;
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> currentpair = layers.at(0);
        std::shared_ptr<const Reco::DiscLayer> currentdisc = std::dynamic_pointer_cast<const Reco::DiscLayer> (currentpair.first);
        if (currentdisc) {
 //           std::cout << "Disc1" << std::endl;
            Alg::Point3D* currentcenter = new Alg::Point3D(currentdisc->transform().Translation().Vect());
            double Zmin           = (center.Z())-halfZ;
            double Zmax           = (center.Z())+halfZ;
            double currenthalfz   = currentdisc->getHalfZ();
            double currentzmin    = (currentcenter->Z())-currenthalfz;
            double currentzmax    = (currentcenter->Z())+currenthalfz;
            std::cout << *currentcenter << std::endl;
            std::cout << center << std::endl;
            std::cout << "halfZ: " << halfZ << ", currenthalfz: " << currenthalfz << ", currentzmin: " << currentzmin << ", currentzmax" << currentzmax << ", zmin:" << Zmin << ", zmax:" << Zmax << std::endl;
            if (Zmin==currentzmin && Zmax==currentzmax) {
 //               std::cout << "Disc2" << std::endl;
                bValues.push_back(Zmin);
                bValues.push_back(Zmax);
                fulllayers.push_back(currentpair);
            }
            if (Zmin==currentzmin && currentzmax<Zmax) {
 //               std::cout << "Disc3" << std::endl;
                //bins
                double navithickness = Zmax-currentzmax;
                bValues.push_back(currentzmin);
                bValues.push_back(currentzmax);
                bValues.push_back(Zmax);
                //Navilayer
                Alg::Point3D navicenter(currentcenter->X(),currentcenter->Y(),currentzmax+0.5*navithickness);
                std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                //set and fill layers
                navilayer->setPreviousLayer(currentdisc);
                fulllayers.push_back(currentpair);
                fulllayers.push_back(navilayer_pair);
            }
            if (Zmin<currentzmin && currentzmax==Zmax) {
  //              std::cout << "Disc4" << std::endl;
                //bins
                double navithickness = fabs(Zmin-currentzmin);
                bValues.push_back(Zmin);
                bValues.push_back(currentzmin);
                bValues.push_back(currentzmax);
                //NaviLayer
                Alg::Point3D navicenter(currentcenter->X(),currentcenter->Y(),currentzmin-0.5*navithickness);
                std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                //set and fill layers
                navilayer->setNextLayer(currentdisc);
                fulllayers.push_back(navilayer_pair);
                fulllayers.push_back(currentpair);
            }
            if (Zmin<currentzmin && currentzmax<Zmax) {
 //               std::cout << "Disc5" << std::endl;
                double navithickness1 = fabs(Zmin-currentzmin);
                double navithickness2 = fabs(Zmax-currentzmax);
                //bins
                bValues.push_back(Zmin);
                bValues.push_back(currentzmin);
                bValues.push_back(currentzmax);
                bValues.push_back(Zmax);
                //Navilayers
                Alg::Point3D navi1center(currentcenter->X(),currentcenter->Y(),currentzmin-0.5*navithickness1);
                std::shared_ptr<const Reco::NavigationLayer> navilayer1(new const Reco::NavigationLayer);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer1_pair (navilayer1,navi1center);
                Alg::Point3D navi2center(currentcenter->X(),currentcenter->Y(),currentzmax+0.5*navithickness2);
                std::shared_ptr<const Reco::NavigationLayer> navilayer2(new const Reco::NavigationLayer);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer2_pair (navilayer2,navi2center);
                //set and fill layers
                navilayer1->setNextLayer(currentdisc);
                navilayer2->setPreviousLayer(currentdisc);
                fulllayers.push_back(navilayer1_pair);
                fulllayers.push_back(currentpair);
                fulllayers.push_back(navilayer2_pair);
            }
        }
    }
    else {
        std::cout << "else" << std::endl;
        for (unsigned i=0; i<(layers.size()-1); i++) {
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> currentpair = layers.at(i);
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> nextpair = layers.at(i+1);
        //Disc
        std::shared_ptr<const Reco::DiscLayer> currentdisc = std::dynamic_pointer_cast<const Reco::DiscLayer> (currentpair.first);
        std::shared_ptr<const Reco::DiscLayer> nextdisc = std::dynamic_pointer_cast<const Reco::DiscLayer> (nextpair.first);
        if (currentdisc && nextdisc) {
   //         std::cout << "Test1" << std::endl;
            Alg::Point3D* currentcenter = new Alg::Point3D(currentdisc->transform().Translation().Vect());
            Alg::Point3D* nextcenter    = new Alg::Point3D(nextdisc->transform().Translation().Vect());
            double Zmin         = (center.Z())-halfZ;
            double Zmax         = (center.Z())+halfZ;
            double currenthalfz = currentdisc->getHalfZ();
            double currentzmin  = (currentcenter->Z())-currenthalfz;
            double currentzmax  = (currentcenter->Z())+currenthalfz;
            double nexthalfz    = nextdisc->getHalfZ();
            double nextzmin     = (nextcenter->Z())-nexthalfz;
            double nextzmax     = (nextcenter->Z())+nexthalfz;
            std::cout << "halfZ: " << halfZ << ", currenthalfz: " << currenthalfz << ", currentzmin: " << currentzmin << ", currentzmax" << currentzmax << ", nexthalfz: " << nexthalfz << ", nextzmin: " << nextzmin << ", nextzmax" << nextzmax << ", zmin:" << Zmin << ", zmax:" << Zmax << std::endl;
            //begin
            if (i==0) {
                if (Zmin<currentzmin) {
                    std::cout << "Test2" << std::endl;
                    double navithickness = fabs(Zmin-currentzmin);
                    bValues.push_back(Zmin);
                    //Navilayer
                    Alg::Point3D navicenter(currentcenter->X(),currentcenter->Y(),currentzmin-0.5*navithickness);
                    std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                    std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                    //set and fill layers
                    navilayer->setNextLayer(currentdisc);
                    fulllayers.push_back(navilayer_pair);
                }
            }
            //middle
            if (nextzmin>currentzmax) {
                std::cout << "Test3" << std::endl;
                //bins
                double navithickness = fabs(nextzmin-currentzmax);
                bValues.push_back(currentzmin);
                bValues.push_back(currentzmax);
                //Navilayer
                Alg::Point3D navicenter(currentcenter->X(),currentcenter->Y(),currentzmax+0.5*navithickness);
                std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                //set and fill layers
                currentdisc->setNextLayer(nextdisc);
                navilayer->setPreviousLayer(currentdisc);
                navilayer->setNextLayer(nextdisc);
                nextdisc->setPreviousLayer(currentdisc);
                fulllayers.push_back(currentpair);
                fulllayers.push_back(navilayer_pair);
            }
            if (nextzmin==currentzmax) {
 //               std::cout << "Test4" << std::endl;
                //bins
                bValues.push_back(currentzmin);
                bValues.push_back(currentzmax);
                //set and fill layers
                currentdisc->setNextLayer(nextdisc);
                nextdisc->setPreviousLayer(currentdisc);
                fulllayers.push_back(currentpair);
            }
            //end
            if (i==(layers.size()-2)) {
                if (nextzmax<Zmax) {
//                    std::cout << "Test5" << std::endl;
                    //bins
                    double navithickness = fabs(Zmax-nextzmax);
                    bValues.push_back(nextzmin);
                    bValues.push_back(nextzmax);
                    bValues.push_back(Zmax);
                    //Navilayer
                    Alg::Point3D navicenter(nextcenter->X(),nextcenter->Y(),nextzmax+0.5*navithickness);
                    std::shared_ptr<const Reco::NavigationLayer> navilayer(new const Reco::NavigationLayer);
                    std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D> navilayer_pair (navilayer,navicenter);
                    //set and fill layers
                    navilayer->setPreviousLayer(nextdisc);
                    fulllayers.push_back(nextpair);
                    fulllayers.push_back(navilayer_pair);
                }
                if (nextzmax==Zmax) {
//                    std::cout << "Test6" << std::endl;
                    bValues.push_back(nextzmin);
                    bValues.push_back(nextzmax);
                    fulllayers.push_back(nextpair);
                }
            }
          
        }
    }
    }

    return StatusCode::SUCCESS;
}
