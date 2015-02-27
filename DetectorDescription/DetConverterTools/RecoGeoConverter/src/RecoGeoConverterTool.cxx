//
//  RecoGeoConverter.cxx
//  
//
//  Created by Julia Hrdinka on 10/11/14.
//
//

#include "RecoGeoConverter/RecoGeoConverterTool.h"

DECLARE_COMPONENT(RecoGeoConverterTool)

// Standard Constructor
RecoGeoConverterTool::RecoGeoConverterTool(const std::string& type,
                                   const std::string& name,
                                   const IInterface* parent)
    : AlgTool( type, name, parent ),
    m_log(msgSvc(), name)
{
        m_out.open("detelements.txt");
        
        
// Declare additional interface
    declareInterface<IGeoConverterTool>(this);

//Declare Properties of the specific tool
        //=> for JobOptions
    if (service("DD4HepDetDesSvc", m_DD4HepSvc).isFailure()) m_log << MSG::ERROR << "Couldn't get DD4HepDetDesSvc" << endmsg;
    if (service("RecoGeoSvc", m_recoGeoSvc, true).isFailure()) m_log << "Couldn't get RecoGeoSvc" << endmsg;
}

double RecoGeoConverterTool::flatrand(double min, double max) {
    std::uniform_real_distribution<double> unif(min,max);
    std::random_device rd;
    std::mt19937 gen(rd());
    
    return (unif(gen));
}

StatusCode RecoGeoConverterTool::convert() {
    
    DD4hep::Geometry::DetElement detworld = m_DD4HepSvc->lcdd()->world();
            translateDetector(detworld);
    m_counter=0;
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

    //multimap of the volumes, ordered by their hierarchy
    std::multimap<int, std::shared_ptr<const Reco::Volume>> volumes;
    std::vector<std::shared_ptr<const Reco::Volume>> previousvol(4,0);
    std::vector<std::shared_ptr<const Reco::Volume>> currentvol(4,0);
    std::vector<float> bZValues(Reco::ContainerCylinderVolume::zLength);
    std::vector<float> bRValues(Reco::ContainerCylinderVolume::rLength);
    std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>* binRvector = new std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>(2);
    std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>* binZvector = new std::vector<std::pair<std::shared_ptr<const Reco::Volume>, Alg::Point3D>>(3);

    //translate the volumes
    translateVolume(detelement, volumes);
    
    //putting keys (hierarchy) in a vector;
    std::vector<std::pair<int,std::shared_ptr<const Reco::Volume>>> keys_dedup;
    std::unique_copy(begin(volumes),
                end(volumes),
                back_inserter(keys_dedup),
                     [](const std::pair<int,std::shared_ptr<const Reco::Volume>> &entry1,
                        const std::pair<int,std::shared_ptr<const Reco::Volume>> &entry2) {
                    return (entry1.first == entry2.first);}
                     );
    m_out << "Number of Hierarchies: " << keys_dedup.size() << std::endl;
    int last_key = keys_dedup.size()-1;
    for (const auto& itkeys : keys_dedup) {
        std::pair <std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator, std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator> ret;
        ret = volumes.equal_range(itkeys.first);
        for (std::multimap<int, std::shared_ptr<const Reco::Volume>>::iterator itvol=ret.first; itvol!=ret.second; ++itvol) {
            std::shared_ptr<const Reco::Volume> volume_ptr(itvol->second);
            currentvol.at(itvol->second->type()) = itvol->second;
        } //for volumevectors fill
        
       //now group the volumes together in containers
       if (itkeys.first==0) {
            previousvol = currentvol;
        }
        
        else if (itkeys.first==1) {
            //falls ich mehr volumes mache, unterscheidung in volume type
            
            //set Volumes for BoundarySurfaces between the volumes
            currentvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
              //beamtube
            if (!previousvol.empty()) {
            currentvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            }
            //Konvention: volumes muessen genau aufeinader treffen, was wenn nicht?->noch ueberpruefungen machen und abfangen
            Alg::Point3D center1(0.,0.,0.);
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
            //barrel
            center1 = currentvol.at(Reco::Volume::barrel)->center();
            halfZ   = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nBarrel) = (center1.Z()-halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pBarrel) = (center1.Z()+halfZ);
            binZvector->at(1) = (make_pair(currentvol.at(Reco::Volume::barrel),center1));
            //pos EndCap
            center2 = currentvol.at(Reco::Volume::pEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::pEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pEnd) = (center2.Z()+halfZ);
            binZvector->at(2) = (make_pair(currentvol.at(Reco::Volume::pEndCap),center2));
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
            m_out << "created containervolume" << std::endl;
            //beampipe
            contvol->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolume(previousvol.at(Reco::Volume::barrel));
            //rest
            contvol->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedvolumes));
            contvol->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedvolumes));
            contvol->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::nEndCap));
            contvol->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(currentvol.at(Reco::Volume::pEndCap));
            contvol->setType(Reco::Volume::container);
            currentvol.at(Reco::Volume::container) = contvol;
            }
            previousvol = currentvol;
            //add containervolume in the end
        } //if status==1
        else {
        //    if(!previousvol.at(Reco::Volume::container)) std::cout << "no container " << std::endl; auch ueberpruefungen machen
            //set Volumes for BoundarySurfaces
        //previous container volume
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
        //previous volumes
                //barrel
                previousvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                previousvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                previousvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setNextVolume(currentvol.at(Reco::Volume::barrel));
                //EndCaps
                previousvol.at(Reco::Volume::nEndCap)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
                previousvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
            
        //current volumes
            //barrel
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setPreviousVolumes(previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->getPreviousVolumes());
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(currentvol.at(Reco::Volume::nEndCap));
            currentvol.at(Reco::Volume::barrel)->getBoundarySurface(Reco::CylinderVolume::posDisc)->setNextVolume(currentvol.at(Reco::Volume::pEndCap));
            //EndCaps
             //Grenzen muessen genau ueberein stimmen!
             //make binned Array
            //binutility
            bRValues.at(Reco::ContainerCylinderVolume::inner) = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmin);
            bRValues.at(Reco::ContainerCylinderVolume::middle) = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmax);
            bRValues.at(Reco::ContainerCylinderVolume::outer) = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmax);
            Trk::BinUtility* binutil = new Trk::BinUtility(bRValues,Trk::open,Trk::binR);
            Alg::Point3D center1(0.,0.,0.);
            Alg::Point3D center2(0.,0.,0.);
            double Rmin = 0;
            double Rmax = 0;
            Rmin = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmin);
            Rmax = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::Rmax);
            center1.SetCoordinates(0.5*(Rmin+Rmax),0.,0.);
            binRvector->at(1)  = (make_pair(currentvol.at(Reco::Volume::barrel),center1));
            //negativ
            Rmin = previousvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::Rmin);
            Rmax = previousvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::Rmax);
            center1.SetCoordinates(0.5*(Rmin+Rmax),0.,0.);
            binRvector->at(0) = (make_pair(previousvol.at(Reco::Volume::nEndCap),center1));
            Trk::BinnedArray1D<Reco::Volume>* binnedRvolumes1 = new Trk::BinnedArray1D<Reco::Volume>(*binRvector,binutil);
            //positiv
            Rmin = previousvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::Rmin);
            Rmax = previousvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::Rmax);
            center1.SetCoordinates(0.5*(Rmin+Rmax),0.,0.);
            binRvector->at(0) = (make_pair(previousvol.at(Reco::Volume::pEndCap),center1));
            Trk::BinnedArray1D<Reco::Volume>* binnedRvolumes2 = new Trk::BinnedArray1D<Reco::Volume>(*binRvector,binutil);
            currentvol.at(Reco::Volume::pEndCap)->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes2));
        //make ContainerVolume1 //momentan wird containervolume an containervolume 1 uebergeben und nicht schon aufgeloest in untervolumen -gut? ABER FUER NEext und previousvolumes is schon aufegloest, also passt
            Rmin = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmin);
            Rmax = previousvol.at(Reco::Volume::container)->getCoordinate(Reco::CylinderVolume::Rmax);
            center1.SetCoordinates(0.5*(Rmin+Rmax),0.,0.);
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
            //set surfacepointers of ContainerVolume1
            contvol1->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->getNextVolumes()); //siehe unten
            contvol1->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(currentvol.at(Reco::Volume::barrel));
            contvol1->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes1));
            contvol1->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*binnedRvolumes2));
            //beamtube
            //make ContainerVolume2
            //fill bZValues
            //neg EndCap
            center1 = currentvol.at(Reco::Volume::nEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::nEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nEnd) = (center1.Z()-halfZ);
            binZvector->at(0) = (make_pair(currentvol.at(Reco::Volume::nEndCap),center1));
            //barrel
            center1 = contvol1->center();
            halfZ   = currentvol.at(Reco::Volume::barrel)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::nBarrel) = (center1.Z()-halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pBarrel) = (center1.Z()+halfZ);
            binZvector->at(1) = (make_pair(std::shared_ptr<const Reco::Volume>(contvol1),center1));
            //pos EndCap
            center2 = currentvol.at(Reco::Volume::pEndCap)->center();
            halfZ   = currentvol.at(Reco::Volume::pEndCap)->getCoordinate(Reco::CylinderVolume::halfZ);
            bZValues.at(Reco::ContainerCylinderVolume::pEnd) = (center2.Z()+halfZ);
            binZvector->at(2) = (make_pair(currentvol.at(Reco::Volume::pEndCap),center2));
            
            Trk::BinUtility* binZutil = new Trk::BinUtility(bZValues, Trk::open, Trk::binZ);
            conthalfZ = sqrt((center2.Z()-center1.Z())*(center2.Z()-center1.Z()))+halfZ;
            Trk::BinnedArray1D<Reco::Volume>* containerbin2 = new Trk::BinnedArray1D<Reco::Volume>(*binZvector,binZutil);
            std::shared_ptr<const Reco::ContainerCylinderVolume> contvol2(new Reco::ContainerCylinderVolume(containerbin2,conttrans,contRmin,contRmax,conthalfZ));
            //make new binnedarray
            center1 = currentvol.at(Reco::Volume::barrel)->center();
            binZvector->at(1) = (make_pair(std::shared_ptr<const Reco::Volume>(currentvol.at(Reco::Volume::barrel)),center1));
            Trk::BinnedArray1D<Reco::Volume>* containerbin3 = new Trk::BinnedArray1D<Reco::Volume>(*binZvector,binZutil);
            //set surfacepointers of ContainerVolume2
            const Trk::BinnedArray1D<Reco::Volume>* containerbin4 = previousvol.at(Reco::Volume::container)->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->getNextVolumes();
            contvol2->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolumes(new Trk::BinnedArray1D<Reco::Volume>(*containerbin4));
            contvol2->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolumes(new Trk::BinnedArray1D<Reco::Volume>(*containerbin3));
            contvol2->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(currentvol.at(Reco::Volume::nEndCap));
            contvol2->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(currentvol.at(Reco::Volume::pEndCap));

            previousvol = currentvol;
            
            //hand worldvolume over to the RecoGeoSvc
            if (itkeys.first==last_key) {
                std::cout << "set WorldVolume!!!!" << std::endl;
                m_recoGeoSvc->setWorldVolume(contvol2);
            }
            
        }
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
        m_out << "detcylindervolume" << std::endl;
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
                            m_out << "size fulllayers: " << fulllayers->size() << std::endl;
                            Trk::BinUtility* binutility = new Trk::BinUtility(bValues,Trk::open,Trk::binR);
                            if(binutility) {
                                Trk::BinnedArray1D<Reco::Layer>* binnedlayers = new Trk::BinnedArray1D<Reco::Layer>(*fulllayers,binutility);
                                std::shared_ptr<const Reco::Volume> cylindervolume(new Reco::CylinderVolume(binnedlayers, geonode, tube));
                                if (cylindervolume) {
                                    m_out  << "##### created cylindervolume" << std::endl;
                                    //set the volumes of the boundary surfaces
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }

                                    cylindervolume->setType(Reco::Volume::barrel);
                                    m_out << "status" << status << std::endl;
                                    volumes.emplace(status,cylindervolume);
                                } //if cylindervolume
                            }
                            else return StatusCode::FAILURE;
                        }
                    }//if !layer.empty()
                    else {
                        std::shared_ptr<const Reco::Volume> cylindervolume(new Reco::CylinderVolume(geonode, tube));
                                if (cylindervolume) {
                                    m_out  << "##### created cylindervolume without layers" << std::endl;
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }
                                    cylindervolume->setType(Reco::Volume::barrel);
                                    m_out << "status" << status << std::endl;
                                    volumes.emplace(status,cylindervolume);
                                } //if cylindervolume
                            else return StatusCode::FAILURE;
                    } // layer.empty()
                }//if tube
            }//if geoshape
        } //if geonode
        
    } //if detcylindervolume
    
    if (detdiscvolume) {
        m_out << "detdiscvolume" << std::endl;
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
            m_out << *transform << std::endl;
            if(geoshape && geoshape->IsA() == TGeoConeSeg::Class()){
                TGeoConeSeg* tube = dynamic_cast<TGeoConeSeg*>(geoshape);
                if(tube){
                    LayerVector layers;          //translateLayers
                    m_counter = 0;
                    translateLayer(detelement,layers, transform);
                    if(!layers.empty()){
                        double halfZ = tube->GetDz();
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
                            m_out << "size fulllayers: " << fulllayers->size() << std::endl;
                            m_out << "bValues size: " << bValues.size() << std::endl;
                            Trk::BinUtility* binutility = new Trk::BinUtility(bValues,Trk::open,Trk::binZ);
                            m_out << "binutility size: " << binutility->dimensions() << std::endl;;
                            if(binutility) {
                                Trk::BinnedArray1D<Reco::Layer>* binnedlayers = new Trk::BinnedArray1D<Reco::Layer>(*fulllayers,binutility);
                                std::shared_ptr<const Reco::CylinderVolume> cylindervolume(new Reco::CylinderVolume(binnedlayers, geonode, tube));
                                if (cylindervolume) {

                                    m_out << cylindervolume->transform() << std::endl;
                                    
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::posDisc)->setPreviousVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::negDisc)->setNextVolume(cylindervolume);
                                    cylindervolume->getBoundarySurface(Reco::CylinderVolume::outerCylinder)->setPreviousVolume(cylindervolume);
                                    if (cylindervolume->NumberOfSurfaces()==4) {
                                        cylindervolume->getBoundarySurface(Reco::CylinderVolume::innerCylinder)->setNextVolume(cylindervolume);
                                    }
                                    if (center.Z()<0.) cylindervolume->setType(Reco::Volume::nEndCap);
                                    else cylindervolume->setType(Reco::Volume::pEndCap);
                                    m_out << "status" << status << std::endl;
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
            m_out << "detcylinderlayer" << std::endl;
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
                            m_out << "surfaces empty" << std::endl;
                        }
                        if (!surfaces.empty()) {
                            Trk::BinnedArray2D<Reco::Surface>* bin = new Trk::BinnedArray2D<Reco::Surface>(surfaces,currentBinUtility);
                            std::shared_ptr<const Reco::CylinderLayer> cyllayer(new Reco::CylinderLayer(transf,tube,bin));
                            if (cyllayer) {
                                m_out << "### Created CylinderLayer ###" << std::endl;

                                double R = 0.5*(cyllayer->getRmax()+cyllayer->getRmin());
                                Alg::Point3D center(R,0.,0.);
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
            m_out << "detdisclayer" << std::endl;
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
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::translateModule(DD4hep::Geometry::DetElement det, std::vector<std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D>>& surfaces, std::shared_ptr<const Alg::Transform3D> transform)
{
    const  DD4hep::Geometry::DetElement::Children& children = det.children();
    int m = 0;
    for (DD4hep::Geometry::DetElement::Children::const_iterator i=children.begin(); i != children.end(); ++i)
    {
        m++;
        DD4hep::Geometry::PlacedVolume pv = (*i).second.placement();
        Det::IExtension* ext = (*i).second.extension<Det::IExtension>();
        Det::DetModule* detm = dynamic_cast<Det::DetModule*>(ext);
        DD4hep::Geometry::DetElement detelement = (*i).second;
        //MODULE
        if  (detm) {
              m_out << "Module" << std::endl;
              m_out << "Child from Layer: " << m_counter << ", " << (*i).first << std::endl;
            TGeoNode* node = pv.ptr();
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
                        std::map<std::pair<int,int>, Reco::Material*>* map = new std::map<std::pair<int,int>, Reco::Material*>();
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
  
                                        TGeoShape* childshape = childnode->GetVolume()->GetShape();
                                        
                                        if (childshape && childshape->IsA()==TGeoBBox::Class()) {
                                            TGeoBBox* childbox = dynamic_cast<TGeoBBox*>(childshape);
                                            if ((fabs(local[0])<=(childbox->GetDX())) && (fabs(local[1])<=(childbox->GetDY()))) {
                                                t          = 2*(childbox->GetDZ());
                                                sumt      += t;
                                                t_x0      += t/mat.radLength();
                                                t_lambda0 += t/mat.intLength();
                                                A         += mat.density()*mat.A();
                                                Z         += mat.density()*mat.Z();
                                                density   += mat.density()*t;
                                                sumdens   += mat.density();
                                            }
                                        }
                                        
                                    }
                                    
                                    if (m==0 && n==0) m_out << " component" << comp_num << " t: " << t << "sumt: " << sumt << std::endl;
                                    ++comp_num;
                                } //for components
                                
                                sensper = sensper/capacity;
                                density = density/sumt;
                                A       = A/sumdens;
                                Z       = Z/sumdens;
                                //fill map with material for each bin
                                map->emplace(std::make_pair(m,n), new Reco::Material(A, Z, density, t_x0, t_lambda0,sensper));
                            } // for binsY
                        } //for binsX
                        //create MaterialMap
                        Reco::MaterialMap* materialmap = new Reco::MaterialMap(binutility, map);
                        //create Surface
                        std::shared_ptr<const Reco::PlaneSurface> plane(new Reco::PlaneSurface(box,materialmap,transf));
                        if(plane) {
                            m_out << "#created plane#" << std::endl;
                            Alg::Point3D center = plane->center();
                            std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D> surface (plane,center);
                            surfaces.push_back(surface);
                        }
                    } //if box
                    
                } //ifshape && isA Box
                
                if(shape && shape->IsA() == TGeoTrd2::Class())
                {
                    TGeoTrd2* trapezoid = dynamic_cast<TGeoTrd2*>(shape);
                    if (trapezoid) {
                        double halfY   = trapezoid->GetDz();
                        double halfX1  = trapezoid->GetDx1();
                        double halfX2  = trapezoid->GetDx2();
                        int binsX = 100;
                        int binsY = 100;
                        Trk::BinUtility* binutility  = new Trk::BinUtility(binsY,-halfY,halfY,Trk::open,Trk::binY);
                        std::vector<Trk::BinUtility*>* binXvector = new std::vector<Trk::BinUtility*>(binsY);
                        std::vector<int> binsXvec(binsY);
                        double k = (2.*halfY)/(halfX2-halfX1);
                        double d = -k*(halfX2+halfX1)*0.5;
                        double stepsizeY = (2.*halfY)/binsY;
                        double x = 0;
                        double y = 0;
                        for (int i=0; i<binsY; i++) {
                            y = -halfY + i*stepsizeY;
                            x = (y-k)/d;
                            Trk::BinUtility* binXutil = new Trk::BinUtility(binsX,-x,x,Trk::open,Trk::binX);
                            binXvector->at(i) = binXutil;
                        }
                        //map for the material
                        std::map<std::pair<int,int>, Reco::Material*>* map = new std::map<std::pair<int,int>, Reco::Material*>();
                        double newx = 0.;
                        double newy = 0.;
                        for (int m=0; m<binsY; m++)
                        {
                            newy = binutility->bincenter(m,0);
                            for (int n=0; n<binsX; n++)
                            {
                                //capacity of all module components
                                double capacity = 0;
                                //senitive percentage
                                double sensper = 0;
                                newx = binXvector->at(m)->bincenter(n,0);
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
 
                                    DD4hep::Geometry::Volume vol   = (*j).second.volume();
                                    DD4hep::Geometry::Material mat = vol.material();
                                    //calculate the whole volume of all components
                                    capacity += vol.ptr()->Capacity();
                                    //calculate the sensitive volume
                                    if(vol.isSensitive()) sensper += vol.ptr()->Capacity();
                                    TGeoNode* childnode = (*j).second.placement().ptr();
                                    if (childnode) {
                                        childnode->MasterToLocal(master,local);
                                        
                                        TGeoShape* childshape = childnode->GetVolume()->GetShape();
                                        
                                        if (childshape && childshape->IsA()==TGeoTrd2::Class()) {
                                            TGeoTrd2* childtrapezoid = dynamic_cast<TGeoTrd2*>(childshape);
                                            if (isInsideTrapezoid(local,childtrapezoid->GetDx1(),childtrapezoid->GetDx2(),childtrapezoid->GetDz())) {
                                                t          = 2*(childtrapezoid->GetDY());
                                                sumt      += t;
                                                t_x0      += t/mat.radLength();
                                                t_lambda0 += t/mat.intLength();
                                                A         += mat.density()*mat.A();
                                                Z         += mat.density()*mat.Z();
                                                density   += mat.density()*t;
                                                sumdens   += mat.density();
                                            }
                                        }
                                        
                                    }
                                    ++comp_num;
                                } //for components
                                
                                sensper = sensper/capacity;
                                density = density/sumt;
                                A       = A/sumdens;
                                Z       = Z/sumdens;
                                //fill map with material for each bin
                                map->emplace(std::make_pair(m,n), new Reco::Material(A, Z, density, t_x0, t_lambda0,sensper));
         
                            } // for binsY
                        } //for binsX
                        //create MaterialMap
                        Reco::MaterialMap* materialmap = new Reco::MaterialMap(binutility, binXvector, map);
                        //create Surface
                        std::shared_ptr<const Reco::TrapezoidSurface> trapez(new Reco::TrapezoidSurface(trapezoid,materialmap,transf));
                        if(trapez) {
                            m_out << "#created trapezoid#" << std::endl;
                            Alg::Point3D center = trapez->center();
                            std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D> surface (trapez,center);
                            surfaces.push_back(surface);
                        }
                    } //if trapezoid
                } //ifshape && isA Trapezoid
                
            } //ifnode
            
            //Histograms

            
        } //if detector module
    
    } //for layer children
    return StatusCode::SUCCESS;
}

StatusCode RecoGeoConverterTool::binCylinderLayers(LayerVector& layers, LayerVector& fulllayers, std::vector<float>& bValues, Alg::Point3D, double Rmax) const
{
 //   m_out << "cylinderlayers.size(): " << layers.size() << std::endl;
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
        sort(layers.begin(),layers.end(),sortCylinderLayers);
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
{
    //disc
//    m_out << "disclayers.size(): " << layers.size() << std::endl;
    //only one layer
    if (layers.size()==1) {
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
//            std::cout << "halfZ: " << halfZ << ", currenthalfz: " << currenthalfz << ", currentzmin: " << currentzmin << ", currentzmax" << currentzmax << ", zmin:" << Zmin << ", zmax:" << Zmax << std::endl;
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
 //       std::cout << "else" << std::endl;
        sort(layers.begin(),layers.end(),sortDiscLayers);
        for (unsigned i=0; i<(layers.size()-1); i++) {
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> currentpair = layers.at(i);
        std::pair<std::shared_ptr<const Reco::Layer>,Alg::Point3D> nextpair = layers.at(i+1);
        //Disc
 //           std::cout << "Point3D currentpair: " << currentpair.second << std::endl;
 //           std::cout << "Point3D nextpair: " << nextpair.second << std::endl;
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
 //           std::cout << "halfZ: " << halfZ << ", currenthalfz: " << currenthalfz << ", currentzmin: " << currentzmin << ", currentzmax" << currentzmax << ", nexthalfz: " << nexthalfz << ", nextzmin: " << nextzmin << ", nextzmax" << nextzmax << ", zmin:" << Zmin << ", zmax:" << Zmax << std::endl;
            //begin
            if (i==0) {
                if (Zmin<currentzmin) {
  //                  std::cout << "Test2" << std::endl;
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
    //            std::cout << "Test3" << std::endl;
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

bool RecoGeoConverterTool::isInsideTrapezoid(double locpos[3], double trapX1, double trapX2, double trapY) const
{
    double xmin = 0;
    double xmax = 0;
    if (trapX2>trapX1) {
        xmin = trapX1;
        xmax = trapX2;
    }
    else {
        xmin = trapX2;
        xmax = trapX1;
    }
    if (fabs(locpos[0])>xmax) return false;
    if (fabs(locpos[1])>trapY) return false;
    if (fabs(locpos[0])<xmin) return true;
    if (trapX1==trapX2) return true;
    
    double sign = ((locpos[0]>0.) ? 1. : -1. );
    double k = sign*(2.*trapY)/(trapX2-trapX1);
    double d = -k*(trapX2+trapX1)*0.5;
    double x = (locpos[1]-d)/k;
    if (fabs(locpos[0])<fabs(x)) return true;
    else return false;
}

bool RecoGeoConverterTool::sortCylinderLayers(const std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D>& a,const std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D>& b)
{
    return (a.second.X()<b.second.X());
}

bool RecoGeoConverterTool::sortDiscLayers(const std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D>& a,const std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D>& b)
{
    return (a.second.Z()<b.second.Z());
}








