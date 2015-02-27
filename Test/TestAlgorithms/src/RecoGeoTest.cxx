//
//  RecoGeoTest.cxx
//  
//
//  Created by Julia Hrdinka on 26/02/15.
//
//

#include "TestAlgorithms/RecoGeoTest.h"


DECLARE_COMPONENT(RecoGeoTest)

RecoGeoTest::RecoGeoTest(const std::string& name, ISvcLocator* svcLoc):
GaudiAlgorithm(name, svcLoc),
m_recoGeoSvc(nullptr)
{}

StatusCode RecoGeoTest::initialize() {
    
    if (GaudiAlgorithm::initialize().isFailure()){
        return StatusCode::FAILURE;
    }
    
    if (service("RecoGeoSvc", m_recoGeoSvc, true).isFailure()) {
        error() << "Couldn't get RecoGeoSvc" << endmsg;
        return StatusCode::FAILURE;
    }
    
    return StatusCode::SUCCESS;
}


StatusCode RecoGeoTest::execute() {
    
    m_recoGeoSvc->buildGeometry();
    std::shared_ptr<const Reco::ContainerVolume> worldVolume(m_recoGeoSvc->getWorldVolume()->clone());
    
    if (worldVolume) {
        std::cout << "retrieved WorldVolume!!!" << std::endl;
    }
    return StatusCode::SUCCESS;
}


StatusCode RecoGeoTest::finalize() {
    
    if (GaudiAlgorithm::finalize().isFailure())
        return StatusCode::FAILURE;
    
    return StatusCode::SUCCESS;
    
}
