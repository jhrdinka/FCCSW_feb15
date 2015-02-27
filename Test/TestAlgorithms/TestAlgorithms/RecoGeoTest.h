//
//  RecoGeoTest.h
//  
//
//  Created by Julia Hrdinka on 26/02/15.
//
//

#ifndef RECOGEOTEST_H
#define RECOGEOTEST_H

#include "GaudiAlg/GaudiAlgorithm.h"
#include "DetDesInterfaces/IRecoGeoSvc.h"
//RecoGeometry
#include "RecoGeometry/ContainerVolume.h"

class RecoGeoTest: public GaudiAlgorithm {
    friend class AlgFactory<RecoGeoTest> ;
    
public:
    /// Constructor.
    RecoGeoTest(const std::string& name, ISvcLocator* svcLoc);
    /// Initialize.
    virtual StatusCode initialize();
    /// Execute.
    virtual StatusCode execute();
    /// Finalize.
    virtual StatusCode finalize();
    
private:
    
    IRecoGeoSvc* m_recoGeoSvc; //Detector Description Service
};


#endif //RECOGEOTEST_H
