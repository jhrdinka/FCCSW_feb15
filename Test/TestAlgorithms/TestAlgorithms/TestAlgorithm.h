#ifndef _TESTALGORITHM_H_
#define _TESTALGORITHM_H_

#include "GaudiAlg/GaudiAlgorithm.h"
#include "DetDesInterfaces/IDetDesSvc.h"

class TestAlgorithm: public GaudiAlgorithm {
  friend class AlgFactory<TestAlgorithm> ;

public:
  /// Constructor.
  TestAlgorithm(const std::string& name, ISvcLocator* svcLoc);
  /// Initialize.
  virtual StatusCode initialize();
  /// Execute.
  virtual StatusCode execute();
  /// Finalize.
  virtual StatusCode finalize();

private:
    
    IDetDesSvc* m_detsvc; //Detector Description Service
};

#endif
