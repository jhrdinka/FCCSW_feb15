//
//  RecoGeoConverterTool.h
//  
//
//  Created by Julia Hrdinka on 10/11/14.
//
//

#ifndef RECO_RECOGEOCONVERTERTOOL_H
#define RECO_RECOGEOCONVERTERTOOL_H

#define NMAX 100 //spaeter in fkt mitgeben lassen
//Gaudi
#include "GaudiKernel/AlgTool.h"
#include "SimulationInterfaces/IGeoConverterTool.h"
#include "GaudiKernel/ITHistSvc.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "GaudiKernel/DataObjectHandle.h"
//DD4hep
#include "DD4hep/LCDD.h"
#include "DD4hep/Detector.h"
//Root
#include "TGeoManager.h"
//RecoGeometry
#include "RecoGeometry/Surface.h"
#include "RecoGeometry/PlaneSurface.h"
#include "RecoGeometry/CylinderSurface.h"
#include "RecoGeometry/Layer.h"
#include "RecoGeometry/CylinderLayer.h"
#include "RecoGeometry/Volume.h"
#include "RecoGeometry/BoundarySurface.h"
#include "RecoGeometry/CylinderVolume.h"
//std
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <random>

class RecoGeoConverterTool : public AlgTool, virtual public IGeoConverterTool {
public:
    
    typedef std::vector<std::pair<std::shared_ptr<const Reco::Layer>, Alg::Point3D>> LayerVector;
    
    RecoGeoConverterTool(const std::string& type, const std::string& name, const IInterface* parent);
    virtual ~RecoGeoConverterTool (){
        m_out.close();
    }
    
    
  //  virtual StatusCode initialize();
  //  virtual StatusCode finalize ();
    
    double flatrand(double min, double max);
    
    virtual StatusCode convert(DD4hep::Geometry::LCDD* lcdd);
    
    StatusCode scanDetector(DD4hep::Geometry::DetElement detelement);
    
    StatusCode translateDetector(DD4hep::Geometry::DetElement det);
    
    StatusCode translateVolume(DD4hep::Geometry::DetElement det, std::multimap<int, std::shared_ptr<const Reco::Volume>>& volumes);
    
    StatusCode translateModule(DD4hep::Geometry::DetElement det, std::vector<std::pair<std::shared_ptr<const Reco::Surface>, Alg::Point3D>>& surfaces, std::shared_ptr<const Alg::Transform3D> transform);
    
    StatusCode translateLayer(DD4hep::Geometry::DetElement det, LayerVector& layers, std::shared_ptr<const Alg::Transform3D> transform);
    StatusCode binCylinderLayers(LayerVector& layers, LayerVector& fulllayers, std::vector<float>& bValues, Alg::Point3D center, double Rmax) const;
    StatusCode binDiscLayers(LayerVector& layers, LayerVector& fulllayers, std::vector<float>& bValues, Alg::Point3D center, double halfZ) const;
 //   StatusCode binLayers(std::pair<std::shared_ptr<const Reco::Layer>,Alg::Vector3D> currentpair, std::pair<std::shared_ptr<const Reco::Layer>,Alg::Vector3D> nextpair, const Reco::Layer current, const Reco::Layer next) const;
//    Reco::Surface& installSurface(TGeoNode* node);
    

    

    
private:
    MsgStream               m_log;
    DD4hep::Geometry::LCDD* m_lcdd;
    std::ofstream           m_out;
    int                     m_counter;
    //Histogram
 /*   ITHistSvc*              m_ths;
    TH2F*                   m_t;
    TH2F*                   m_sens;
    TH2F*                   m_dens;
    TH2F*                   m_tx0;
    TH2F*                   m_tlambda0;
    TH2F*                   m_A;
    TH2F*                   m_Z;
  
*/
};

#endif //RECO_RECOGEOCONVERTERTOOL_H
