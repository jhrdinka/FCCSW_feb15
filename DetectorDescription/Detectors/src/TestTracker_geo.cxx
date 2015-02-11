//
//  TestTracker_geo.cc
//  
//
//  Created by Julia Hrdinka on 21/10/14.
//
//

/*#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/TGeoUnits.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;


static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {

    xml_det_t   x_det = e;
    string      name = x_det.nameStr();
    DetElement  tracker(name, x_det.id());
    
    Volume experimentalHall_log =  lcdd.pickMotherVolume(tracker);
    
    double radius   = x_det.radius();
    
    xml_comp_t TrkElement (x_det.child(_U(box)));
    
    double deltaphi = TrkElement.deltaphi();
    
    Box TrkBox(TrkElement.width(),TrkElement.thickness(),TrkElement.length());
    Volume TrkBoxVolume("TrkBoxVolume", TrkBox, lcdd.material(TrkElement.materialStr()));
    
    //TrkBoxVolume.setSensitiveDetector(sens);
    
    for (int i = 0, n = TrkElement.repeat(); i < n; ++i){
        
        double phi = deltaphi/dd4hep::rad * i;
        
        Position trans(radius * sin(phi),
                       radius * cos(phi),
                       0.);
        
        PlacedVolume placedBox = experimentalHall_log.placeVolume(TrkBoxVolume, Transform3D(RotationZ(-phi-0.25*M_PI), trans));
        
         placedBox.addPhysVolID("system",x_det.id());
         placedBox.addPhysVolID("module", n);
        
        
        tracker.setPlacement(placedBox);
        
    }

   
    tracker.setVisAttributes(lcdd, x_det.visStr(), TrkBoxVolume);
    
    
    return tracker;
}

DECLARE_DETELEMENT( TestTracker, create_element )*/


#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/TGeoUnits.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;


static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)  {
    
    xml_det_t   x_det = e;
    string      name = x_det.nameStr();
    DetElement  tracker(name, x_det.id());
    
    Volume experimentalHall_log =  lcdd.pickMotherVolume(tracker);
    
    xml_comp_t TrkElement (x_det.child(_U(box)));
    int counter = 0;
    
    for (xml_coll_t j(e,_U(layer)); j; ++j ) {
    xml_comp_t x_layer = j;
        ++counter;
    
    double radius = x_layer.radius();
    
    double deltaphi = TrkElement.deltaphi();
    
    Box TrkBox(TrkElement.width(),TrkElement.thickness(),TrkElement.length());
    Volume TrkBoxVolume("TrkBoxVolume", TrkBox, lcdd.material(x_layer.materialStr()));
    
    //TrkBoxVolume.setSensitiveDetector(sens);
    
        double z = x_layer.dz();
        //double z = 10.*counter;
        
    for (int i = 0, n = TrkElement.repeat(); i < n; ++i){
        
        int layerID = 0;
        
        double phi = deltaphi/dd4hep::rad * i;
        
        Position trans(radius * sin(phi),
                       radius * cos(phi),
                       z);
        
        PlacedVolume placedBox = experimentalHall_log.placeVolume(TrkBoxVolume, Transform3D(RotationZ(-phi-0.11*M_PI), trans));
        
        placedBox.addPhysVolID("system",x_det.id());
        placedBox.addPhysVolID("layer",layerID);
        placedBox.addPhysVolID("module", n);
        
        
        
        tracker.setPlacement(placedBox);
        
    }
    
    tracker.setVisAttributes(lcdd, x_layer.visStr(), TrkBoxVolume);
    
  }
    return tracker;
}

DECLARE_DETELEMENT( TestTracker1, create_element )
