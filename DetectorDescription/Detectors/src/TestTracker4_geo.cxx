#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/TGeoUnits.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;


static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)
{
    
    xml_det_t   x_det = e;

    DetElement  tracker(x_det.nameStr(), x_det.id());
    
    Volume experimentalHall_log =  lcdd.pickMotherVolume(tracker);
    
    xml_comp_t TrkElement (x_det.child(_U(box)));
    int counter = 0;
    
    int layerID = 0;
    
    for (xml_coll_t j(e,_U(layer)); j; ++j )
    {
    xml_comp_t x_layer = j;
        ++counter;
    
    double radius = x_layer.radius();
    
    double deltaphi = TrkElement.deltaphi();
    
    Box TrkBox(TrkElement.width(),TrkElement.thickness(),TrkElement.length());
    Volume TrkBoxVolume("TrkBoxVolume", TrkBox, lcdd.material(x_layer.materialStr()));
    
    //TrkBoxVolume.setSensitiveDetector(sens);
    
        double z = x_layer.dz();
        //double z = 10.*counter;
        
        for (int k = -6; k<7; k++)
        {
            
            for (int i = 0, n = TrkElement.repeat(); i < n; ++i)
            {
                double phi = deltaphi/dd4hep::rad * i;
                
                Position trans(radius * sin(phi),
                               radius * cos(phi),
                               k*z);
                
                PlacedVolume placedBox = experimentalHall_log.placeVolume(TrkBoxVolume, Transform3D(RotationZ(-phi-0.11*M_PI), trans));
                
                placedBox.addPhysVolID("system",x_det.id());
                placedBox.addPhysVolID("layer",layerID);
                placedBox.addPhysVolID("module", n);
                
                tracker.setPlacement(placedBox);
                
            }
            
        }
    
    tracker.setVisAttributes(lcdd, x_layer.visStr(), TrkBoxVolume);
    
  }
    return tracker;
}

DECLARE_DETELEMENT( TestTracker4, create_element )
