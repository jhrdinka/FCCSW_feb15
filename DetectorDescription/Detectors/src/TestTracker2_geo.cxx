#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/TGeoUnits.h"
#include "DetExtensions/DetCylinderLayer.h"
#include "DetExtensions/Extension.h"

using namespace std;
using namespace DD4hep;
using namespace DD4hep::Geometry;

//Number of Elements in pos z
#define N_ELEMENTS 6


static Ref_t create_element(LCDD& lcdd, xml_h e, SensitiveDetector sens)
{
    
    xml_det_t   x_det = e;
    string det_name = x_det.nameStr();
    Material air = lcdd.air();
    
    //Detector mother volume
    DetElement tracker(det_name, x_det.id());
    //add Extension to Detlement for the RecoGeometry
    Det::Extension* ex = new Det::Extension();
    tracker.addExtension<Det::IExtension> (ex);
    //Place subdetector envelope into its mother/world volume
    Volume mother_vol =  lcdd.pickMotherVolume(tracker);
    
    //Create the Detector mother volume (detector envelope)
    DD4hep::XML::Dimension x_det_dim(x_det.dimensions());
    double z = x_det_dim.z();
    Tube tracker_shape(x_det_dim.rmin(),x_det_dim.rmax(),z);
    Volume tracker_vol(x_det.nameStr()+"_envelope",tracker_shape, air);
    
    tracker_vol.setVisAttributes(lcdd.invisible());
    
    //Create Box
    xml_comp_t TrkElement (x_det.child(_U(box)));
    Box TrkBox(TrkElement.width(),TrkElement.thickness(),TrkElement.length());
    sens.setType("tracker");
    
    int layer_num = 0;
    
    //Go through layers
    for (xml_coll_t j(e,_U(layer)); j; ++j )
    {
        xml_comp_t x_layer = j;
        double rmin = x_layer.inner_r();
        double rmax = x_layer.outer_r();
        double radius = (rmax+rmin)*0.5;
        double layer_z   = x_layer.z();
        
        //Create Volume and DetElement for Layer
        string layer_name  = det_name + _toString(layer_num,"layer%d");
        Volume layer_vol(layer_name,Tube(rmin,rmax,layer_z), lcdd.material(x_layer.materialStr()));
        DetElement lay_det (tracker,layer_name,layer_num);
        
        //Visualization
        layer_vol.setVisAttributes(lcdd.invisible());
        
        xml_comp_t x_module = x_layer.child(_U(module));
        xml_comp_t x_slice = x_layer.child(_U(slice));
        
        int repeat = x_module.repeat();
        double deltaphi = 2.*M_PI/repeat;
        int zrepeat = x_slice.repeat();
        double dz = x_slice.z();
        
        //add Extension to Detlement for the RecoGeometry
        Det::DetCylinderLayer* detcylinderlayer = new Det::DetCylinderLayer(repeat,2*zrepeat);
        lay_det.addExtension<Det::IExtension>(detcylinderlayer);
        
        //Create Box Volume
        Volume TrkBoxVolume("TrkBoxVolume", TrkBox, lcdd.material(x_module.materialStr()));
        //Set Box Volume sensitive
        if (x_module.isSensitive()) {
            TrkBoxVolume.setSensitiveDetector(sens);
        }
        
        int module_num = 0;
        //Place the Boxes in z and phi
        for (int k = -zrepeat; k<=zrepeat; k++)
        {
            string zname = _toString(k,"z%d");
            
            for (int i = 0; i < repeat; ++i)
            {
                double phi = deltaphi/dd4hep::rad * i;
                string module_name = zname + _toString(i,"module%d");
                
                Position trans(radius * sin(phi),
                               radius * cos(phi),
                               k*dz);
                //Place Box Volumes in layer -- question: Transform global or to layer??
                PlacedVolume placedBox = layer_vol.placeVolume(TrkBoxVolume, Transform3D(RotationZ(-phi-0.11*M_PI), trans));
                
                placedBox.addPhysVolID("module", module_num);
                //Create Box (Module) DetElement and assign it to the placed Box volume
                DetElement mod_det(lay_det,module_name,module_num);
                mod_det.setPlacement(placedBox);
                //add Extension to Detlement for the RecoGeometry
                mod_det.addExtension<Det::IExtension> (ex);
                
                ++module_num;
                
            }
            
            ++module_num;
            
        }
        //Place Layervolume

        PlacedVolume placedLayer = tracker_vol.placeVolume(layer_vol);
        placedLayer.addPhysVolID("layer",layer_num);
        placedLayer.addPhysVolID("system",x_det.id());
        //Assign Layer DetElement to LayerVolume
        lay_det.setPlacement(placedLayer);
        ++layer_num;
        
        tracker.setVisAttributes(lcdd, x_module.visStr(), TrkBoxVolume);
        
    }
    //Place envelopevolume in mothervolume
    PlacedVolume placed_env = mother_vol.placeVolume(tracker_vol);
    //assign tracker DetElement to tracker volume
    tracker.setPlacement(placed_env); //fuer envelope moeglich
    
    
    return tracker;
}

DECLARE_DETELEMENT( TestTracker2, create_element )
