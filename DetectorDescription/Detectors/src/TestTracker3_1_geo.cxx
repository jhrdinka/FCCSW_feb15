//
//  TestTracker3_geo.cxx
//  
//
//  Created by Julia Hrdinka on 15/12/14.
//
//

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
         tracker_vol.setVisAttributes(lcdd.invisible());
        
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
        Volume mod_vol("module", Box(x_module.width(), x_module.thickness(),x_module.length()), air);
        mod_vol.setVisAttributes(lcdd.invisible());
        
        int comp_num = 0;
        //get components of the module
        for (xml_coll_t n(x_module,_U(module_component)); n; ++n) {
            xml_comp_t x_comp = n;
            
            Volume comp_vol("component " + x_comp.materialStr(), Box(x_comp.width(),x_comp.thickness(), x_comp.length()),lcdd.material(x_comp.materialStr()));
            comp_vol.setVisAttributes(lcdd, x_comp.visStr());
            //Set Sensitive Volmes sensitive
            if (x_comp.isSensitive()) {
                comp_vol.setSensitiveDetector(sens);
            }
            //Create DetElement
            DetElement comp_det (lay_det, "component " + x_comp.materialStr(),comp_num);
            //add Extension
            comp_det.addExtension<Det::Extension> (ex);
            //place component in Module
            Position trans (0.,x_comp.z(),0.);
            PlacedVolume placedcomp = mod_vol.placeVolume(comp_vol,trans);
            //assign the placed Volume to the DetElement
            comp_det.setPlacement(placedcomp);
            placedcomp.addPhysVolID("component",comp_num);
            ++comp_num;
        }
        int module_num = 0;
        //Place the Modules in z and phi
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
                PlacedVolume placedmodule = layer_vol.placeVolume(mod_vol, Transform3D(RotationZ(-phi-0.11*M_PI), trans));
                
                placedmodule.addPhysVolID("module", module_num);
                //Create Box (Module) DetElement and assign it to the placed Box volume
                DetElement mod_det(lay_det,module_name,module_num);
                mod_det.setPlacement(placedmodule);
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
        
    }
    //Place envelopevolume in mothervolume
    PlacedVolume placed_env = mother_vol.placeVolume(tracker_vol);
    //assign tracker DetElement to tracker volume
    tracker.setPlacement(placed_env); //fuer envelope moeglich
    
    
    return tracker;
}

DECLARE_DETELEMENT( TestTracker31, create_element )
