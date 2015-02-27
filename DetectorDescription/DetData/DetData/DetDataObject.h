//
//  DetDataObject.h
//  
//
//  Created by Julia Hrdinka on 22/10/14.
//
//

#ifndef DET_DETDATAOBJECT_H
#define DET_DETDATAOBJECT_H

#include "GaudiKernel/DataObject.h"
#inlcude "RecoGeometry/ContainerVolume.h"

namespace Det {
    
    class GAUDI_API DetDataObject: public DataObject {
        
    public:
        
        DetDataObject() : m_worldVolume(0)
        {}
        
        virtual ~DetDataObject()
        {
            if (m_worldVolume!=0) {
                delete m_worldVolume;
            }
        }
        
        Reco::ContainerVolume* getWorldVolume()
        {
            return m_worldVolume;
        }
        
        void setWorldVolume (Reco::ContainerVolume worldVolume)
        {
            m_worldVolume = worldVolume;
        }
        
    private:
        
        Reco::ContainerVolume* m_worldVolume;
    };
}

#endif //DET_DETDATAOBJECT_H
