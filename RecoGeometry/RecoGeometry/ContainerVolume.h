//
//  ContainerVolume.h
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#ifndef RECO_CONTAINERVOLUME_H
#define RECO_CONTAINERVOLUME_H

#include "Algebra/AlgPrimitives.h"
#include "TrkGeometryUtils/BinnedArray.h"
#include "RecoGeometry/Volume.h"

namespace Reco {
    
    class ContainerVolume : public virtual Volume {
        
    public:
        
        typedef Trk::BinnedArray<Volume> VolumeArray;
        
        ContainerVolume(VolumeArray* volumes);
        virtual bool isInside(const Alg::Point3D& glopos, double tol=0.) const = 0;
        virtual ~ContainerVolume();
        
    private:
        VolumeArray* m_volumes;
    };
}


#endif //RECO_CONTAINERVOLUME_H
