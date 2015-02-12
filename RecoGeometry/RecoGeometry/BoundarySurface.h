//
//  BoundarySurface.h
//  
//
//  Created by Julia Hrdinka on 26/01/15.
//
//

#ifndef RECO_BOUNDARYSURFACE_H
#define RECO_BOUNDARYSURFACE_H

#include "TrkGeometryUtils/BinnedArray1D.h"

namespace Reco {

    class Volume;
    class Surface;
    
    class BoundarySurface  {
        
    public:
    
        typedef Trk::BinnedArray1D<Volume> VolumeArray;
        
  //      BoundarySurface(Volume* nextVolume, Volume* previousVolume);
 //       BoundarySurface(VolumeArray* nextVolumes, VolumeArray* previousVolumes);
        virtual ~BoundarySurface();
        BoundarySurface& operator=(const BoundarySurface& boundarysurface);
        
        virtual void setNextVolume(std::shared_ptr<const Volume> nextVolume) const;
        virtual void setPreviousVolume(std::shared_ptr<const Volume> previousVolume) const;
        virtual void setNextVolumes(const VolumeArray* nextVolumes) const;
        virtual void setPreviousVolumes(const VolumeArray* previousVolumes) const;
        virtual const Volume* getNextVolume() const;
        virtual const Volume* getPreviousVolume() const;
        virtual const VolumeArray* getNextVolumes() const;
        virtual const VolumeArray* getPreviousVolumes() const;
        
  //      virtual const Volume* getNextVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const = 0;
  //      virtual const Volume* getPreviousVolume(const Alg::Point3D& glopos, const Alg::Vector3D& dir) const = 0;
        virtual const Surface* surfaceRepresentation() const = 0;
//        BoundarySurfaceType type() const = 0;
        
    protected:
        
        BoundarySurface();
        //surrounding Volumes
//        VolumeArray* m_volumes;
        //binnedarray von next und previous volumes
        mutable std::shared_ptr<const Volume>   m_nextVolume;
        mutable std::shared_ptr<const Volume>   m_previousVolume;
        mutable const VolumeArray*              m_nextVolumes;
        mutable const VolumeArray*              m_previousVolumes;
    };
}
#endif //RECO_BOUNDARYSURFACE_H
