//
//  Surface.h
//  
//
//  Created by Julia Hrdinka on 25/11/14.
//
//

#ifndef RECO_SURFACE_H
#define RECO_SURFACE_H

#include "Algebra/AlgPrimitives.h"
#include "TGeoNode.h"

#include "TrkParametersBase/ParametersBase.h"
#include "TrkParametersBase/Charged.h"
#include "TrkParametersBase/Neutral.h"

#include "TrkGeometryUtils/Intersection.h"

#include "RecoGeometry/Material.h"
#include "RecoGeometry/MaterialMap.h"

namespace Reco {
    
    class Surface {
        
    public:
        //constructor only with TGeoNode, creates a Surface without material, is used only for layers
        Surface (TGeoNode* node);
        //constructor with TGeoNode and transformation matrix
 //       Surface (TGeoNode* node, Alg::Transform3D* transf);
        //standard constructor with TGeoNode and MateriaMap
        Surface (TGeoNode* node, MaterialMap* materialmap);
        Surface (MaterialMap* materialmap, std::shared_ptr<const Alg::Transform3D> transf);
        Surface (std::shared_ptr<const Alg::Transform3D> transf);
        
        //destructor
        virtual ~Surface();
        //assignment operator
        Surface& operator=(const Surface& surface);
        //get transform matrix
        const Alg::Transform3D& transform() const;
        //get the center of the surface
        const Alg::Point3D& center() const;
        //get the normal vector of the surface
        virtual const Alg::Vector3D& normal() const;
        virtual const Alg::Vector3D* normal(const Alg::Point2D& locpos) const;
        //get the MaterialMap of the Surface
        const MaterialMap& materialmap() const;
        //get material on at specific local point of the Surface
        const Material& material(Alg::Point2D& locpos) const;
        
        virtual bool isInside(const Alg::Point2D& locpos, double tol1, double tol2) const = 0;
        
        virtual void localToGlobal(const Alg::Point2D& locpos, const Alg::Vector3D& mom, Alg::Point3D& glopos) const = 0;
        virtual bool globalToLocal(const Alg::Point3D& glopos, const Alg::Vector3D& mom, Alg::Point2D& locpos) const = 0;
        
        /** Use the Surface as a ParametersBase constructor, from local parameters - charged */
        virtual const Trk::ParametersBase<5, Trk::Charged>* createTrackParameters(double, double, double, double, double, Alg::AmgSymMatrix<5>* cov = 0) const = 0;
        
        /** Use the Surface as a ParametersBase constructor, from global parameters - charged*/
        virtual const Trk::ParametersBase<5, Trk::Charged>* createTrackParameters(const Alg::Point3D&, const Alg::Vector3D&, double, Alg::AmgSymMatrix<5>* cov = 0) const = 0;
        
        /** Use the Surface as a ParametersBase constructor, from local parameters - neutral */
        virtual const Trk::ParametersBase<5, Trk::Neutral>* createNeutralParameters(double, double, double, double, double, Alg::AmgSymMatrix<5>* cov = 0) const = 0;
        
        /** Use the Surface as a ParametersBase constructor, from global parameters - neutral */
        virtual const Trk::ParametersBase<5, Trk::Neutral>* createNeutralParameters(const Alg::Point3D&, const Alg::Vector3D&, double charge=0., Alg::AmgSymMatrix<5>* cov = 0) const = 0;
        
        /** fast straight line intersection schema - standard: provides closest intersection and (signed) path length
         forceFwd is to provide the closest forward solution
         */
        virtual Trk::Intersection straightLineIntersection(const Alg::Point3D& pos,
         const Alg::Vector3D& dir,
         bool forceDir = false) const = 0;
        
                
    protected:
        
        Surface ();
        
        mutable Alg::Point3D*               m_center;
        mutable Alg::Vector3D*              m_normal;
        std::shared_ptr<const Alg::Transform3D>  m_transform;
        
        TGeoNode*                   m_node;
        MaterialMap*                m_materialmap;
        
       //constants
        static double               s_onSurfaceTolerance;
        static Alg::Point3D         s_origin;
        static Alg::Transform3D     s_idTransform;
        
    };
}

#endif //RECO_SURFACE_H)