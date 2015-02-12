//
//  PlaneSurface.h
//  
//
//  Created by Julia Hrdinka on 25/11/14.
//
//

#ifndef RECO_PLANESURFACE_H
#define RECO_PLANESURFACE_H

#include "RecoGeometry/Surface.h"
#include "TGeoBBox.h"

#include "TrkParametersBase/ParametersT.h"

namespace Reco {
   
    class PlaneSurface : public Surface {
    
    public:
        
        PlaneSurface(TGeoNode* node, TGeoBBox* box);
        PlaneSurface(TGeoBBox* box, std::shared_ptr<const Alg::Transform3D> transf);
        PlaneSurface(TGeoNode* node, TGeoBBox* box, MaterialMap* materialmap);
        PlaneSurface(TGeoBBox* box, MaterialMap* materialmap, std::shared_ptr<const Alg::Transform3D> transf);
        PlaneSurface(std::shared_ptr<const Alg::Transform3D> transf, double halfX, double halfY);
        virtual ~PlaneSurface();
        
        PlaneSurface& operator=(const PlaneSurface& planesurface);
        
        double getHalfX() const;
        double getHalfY() const;
        
        //get the normal vector of the surface
        virtual const Alg::Vector3D& normal() const override;
        virtual const Alg::Vector3D* normal(const Alg::Point2D& locpos) const override;
        
        
        virtual bool isInside(const Alg::Point2D& locpos, double tol1, double tol2) const override;
        
        virtual void localToGlobal(const Alg::Point2D& locpos, const Alg::Vector3D& mom, Alg::Point3D& glopos) const override;
        virtual bool globalToLocal(const Alg::Point3D& glopos, const Alg::Vector3D& mom, Alg::Point2D& locpos) const override;
        
        /** Use the Surface as a ParametersBase constructor, from local parameters - charged */
        virtual const Trk::ParametersT<5, Trk::Charged, PlaneSurface>* createTrackParameters(double l1,
                                                                                             double l2,
                                                                                             double phi,
                                                                                             double theta,
                                                                                             double qop,
                                                                                             Alg::AmgSymMatrix<5>* cov = 0) const override
        { return new Trk::ParametersT<5, Trk::Charged, PlaneSurface>(l1, l2, phi, theta, qop, *this, cov); }
        
        
        /** Use the Surface as a ParametersBase constructor, from global parameters - charged*/
        virtual const Trk::ParametersT<5, Trk::Charged, PlaneSurface>* createTrackParameters(const Alg::Point3D& position,
                                                                                             const Alg::Vector3D& momentum,
                                                                                             double charge,
                                                                                             Alg::AmgSymMatrix<5>* cov = 0) const override
        { return new Trk::ParametersT<5, Trk::Charged, PlaneSurface>(position, momentum, charge, *this, cov); }
        
        /** Use the Surface as a ParametersBase constructor, from local parameters - neutral */
        virtual const Trk::ParametersT<5, Trk::Neutral, PlaneSurface>* createNeutralParameters(double l1,
                                                                                               double l2,
                                                                                               double phi,
                                                                                               double theta,
                                                                                               double oop,
                                                                                               Alg::AmgSymMatrix<5>* cov = 0) const override
        { return new Trk::ParametersT<5, Trk::Neutral, PlaneSurface>(l1, l2, phi, theta, oop, *this, cov); }
        
        /** Use the Surface as a ParametersBase constructor, from global parameters - neutral */
        virtual const Trk::ParametersT<5, Trk::Neutral, PlaneSurface>* createNeutralParameters(const Alg::Point3D& position,
                                                                                               const Alg::Vector3D& momentum,
                                                                                               double charge = 0.,
                                                                                               Alg::AmgSymMatrix<5>* cov = 0) const override
        { return new Trk::ParametersT<5, Trk::Neutral, PlaneSurface>(position, momentum, charge, *this, cov); }
        
        /** Use the Surface as a ParametersBase constructor, from local parameters */
        template <int DIM, class T> const Trk::ParametersT<DIM, T, PlaneSurface>* createParameters(double l1,
                                                                                                   double l2,
                                                                                                   double phi,
                                                                                                   double theta,
                                                                                                   double qop,
                                                                                                   Alg::AmgSymMatrix<DIM>* cov = 0) const
        { return new Trk::ParametersT<DIM, T, PlaneSurface>(l1, l2, phi, theta, qop, *this, cov); }
        
        /** Use the Surface as a ParametersBase constructor, from global parameters */
        template <int DIM, class T> const Trk::ParametersT<DIM, T, PlaneSurface>* createParameters(const Alg::Point3D& position,
                                                                                                   const Alg::Vector3D& momentum,
                                                                                                   double charge,
                                                                                                   Alg::AmgSymMatrix<DIM>* cov = 0) const
        { return new Trk::ParametersT<DIM, T, PlaneSurface>(position, momentum, charge, *this, cov); }
        
        /** fast straight line intersection schema - standard: provides closest intersection and (signed) path length
         forceDir is to provide the closest forward solution
         
         <b>mathematical motivation:</b>
         
         the equation of the plane is given by: <br>
         @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
         where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of the plane,
         @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on the plane and @f$ \vec x = (x,y,z) @f$ all possible points
         on the plane.<br>
         Given a line with:<br>
         @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
         the solution for @f$ u @f$ can be written:
         @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
         If the denominator is 0 then the line lies:
         - either in the plane
         - perpenticular to the normal of the plane
         
         */
        virtual Trk::Intersection straightLineIntersection(const Alg::Point3D& pos,
         const Alg::Vector3D& dir,
         bool forceDir) const override;
        
    protected:
        
        PlaneSurface();
        
        //halflength:
        double m_halfX;
        double m_halfY;
    };
    
    inline Trk::Intersection PlaneSurface::straightLineIntersection(const Alg::Point3D& pos,
     const Alg::Vector3D& dir,
     bool forceDir) const
     {
     double denom = dir.Dot(normal());
     if (denom){
     double u = (normal().Dot((center()-pos)))/(denom);
     Alg::Point3D intersectPoint(pos + u * dir);
     // evaluate the intersection in terms of direction
     bool isValid = forceDir ?  ( u > 0.) : true;
     //check, if on surface
     Alg::Point2D loc;
     if (isValid) isValid = (globalToLocal(pos,dir,loc) && isInside(loc,s_onSurfaceTolerance,s_onSurfaceTolerance));
     // return the result
         return Trk::Intersection(intersectPoint,u,isValid);
     }
         return Trk::Intersection(pos,0.,false);
     }
}


#endif //RECO_PLANESURFACE_H
