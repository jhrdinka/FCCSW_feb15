//
//  Intersection.h
//  
//
//  Created by Julia Hrdinka on 26/11/14.
//
//

#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "Algebra/AlgPrimitives.h"

namespace Trk {

struct Intersection {
    
    Alg::Point3D position;
    double pathlength;
    bool onSurface;
    
    Intersection(Alg::Point3D pos, double path, bool onS) {
        position    = pos;
        pathlength  = path;
        onSurface   = onS;
    }

    };
    
}
#endif //INTERSECTION_H
