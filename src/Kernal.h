//
//  Kernal.h
//  Proj
//
//  Created by 余爽 on 2023/11/19.
//

#ifndef Kernal_h
#define Kernal_h

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace Eigen;


class Kernel {
public:
    static constexpr double RADIUS = 0.6; // neighbor radius
    static constexpr double RADIUSPOW3 = RADIUS * RADIUS * RADIUS;
    static constexpr double POLY6_COEF = 315 / (64 * M_PI * Kernel::RADIUSPOW3 * Kernel::RADIUSPOW3 * Kernel::RADIUSPOW3);
    static constexpr double SPIKY_GRADIENT_COEF = 45 / (M_PI * Kernel::RADIUSPOW3 * Kernel::RADIUSPOW3);
    
    static double poly6(double r) {
        double h = RADIUS;
        if (r <= h) {
            return std::pow(h * h - r * r, 3) * POLY6_COEF;
        } else {
            return 0;
        }
        
    }
    static double poly6dr(double r) {
        double h = RADIUS;
        if (r < h) {
            double sqrDiff = pow(h-r,2);
            return sqrDiff * SPIKY_GRADIENT_COEF;
        }
        else {
            return 0;
        }
    }
    static Vector3d poly6dp(Vector3d dp) {
        double h = RADIUS;
        double r = dp.norm();
        if (r < h) {
            double sqrDiff = pow(h - r, 2);
            double dWdr = sqrDiff * SPIKY_GRADIENT_COEF;
            dp.normalize();
            return -dp * dWdr;
        }
        else {
            return Vector3d(0, 0, 0);
        }
    }
};


#endif /* Kernal_h */
