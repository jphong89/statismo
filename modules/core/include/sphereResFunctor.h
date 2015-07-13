#ifndef __SPHERERESFUNCTOR_H__
#define __SPHERERESFUNCTOR_H__
#include <Eigen/Dense>
#include <cmath>
#include "GenericFunctor.h"

class sphereResFunctor : public GenericFunctor<double> {
    public:
        typedef enum mode_t{
            ST=0,
            SC=1,
            GC=2
        } mode_t;
        mode_t m_mode;
        sphereResFunctor(int sizeX, int sizeY, const ValueType& y, int mode); 
        int operator() (const VectorXd& x, VectorXd& fvec) const;
        int df(const VectorXd& x, MatrixXd& fjac) const;
};


#include "sphereResFunctor.hxx"

#endif /*__SPHERERESFUNCTOR_H__*/
