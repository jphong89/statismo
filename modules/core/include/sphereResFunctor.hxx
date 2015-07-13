#include "sphereResFunctor.h"

sphereResFunctor::sphereResFunctor(int sizeX, int sizeY, const ValueType& y, int mode) :GenericFunctor<double>(sizeX, sizeY, y){
    m_mode = (mode_t)mode;
}

int sphereResFunctor::operator() (const VectorXd& x, VectorXd& fvec) const {
    assert(x.size()==inputs());
    assert(fvec.size()==values());

    MatrixXd inputMatrix(inputs(), values());
    MatrixXd auxMatrix(inputs(), values());
    double radius(0);

    inputMatrix = x.replicate(1, values());
    auxMatrix   = ( (m_y.array() - inputMatrix.array()).square() ).matrix();

    VectorXd di = (auxMatrix.colwise().sum()).array().sqrt().matrix();
    if(m_mode != GC) {
        radius = di.array().mean();
    }
    else {
        radius = M_PI / 2;
    }

    fvec = ( di.array() - radius ).matrix();
    return 0;
}


int sphereResFunctor::df(const VectorXd &x, MatrixXd &fjac) const
{
    assert(fjac.rows() == values() && fjac.cols() == inputs());
    VectorXd rC(values());
    (*this)(x,rC);
    // TODO: I don't like these magic numbers!
    double TOL = 1e-10;
    double dx = 0.25*TOL;
    VectorXd rd(values());
    VectorXd xd(inputs());

    for(int i=0;i<inputs();++i) { 
        xd = x;
        xd(i) += dx;
        rd.setZero();
        (*this)(xd, rd);
        fjac.col(i) = (rd.array() - rC.array())/dx;
    }
    return 0;
}
