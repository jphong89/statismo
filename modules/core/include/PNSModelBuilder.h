#ifndef __PNSMODELBUILDER_H_
#define __PNSMODELBUILDER_H_

#include <memory>
#include <vector>

#include "CommonTypes.h"
#include "Config.h"
#include "DataManager.h"
#include "ModelBuilder.h"
#include "ModelInfo.h"
#include "StatisticalModel.h"
// To use nonlinear least square solver
// TODO: Check the name as well as its API
#include <Eigen/unsupported>

namespace statismo {

template <typename T>
class PNSModelBuilder : public ModelBuilder<T> {


            public:

  public:

    typedef ModelBuilder<T> Superclass;
    typedef typename Superclass::DataManagerType DataManagerType;
    typedef typename Superclass::StatisticalModelType StatisticalModelType;
    typedef typename DataManagerType::DataItemListType DataItemListType;
    typedef typename Eigen::MatrixXd MatrixXd;
    typedef typename Eigen::VectorXd VectorXd;
    typedef typename Eigen::RowVectorXd RowVectorXd;

    typedef enum { JacobiSVD, SelfAdjointEigenSolver } EigenValueMethod;
    StatisticalModelType* BuildNewModelInternal(const Representer<T>* representer, const MatrixType& X, double noiseVariance, EigenValueMethod method = JacobiSVD) const;
    // PNS stuff
    // Do I want to return pointer instead of a value??
    MatrixXd computeRotMat( VectorXd vec ) const;
    MatrixXd computeRiemannianExpMap( MatrixXd mat ) const;


};

} // namespace statismo

#include "PNSModelBuilder.hxx"

#endif /* __PNSMODELBUILDER_H_ */
