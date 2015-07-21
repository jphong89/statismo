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
#include <unsupported/Eigen/NonLinearOptimization>

namespace statismo {

template <typename T>
class PNSModelBuilder : public ModelBuilder<T> {


            public:

  public:


                /**
                 * The desctructor
                 */
                virtual ~PNSModelBuilder() {}

                /**
                 * Build a new model from the training data provided in the dataManager.
                 * \param samples A sampleSet holding the data
                 * \param noiseVariance The variance of N(0, noiseVariance) distributed noise on the points.
                 * If this parameter is set to 0, we have a standard PNS model. For values > 0 we have a PPNS model.
                 * \param computeScores Determines whether the scores (the pca coefficients of the examples) are computed and stored as model info
                 * (computing the scores may take a long time for large models).
                 * \param method Specifies the method which is used for the decomposition resp. eigenvalue solver.
                 *
                 * \return A new Statistical model
                 * \warning The method allocates a new Statistical Model object, that needs to be deleted by the user.
                 */
                StatisticalModelType* BuildNewModel(const DataItemListType& samples, double noiseVariance, bool computeScores = true, EigenValueMethod method = JacobiSVD) const;

            private:
                // to prevent use
                PNSModelBuilder();
                PNSModelBuilder(const PNSModelBuilder& orig);
                PNSModelBuilder& operator=(const PNSModelBuilder& rhs);
                StatisticalModelType* BuildNewModelInternal(const Representer<T>* representer, const MatrixType& X, double noiseVariance, EigenValueMethod method = JacobiSVD) const;


                // PNS stuff
                // Things to consider: Do I want to return pointer instead of a value??
                // If I were to return as a pointer, then managing the memory would be a pain
                // Probably good to take a look at smart pointer stuff and understand how that works


                typedef NestedSphere_t {
                    double radius;
                } NestedSphere;

                MatrixXd computeRotMat( const VectorXd& vec ) const;
                MatrixXd computeRiemannianExpMap( const MatrixXd& mat ) const;
                MatrixXd computeRiemannianLogMap( const MatrixXd& mat ) const;
                double computeGeodesicMeanS1( const VectorXd& angles ) const;
                double modBy2PI( const double& x ) const;
                // We may want to have LM optimizer as a member
                // TODO: Try to change magic numbers to enums I defined in fuctors
                double LMsphereFit( const MatrixXd& data, VectorXd& x, const int itype = 1) const;
                // internal objective function to be used inside getSubSphere
                double objFn( const VectorXd& center, const MatrixXd& data, const double r ) const;
                double computeSubSphere( VectorXd& center, const MatrixXd& data, const int itype = 1, double& error ) const;
                double getSubSphere( VectorXd& center, const MatrixXd& data, const int itype = 1) const;
                MatrixXd PNSmain( const MatrixXd& data, const int itype = 1, MatrixXd& residualMat ) const;
        };

} // namespace statismo

#include "PNSModelBuilder.hxx"

#endif /* __PNSMODELBUILDER_H_ */
