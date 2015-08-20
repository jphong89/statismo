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
#include "PNS.h"

namespace statismo {

    template <typename T>
        class PNSModelBuilder : public ModelBuilder<T> {


            public:

                typedef ModelBuilder<T> Superclass;
                typedef typename Superclass::DataManagerType DataManagerType;
                typedef typename Superclass::StatisticalModelType StatisticalModelType;
                typedef typename DataManagerType::DataItemListType DataItemListType;

                typedef enum { JacobiSVD, SelfAdjointEigenSolver } EigenValueMethod;

                /**
                 * Factory method to create a new PNSModelBuilder
                 */
                static PNSModelBuilder* Create() {
                    return new PNSModelBuilder();
                }

                /**
                 * Destroy the object.
                 * The same effect can be achieved by deleting the object in the usual
                 * way using the c++ delete keyword.
                 */
                void Delete() {
                    delete this;
                }


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
                StatisticalModelType* BuildNewModel(const DataItemListType& samples, double noiseVariance, bool computeScores = true, EigenValueMethod method = JacobiSVD, itype flag = static_cast<itype>(2) ) const;

            private:
                // to prevent use
                PNSModelBuilder();
                PNSModelBuilder(const PNSModelBuilder& orig);
                PNSModelBuilder& operator=(const PNSModelBuilder& rhs);
                StatisticalModelType* BuildNewModelInternal(const Representer<T>* representer, const Eigen::MatrixXd& X, double noiseVariance, EigenValueMethod method = JacobiSVD) const;
        };

} // namespace statismo

#include "PNSModelBuilder.hxx"

#endif /* __PNSMODELBUILDER_H_ */
