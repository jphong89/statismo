/*
 * This file is part of the statismo library.
 *
 * Author: Marcel Luethi (marcel.luethi@unibas.ch)
 *
 * Copyright (c) 2011 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef __PNSModelBuilder_TXX
#define __PNSModelBuilder_TXX

#include "PNSModelBuilder.h"

#include <iostream>

#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

#include "CommonTypes.h"
#include "Exceptions.h"

namespace statismo {

    template <typename T>
        PNSModelBuilder<T>::PNSModelBuilder()
        : Superclass() {
        }


    template <typename T>
        typename PNSModelBuilder<T>::StatisticalModelType*
        PNSModelBuilder<T>::BuildNewModel(const DataItemListType& sampleDataList, double noiseVariance, bool computeScores, EigenValueMethod method) const {

            unsigned n = sampleDataList.size();
            if (n <= 0) {
                throw StatisticalModelException("Provided empty sample set. Cannot build the sample matrix");
            }

            unsigned p = sampleDataList.front()->GetSampleVector().rows();
            const Representer<T>* representer = sampleDataList.front()->GetRepresenter();

            // Build the sample matrix X
            MatrixType X(n, p);

            unsigned i = 0;
            for (typename DataItemListType::const_iterator it = sampleDataList.begin();
                    it != sampleDataList.end();
                    ++it) {
                assert ((*it)->GetSampleVector().rows() == p); // all samples must have same number of rows
                assert ((*it)->GetRepresenter() == representer); // all samples have the same representer
                X.row(i++) = (*it)->GetSampleVector();
            }


            // build the model
            StatisticalModelType* model = BuildNewModelInternal(representer, X, noiseVariance, method);
            MatrixType scores;
            if (computeScores) {
                scores = this->ComputeScores(X, model);
            }


            typename BuilderInfo::ParameterInfoList bi;
            bi.push_back(BuilderInfo::KeyValuePair("NoiseVariance ", Utils::toString(noiseVariance)));

            typename BuilderInfo::DataInfoList dataInfo;
            i = 0;
            for (typename DataItemListType::const_iterator it = sampleDataList.begin();
                    it != sampleDataList.end();
                    ++it, i++) {
                std::ostringstream os;
                os << "URI_" << i;
                dataInfo.push_back(BuilderInfo::KeyValuePair(os.str().c_str(),(*it)->GetDatasetURI()));
            }


            // finally add meta data to the model info
            BuilderInfo builderInfo("PNSModelBuilder", dataInfo, bi);

            ModelInfo::BuilderInfoList biList;
            biList.push_back(builderInfo);

            ModelInfo info(scores, biList);
            model->SetModelInfo(info);

            return model;
        }


    template <typename T>
        typename PNSModelBuilder<T>::StatisticalModelType*
        PNSModelBuilder<T>::BuildNewModelInternal(const Representer<T>* representer, const MatrixType& X, double noiseVariance, EigenValueMethod method) const {

            unsigned n = X.rows();
            unsigned p = X.cols();

            RowVectorType mu = X.colwise().mean(); // needs to be row vector
            MatrixType X0 = X.rowwise() - mu;

            switch(method) {
                case JacobiSVD:

                    typedef Eigen::JacobiSVD<MatrixType> SVDType;
                    typedef Eigen::JacobiSVD<MatrixTypeDoublePrecision> SVDDoublePrecisionType;

                    // We destinguish the case where we have more variables than samples and
                    // the case where we have more samples than variable.
                    // In the first case we compute the (smaller) inner product matrix instead of the full covariance matrix.
                    // It is known that this has the same non-zero singular values as the covariance matrix.
                    // Furthermore, it is possible to compute the corresponding eigenvectors of the covariance matrix from the
                    // decomposition.

                    if (n < p) {
                        // we compute the eigenvectors of the covariance matrix by computing an SVD of the
                        // n x n inner product matrix 1/(n-1) X0X0^T
                        MatrixType Cov = X0 * X0.transpose() * 1.0/(n-1);
                        SVDDoublePrecisionType SVD(Cov.cast<double>(), Eigen::ComputeThinV);
                        VectorType singularValues = SVD.singularValues().cast<ScalarType>();
                        MatrixType V = SVD.matrixV().cast<ScalarType>();

                        unsigned numComponentsAboveTolerance = ((singularValues.array() - noiseVariance - Superclass::TOLERANCE) > 0).count();

                        // there can be at most n-1 nonzero singular values in this case. Everything else must be due to numerical inaccuracies
                        unsigned numComponentsToKeep = std::min(numComponentsAboveTolerance, n - 1);
                        // compute the pseudo inverse of the square root of the singular values
                        // which is then needed to recompute the PNS basis
                        VectorType singSqrt = singularValues.array().sqrt();
                        VectorType singSqrtInv = VectorType::Zero(singSqrt.rows());
                        for (unsigned i = 0; i < numComponentsToKeep; i++) {
                            assert(singSqrt(i) > Superclass::TOLERANCE);
                            singSqrtInv(i) = 1.0 / singSqrt(i);
                        }

                        // we recover the eigenvectors U of the full covariance matrix from the eigenvectors V of the inner product matrix.
                        // We use the fact that if we decompose X as X=UDV^T, then we get X^TX = UD^2U^T and XX^T = VD^2V^T (exploiting the orthogonormality
                        // of the matrix U and V from the SVD). The additional factor sqrt(n-1) is to compensate for the 1/sqrt(n-1) in the formula
                        // for the covariance matrix.
                        MatrixType pcaBasis = (X0.transpose() * V * singSqrtInv.asDiagonal() / sqrt(n-1.0)).topLeftCorner(p, numComponentsToKeep);;

                        if (numComponentsToKeep == 0) {
                            throw StatisticalModelException("All the eigenvalues are below the given tolerance. Model cannot be built.");
                        }

                        VectorType sampleVarianceVector = singularValues.topRows(numComponentsToKeep);
                        VectorType pcaVariance = (sampleVarianceVector - VectorType::Ones(numComponentsToKeep) * noiseVariance);

                        StatisticalModelType* model = StatisticalModelType::Create(representer, mu, pcaBasis, pcaVariance, noiseVariance);

                        return model;
                    } else {
                        // we compute an SVD of the full p x p  covariance matrix 1/(n-1) X0^TX0 directly
                        SVDType SVD(X0.transpose() * X0 * 1.0/(n-1), Eigen::ComputeThinU);
                        VectorType singularValues = SVD.singularValues();
                        unsigned numComponentsToKeep = ((singularValues.array() - noiseVariance - Superclass::TOLERANCE) > 0).count();
                        MatrixType pcaBasis = SVD.matrixU().topLeftCorner(p, numComponentsToKeep);

                        if (numComponentsToKeep == 0) {
                            throw StatisticalModelException("All the eigenvalues are below the given tolerance. Model cannot be built.");
                        }

                        VectorType sampleVarianceVector = singularValues.topRows(numComponentsToKeep);
                        VectorType pcaVariance = (sampleVarianceVector - VectorType::Ones(numComponentsToKeep) * noiseVariance);
                        StatisticalModelType* model = StatisticalModelType::Create(representer, mu, pcaBasis, pcaVariance, noiseVariance);
                        return model;
                    }
                    break;

                case SelfAdjointEigenSolver: {
                                                 // we compute the eigenvalues/eigenvectors of the full p x p  covariance matrix 1/(n-1) X0^TX0 directly

                                                 typedef Eigen::SelfAdjointEigenSolver<MatrixType> SelfAdjointEigenSolver;
                                                 SelfAdjointEigenSolver es;
                                                 es.compute(X0.transpose() * X0 * 1.0/(n-1));
                                                 VectorType eigenValues = es.eigenvalues().reverse(); // SelfAdjointEigenSolver orders the eigenvalues in increasing order


                                                 unsigned numComponentsToKeep = ((eigenValues.array() - noiseVariance - Superclass::TOLERANCE) > 0).count();
                                                 MatrixType pcaBasis = es.eigenvectors().rowwise().reverse().topLeftCorner(p, numComponentsToKeep);

                                                 if (numComponentsToKeep == 0) {
                                                     throw StatisticalModelException("All the eigenvalues are below the given tolerance. Model cannot be built.");
                                                 }

                                                 VectorType sampleVarianceVector = eigenValues.topRows(numComponentsToKeep);
                                                 VectorType pcaVariance = (sampleVarianceVector - VectorType::Ones(numComponentsToKeep) * noiseVariance);
                                                 StatisticalModelType* model = StatisticalModelType::Create(representer, mu, pcaBasis, pcaVariance, noiseVariance);
                                                 return model;
                                             }
                                             break;

                default:
                                             throw StatisticalModelException("Unrecognized decomposition/eigenvalue solver method.");
                                             return 0;
                                             break;
            }
        }

    template <typename T>
        typename PNSModelBuilder<T>::MatrixXd
        PNSModelBuilder<T>::computeRotMat( const VectorXd& vec ) const {

            Eigen::MatrixXd rotMat( vec.size(), vec.size() );
            roMat.setIdentity();

            Eigen::VectorXd northPole = Eigen::Zero( vec.size() );
            northPole( vec.size() - 1 ) = 1;

            Eigen::VectorXd vecUnit = vec;
            vecUnit.normalize();

            double mult = northPole.adjoint() * vecUnit;
            double alpha = acos( mult );

            if ( abs( mult -1) < 1e-15 ) {
                return rotMat;
            }
            if ( abs( mult +1) < 1e-15 ) {
                rotMat*= -1;
                return rotMat;
            }
            Eigen::VectorXd auxVec = vecUnit - mult*northPole;
            aux.normalize();
            Eigen::MatrixXd auxMat = vecUnit*aux.transpose() - aux*vecUnit.transpose();
            rotMat += sin(alpha)*auxMat + (cos(alpha) - 1)*(vecUnit*vecUnit.transpose() + auxVec*auxVec.transpose());
            return rotMat;
        }
    template <typename T>
        typename PNSModelBuilder<T>::MatrixXd 
        PNSModelBuilder<T>::computeRiemannianExpMap( const MatrixXd& mat ) const {
            
            MatrixXd result( mat.rows()+1, mat.cols() );
            
            RowVectorXd normVec     = mat.colwise().norm();
            RowVectorXd normVec2    = (normVec.array() < 1e-16).select(0,normVec);
            RowVectorXd sineVec     = normVec2.array().sin();
            RowVectorXd cosVec      = normVec2.array().cos();

            RowVectorXd auxVec      = sineVec.array() / normVec.array();            
            MatrixXd    auxMat      = auxVec.replicate( mat.rows(), 1 );
            auxMat                  = (auxMat.array() * mat.array()).matrix();
            result << auxMat, cosVec;
            
            return result;

        }
    template <typename T>
        typename PNSModelBuilder<T>::MatrixXd
        PNSModelBuilder<T>::computeRiemannianLogMap( const MatrixXd& mat ) const {
            
            MatrixXd result(mat.rows()-1, mat.cols());
            
            RowVectorXd lastRow = mat.row(mat.rows() - 1);
            RowVectorXd auxV1   = (lastRow.array() * lastRow.array()).matrix();
                        auxV1   = (1 - auxV1.array()).matrix();
            RowVectorXd auxV2   = (lastRow.array().acos()).matrix();
            RowVectorXd auxV3   = auxV1.array() / auxV2.array();
            RowVectorXd scale   = (auxV1.array() < 1e-64 && auxV2.array() < 1e-64).select(1,auxV3);
            MatrixXd    auxM1   = scale.replicate( mat.rows() - 1, 1 );
            MatrixXd    auxM2   = (mat.topRows( mat.rows() - 1 ));
                        result  = auxM1.array() * auxM2.array();

        }
} // namespace statismo
#endif