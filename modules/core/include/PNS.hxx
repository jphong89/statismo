#ifndef __PNS_TXX__
#define __PNS_TXX__

#include "PNS.h"

namespace statismo {
    template <typename T> 
        MatrixXd PNS<T>::computeRotMat( const VectorXd& vec ) {

            MatrixXd rotMat( vec.size(), vec.size() );
            roMat.setIdentity();

            VectorXd northPole = Eigen::Zero( vec.size() );
            northPole( vec.size() - 1 ) = 1;

            VectorXd vecUnit = vec;
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
            VectorXd auxVec = vecUnit - mult*northPole;
            aux.normalize();
            MatrixXd auxMat = northPole*aux.transpose() - aux*northPole.transpose();
            rotMat += sin(alpha)*auxMat + (cos(alpha) - 1)*(northPole*northPole.transpose() + auxVec*auxVec.transpose());

            return rotMat;
        }

    template <typename T>
        MatrixXd PNS<T>::computeRiemannianExpMap( const MatrixXd& mat ) {

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
        MatrixXd PNS<T>::computeRiemannianLogMap( const MatrixXd& mat ) {

            MatrixXd result(mat.rows()-1, mat.cols());

            RowVectorXd lastRow = mat.row(mat.rows() - 1);
            RowVectorXd auxV1   = (lastRow.array() * lastRow.array()).matrix();
            auxV1   = (1 - auxV1.array()).matrix();
            RowVectorXd auxV2   = (lastRow.array().acos()).matrix();
            RowVectorXd auxV3   = auxV1.array() / auxV2.array();
            // This line is to check for NaN occurrences.
            // If any entry in auxV3 is NaN then replace it with 1.
            RowVectorXd scale   = (auxV1.array() < 1e-64 && auxV2.array() < 1e-64).select(1,auxV3); 
            MatrixXd    auxM1   = scale.replicate( mat.rows() - 1, 1 );
            MatrixXd    auxM2   = (mat.topRows( mat.rows() - 1 ));
            result  = auxM1.array() * auxM2.array();
        }

    template <typename T>
        double PNS<T>::modBy2PI( const double& x ) {
            // helper function to be used.
            // Maybe consider inlining?
            return ( x - (2*PI)*floor( x / (2*PI) ) );
        }

    template <typename T>
        double PNS<T>::computeGeodesicMeanS1( const VectorXd& angles ) const {
            VectorXd meanCandidate( angles.size() );
            VectorXd auxV1( angles.size() );
            VectorXd theta( angles.size() );
            MatrixXd distMatrix( angles.size(), 3 );
            VectorXd geodVariance( angles.size() );
            double currCandidate(0);
            double currGeodVariance(0);
            int idxToGeodMean(0);
            double geodMean(0);


            // same as theta = mod( angles, 2*pi ) in MATLAB
            theta = angles.unaryExpr( std::ptr_fun(modBy2PI) );
            // Generating mean candidates
            auxV1.setLinSpaced(angles.size(), 0, angles.size()-1);
            auxV1 = (auxV1.array() / auxV1.size()).matrix();

            meanCandidate = (angles.mean() + 2*PI*auxV1.array()).matrix();
            meanCandidate.unaryExpr( std::ptr_fun(modBy2PI));

            for(int i=0; i<angles.size(); ++i) {
                double currCandidate = meanCandidate(i);
                distMatrix.col(0) = ((theta.array() - currCandidate).square()).matrix();
                distMatrix.col(1) = ((theta.array() - currCandidate + 2*PI).square()).matrix();
                distMatrix.col(2) = ((currCandidate - theta.array() + 2*PI).square()).matrix();
                currGeodVariance = (distMatrix.rowwise().minCoeff()).sum();
                geodVariance(i) = currGeodVariance;
            }
            double minGeodVarianceValue = geodVariance.minCoeff(&idxToGeodMean);
            geodMean = modBy2PI( meanCandidate( idxToGeodMean ) );
        }
} // namespace statismo
#endif // __PNS_TXX__
