#ifndef __PNS_TXX__
#define __PNS_TXX__

#include "PNS.h"
#include <iostream>
using std::cout;
using std::endl;

namespace pns {
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




#endif // __PNS_TXX__
