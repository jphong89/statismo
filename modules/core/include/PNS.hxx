#ifndef __PNS_TXX__
#define __PNS_TXX__

#include "PNS.h"

namespace statismo {
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
