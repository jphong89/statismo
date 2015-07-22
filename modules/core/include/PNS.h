#ifndef __PNS_H__
#define __PNS_H__

#include <memory>
#include <vector>
using std::vector;

#include "CommonTypes.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

#include <unsupported/Eigen/NonLinearOptimization>


// TODO: 
// 1. Note that data is stored as row major matrix
// 2. Make sure to cast as double
// 3. Use column major matrix for computation but as for interfacing with the outside, return as a row major matrix
// 4. Set up the class to use GTest framework
namespace statismo {
    template <typename T>
        class PNS{
            private:
                // enum for iteration type
                // TODO: SEQTEST to be implemented later
                typedef enum itype_t {
                    SEQTEST,
                    SMALL,
                    GREAT
                } itype;

                itype flag_;
                MatrixXd data_;
                // PNS data structures
                vector < VectorXd > orthaxis;
                VectorXd radii;
                VectorXd dist;
                //                 VectorXd pvalues; // TODO: To be used later

                MatrixXd computeRotMat( const VectorXd& vec );
                MatrixXd computeRiemannianExpMap( const MatrixXd& mat );
                MatrixXd computeRiemannianLogMap( const MatrixXd& mat );
                double modBy2PI( const double& x );
                double computeGeodesicMeanS1( const VectorXd& angles );
            public:
                // NOTE: punted the job to convert data into column major matrix to the caller
                PNS( const MatrixXd& data, const unsigned int flag = 2 ) : flag_( staic_cast<itype>( flag ) ), data_( data ) { } ;
                ~PNS();
                // TODO: Try to change magic numbers to enums I defined in fuctors
                double LMsphereFit( const MatrixXd& data, VectorXd& x, const int itype = 1);
                // internal objective function to be used inside getSubSphere
                double objFn( const VectorXd& center );
                double computeSubSphere( VectorXd& center);
                void compute(); // this corresponds to PNSmain.m
        };

};

} // namespace statismo

#include "PNS.hxx"

#endif /* __PNS_H_ */
