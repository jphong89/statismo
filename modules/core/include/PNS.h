#ifndef __PNS_H__
#define __PNS_H__

#include <memory>
#include <vector>
using std::vector;

#include "CommonTypes.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::JacobiSVD;
using Eigen::ComputeThinU;

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
                vector < VectorXd > orthaxis_;
                vector < double > radii_;
                vector < double > dist_;
                vector < double > pvalues_; // TODO: To be used with SEQTEST
                MatrixXd basisu_;
                bool basisu_flag_;
            public:

                MatrixXd computeRotMat( const VectorXd& vec ) const;
                MatrixXd computeRiemannianExpMap( const MatrixXd& mat ) const;
                MatrixXd computeRiemannianLogMap( const MatrixXd& mat ) const;

                double computeGeodesicMeanS1( const VectorXd& angles );

                PNS( const MatrixXd& data, const unsigned int flag = 2 ) : flag_( static_cast<itype>( flag ) ), data_( data ), basisu_flag_(false) { } ;
                ~PNS(){};

                double LMsphereFit( const MatrixXd& data, const itype flag, VectorXd& x ) const;
                // internal objective function to be used inside getSubSphere
                double objFn( const VectorXd& center );
                double computeSubSphere( VectorXd& center);
                void compute(); // this corresponds to PNSmain.m
        };


} // namespace statismo

#include "PNS.hxx"

#endif /* __PNS_H_ */
