#ifndef __PNS_H__
#define __PNS_H__

#include <memory>
#include <vector>
using std::vector;
//#include "CommonTypes.h"
#include <unsupported/Eigen/NonLinearOptimization>
using Eigen::LevenbergMarquardt;
#include <Eigen/Dense>
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

#include "sphereResFunctor.h"
#include "PNSutil.h"
// TODO: 
// 1. Note that data is stored as row major matrix
// 2. Make sure to cast as double
// 3. Use column major matrix for computation but as for interfacing with the outside, return as a row major matrix
// 4. Make sure to check all the codes to see if there could be aliasing issue
// If so, use .eval fisrt
namespace statismo {
    template <typename T>
        class PNS{
            private:
                // enum for iteration type
                // TODO: SEQTEST to be implemented later
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
                double modBy2PI( const double& x ) { return  ( x - (2*PI)*floor( x / (2*PI) ) ); }; // Apparently compiler internally inline this function if it can see in the header.
                double computeGeodesicMeanS1( const VectorXd& angles );
            public:
                // NOTE: punted the job to convert data into column major matrix to the caller
                PNS( const MatrixXd& data, const unsigned int flag = 2 ) : flag_( staic_cast<itype>( flag ) ), data_( data ) { } ;
                ~PNS(){ };
                // TODO: Try to change magic numbers to enums I defined in fuctors
                double LMsphereFit( const MatrixXd& data, VectorXd& x, const int itype = 2);
                // internal objective function to be used inside getSubSphere
                double objFn( const VectorXd& center, const MatrixXd& data, const double r) const;
                double computeSubSphere( const MatrixXd& data, const itype flag, VectorXd& center, double& error) const;
                double getSubSphere( const MatrixXd& data, const itype flag, VectorXd& center ) const;
                MatrixXd compute(const double& alpha = 0.05, const double& R = 100); // this corresponds to PNSmain.m

                // Getters for PNS data structures
                VectorXd getOrthAxis( size_t index ) const { return orthaxis_[index]; } ;
                double getRadii( size_t index ) const { return radii_[index]; };
                double getDist( size_t index ) const { return dist_[index]; };
                double getPvalue( size_t index ) const { return pvalues_[index]; };
                MatrixXd getBasisU() const { return basisu_; };

                // PNS transformation methods
                MatrixXd computeS2E( const MatrixXd& sphereData ) const;
                MatrixXd computeE2S( const MatrixXd& euclideanData ) const;
        };


} // namespace statismo

#include "PNS.hxx"

#endif /* __PNS_H_ */
