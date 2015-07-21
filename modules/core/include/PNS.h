#ifndef __PNS_H__
#define __PNS_H__

#include <memory>
#include <vector>
#include "CommonTypes.h"
#include "GenericFunctor.h"
#include <unsupported/Eigen/NonLinearOptimization>
// TODO: 
// 1. Note that data is stored as row major matrix
// 2. Make sure to cast as double
namespace statismo {

    template <typename T>
        class PNS{
            // class to include transformation between s and E
            //
            //
            private:
                typedef struct NestedSphere_t {
                    VectorTypeDoublePrecision axis;
                    VectorTypeDoublePrecision residual;
                    double radius;
                } NestedSphere;

                vector <NestedSphere*> NestedSpheres;
                VectorTypeDoublePrecision* getAxis( size_t pos ) const;
                VectorTypeDoublePrecision* getResidual( size_t pos ) const;
            public:
                PNS();
                ~PNS();

        };

} // namespace statismo

#include "PNS.hxx"

#endif /* __PNS_H_ */
