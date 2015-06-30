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

namespace statismo {

    template <typename T>
        class PNSModelBuilder : public ModelBuilder<T> {


            public:

};

} // namespace statismo

#include "PNSModelBuilder.hxx"

#endif /* __PNSMODELBUILDER_H_ */
