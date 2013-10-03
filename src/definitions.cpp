#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
     DiscreteTimeMatrix::DiscreteTimeMatrix(const Matrix& v,unsigned k):
            isSet_(true),
            k_(k),
            v_(v)
    {
    }

    DiscreteTimeMatrix::DiscreteTimeMatrix():
            isSet_(false),
            k_(0)
    {
    }

}
