#include <state-observation/flexibility-estimation/flexibility-estimator-base.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{
    void FlexibilityEstimatorBase::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    FlexibilityEstimatorBase::FlexibilityEstimatorBase
            (double dt):
            dt_(dt)
    {

    }
}
}
