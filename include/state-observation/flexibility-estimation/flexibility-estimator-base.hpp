/**
 * \file     flexibility-estimator-base.hpp
 * \author   Mehdi Benallegue
 * \date     2013
 * \brief    Definitions of types and some structures.
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP
#define FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{
    class FlexibilityEstimatorBase
    {
    public:
        virtual ~FlexibilityEstimatorBase(){}

        explicit FlexibilityEstimatorBase
                            (double dt=0.005);

        virtual unsigned getStateSize() const =0;

        virtual unsigned getInputSize() const =0;

        virtual void setFlexibilityGuess(const Matrix &)=0;

        virtual void setMeasurement(const Vector &)=0;

        virtual Matrix4 getFlexibility()=0;

        virtual void setSamplingPeriod(double);

    protected:

        double dt_;//sampling period

    };

} //namespace flexibilityEstimation
}//namespace stateobservation

#endif //FLEXIBILITY_ESTIMATION_FLEXIBILITY_ESTIMATOR_BASE_HPP
