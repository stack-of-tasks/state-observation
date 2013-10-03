#ifndef STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
#define STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H

#include <vector>

#include <state-observation/noise/noise-base.hpp>

#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>

namespace stateObservation
{

    class DynamicalSystemSimulator
    {
    public:
        DynamicalSystemSimulator();
        virtual ~DynamicalSystemSimulator();

        virtual void setDynamicsFunctor(DynamicalSystemFunctorBase * );

        virtual void setState( const Vector & x, unsigned k);

        virtual void setInput( const Vector & x, unsigned k);

        virtual Vector getInput( unsigned k );

        virtual std::vector<Vector> getMeasurementVector
                    (unsigned startingTime, unsigned duration);

        virtual void setMeasurementNoise(NoiseBase * );

        virtual void setProcessNoise(NoiseBase *);

    protected:
        DynamicalSystemFunctorBase * p_;

        NoiseBase * measurementNoise_;

        NoiseBase * processNoise_;


    private:
    };
}
#endif // STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
