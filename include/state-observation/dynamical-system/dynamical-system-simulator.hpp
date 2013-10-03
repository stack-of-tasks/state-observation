#ifndef STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
#define STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H

#include <map>

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

        virtual void setInput( const Vector & u, unsigned k);

        virtual Vector getCurrentState() const;

        virtual unsigned getCurrentTime() const;

        virtual void interateDynamics();

        virtual void simulateDynamicsTo(unsigned k);

        virtual Vector getInput(unsigned k) const;

        virtual Vector getMeasurement( unsigned k );

        virtual Vector getState( unsigned k );

        virtual DiscreteTimeArray getMeasurementArray
                    (unsigned startingTime, unsigned duration);

        virtual DiscreteTimeArray getStateArray
                    (unsigned startingTime, unsigned duration);

        virtual void resetDynamics();

        virtual void resetSimulator();

    protected:
        DynamicalSystemFunctorBase * f_;

        DiscreteTimeArray x_;

        DiscreteTimeArray y_;

        std::map<unsigned, Vector> u_;

    };
}
#endif // STATEOBSERVATIONDYNAMICALSYSTEMSIMULATOR_H
