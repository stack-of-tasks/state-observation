#ifndef STATEOBSERVERDYNAMICALSYSTEMFUNCTORBASE_H
#define STATEOBSERVERDYNAMICALSYSTEMFUNCTORBASE_H

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{

    /**
    * \class  DynamicsFunctorBase
    * \brief
    *        This is the base class of any functor that describes the dynamics
    *        of the state and the measurement.
    *        This class is to be derived in order to be given
    *        to the Extended Kalman Filter.
    *
    */

    class DynamicalSystemFunctorBase
    {
    public:
        DynamicalSystemFunctorBase();
        virtual ~DynamicalSystemFunctorBase();

        ///The function to oberload to describe the dynamics of the state
        virtual Vector stateDynamics
        (const Vector& x, const Vector& u, unsigned k)=0;

        ///The function to overloas to describe the dynamics of the sensor (measurements)
        virtual Vector measureDynamics
        (const Vector& x, const Vector& u, unsigned k)=0;

        ///The method to overload if the functor needs to be reset when the
        ///Exteded Kalman filter is reset itself
        virtual void reset(){}
    protected:
    private:
    };
}

#endif // STATEOBSERVERDYNAMICSYSTEMFUNCTORBASE_H
