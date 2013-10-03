#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>

namespace stateObservation
{
    DynamicalSystemSimulator::DynamicalSystemSimulator():
            p_(0)
    {
        //ctor
    }

    DynamicalSystemSimulator::~DynamicalSystemSimulator()
    {
        //dtor
    }

    void DynamicalSystemSimulator::setDynamicsFunctor(DynamicalSystemFunctorBase * )
    {
    }

    void DynamicalSystemSimulator::setState( const Vector & x, unsigned k)
    {
    }

    void DynamicalSystemSimulator::setInput( const Vector & x, unsigned k)
    {
    }

    Vector DynamicalSystemSimulator::getInput( unsigned k )
    {
    }

    std::vector<Vector> DynamicalSystemSimulator::getMeasurementVector
                    (unsigned startingTime, unsigned duration)
    {
    }

    void DynamicalSystemSimulator::setMeasurementNoise(NoiseBase * )
    {
    }

    void DynamicalSystemSimulator::setProcessNoise(NoiseBase *)
    {
    }


}
