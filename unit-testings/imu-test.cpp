#include <iostream>

#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/dynamical-system/imu-dynamical-system.hpp>
#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>

using namespace stateObservation;

int test()
{
    IMUDynamicalSystem imu;

    DynamicalSystemSimulator sym;

    //GaussianWhiteNoise vk(imu.getMeasurementSize());

    sym.setDynamicsFunctor(&imu);

    Vector x=Vector::Zero(imu.getStateSize(),1);

    sym.setState(x,0);

    std::vector<Vector> u;

    for (int i=0;i<50;++i)
    {
        Vector uk=Vector::Random(imu.getInputSize(),1);
        for (int j=0;j<10;++j)
        {
            u.push_back(uk);
        }
        sym.setInput(uk,10*i);

    }

    sym.simulateDynamicsTo(501);

    DiscreteTimeArray y = sym.getMeasurementArray(1,500);








    return 0;
}

int main()
{


    return test();

}
