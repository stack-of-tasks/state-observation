#include <iostream>

#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>

using namespace stateObservation;

int test()
{
    AccelerometerGyrometer imu;

    Vector x=imu.stateVectorZero();

    x[0]=1;

    GaussianWhiteNoise vk(imu.getMeasurementSize());


    imu.setState(x,0);

    imu.setNoise(&vk);

    std::cout<<"State \n"<<x<<std::endl<<std::endl;


    //x=Vector::Random(x.rows(),1);

    std::cout<<"Measurement \n"<<imu.getMeasurements()<<std::endl;
    return 0;
}

int main()
{


    return test();

}
