#include <iostream>
#include <fstream>

#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/examples/offline-ekf-flexibility-estimation.hpp>
#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <time.h>

//Tracking NaN
#include <fenv.h>


using namespace stateObservation;

    double rotationMatrixFromContactsPositionKine(const Vector3 vLFootPos, const Vector3 vRFootPos, Matrix3& R )
    {

        Vector3 Vrl, axis, theta;
        AngleAxis j;
        double h, thetay, thetaz;

        // Definition of the length and the unit vector between the two foots.
        Vrl=vLFootPos-vRFootPos;
        j=kine::rotationVectorToAngleAxis(Vrl);
        h = j.angle(); // length between the ankles
        axis = j.axis(); // unit vector between the ankles

         //Definition of the transformation (rotation) between (x,y,z) and (perpendicular of j, j, z).
        theta=kine::unitVectorToRotationVector(axis);
        thetay = theta[1];
        thetaz = theta[2];


        R <<    cos(thetay)*cos(thetaz), -cos(thetay)*sin(thetaz), sin(thetay),
                sin(thetaz), cos(thetaz), 0,
                -sin(thetay)*cos(thetaz), sin(thetay)*sin(thetaz), cos(thetay);

        return h;

    }


int test()
{

    int t;
    t = time(NULL);
    srand(t);

    Matrix3 R0;
    Vector3 Vrl0, axis0, ori, theta0;
    AngleAxis  j0;
    double h0;

    //feenableexcept(FE_OVERFLOW); // FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW,  FE_ALL_EXCEPT
    //feenableexcept(FE_DIVBYZERO || FE_INVALID || FE_OVERFLOW || FE_UNDERFLOW || FE_INEXACT); // FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW,  FE_ALL_EXCEPT
    //fedisableexcept(FE_ALL_EXCEPT);
    //std::cout << Tc << "\n\n" << std::endl;
    //if(Tc!=Tc){ std::cout << "totoNaN\n\n" << std::endl; }

    /// The number of samples
    const unsigned kmax=2000;

    ///sampling period
    const double dt=5e-3;

    ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=54;


    ///The array containing all the states, the measurements and the inputs
    DiscreteTimeArray x;
    DiscreteTimeArray u;

    ///simulation of the signal
    /// the IMU dynamical system functor
    flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem imu(dt);
    imu.setInputSize(inputSize);

    ///the simulator initalization
    DynamicalSystemSimulator sim;
    sim.setDynamicsFunctor(& imu); // choice of the dynamical system

    double h;

    /// Input initialization

    Vector6 Inert;
    Vector uk=Vector::Zero(imu.getInputSize(),1);

    // Inertia, the robot is considered as a full cylindar
    Inert << 0.25*hrp2::m*hrp2::R*hrp2::R+0.33*hrp2::m*hrp2::H*hrp2::H, 0.25*hrp2::m*hrp2::R*hrp2::R+0.33*hrp2::m*hrp2::H*hrp2::H, 0.5*hrp2::m*hrp2::R*hrp2::R, 0.0, 0.0, 0.0;
    uk.segment(input::Inertia,6)=Inert;

    // std::cout << "Inertia" << Inert << "\n\n" << std::endl;

    // Position of the contacts
//    uk.segment(input::LFootPos,3) <<    0,
//                                        0.155,
//                                        0; // Z compoent have to be 0
    uk.segment(input::contacts,6) << 0.0094904630937003645, -0.095000000000000001, 1.9819700013135044e-07,0,0,0;

//    uk.segment(input::RFootPos,3) <<    0,
//                                        -0.155,
//                                        0; // Z compoent have to be 0
    uk.segment(input::contacts+3,6) << 0.0094904630936998632, 0.095000000000000001, 1.9819700018686159e-07,0,0,0;

    imu.setContactsNumber(2);
    //imu.setContactPosition(0,uk.segment(input::contacts,3));
    //imu.setContactPosition(1,uk.segment(input::contacts+3,3));


    // Linear position of the Com IN THE LOCAL FRAME (here for a full cylindar)
    uk.segment(input::posCom,3) <<  0,
                                    0,
                                    hrp2::H*0.5;

    /// State initialization

    Vector x0=(Vector::Zero(stateSize,1));

    x0.segment(kine::pos,3) <<  0,
                                0,
                                -9.8*hrp2::m/hrp2::linKe; // flexibility excited by the weight

    // To put a well interpretable excitation despite of a rotated axe for the abnkles.
    ori    <<   0,//-2*PI/180,
                2*PI/180,
                0;
   // h=rotationMatrixFromContactsPositionKine(uk.segment(input::LFootPos,3),uk.segment(input::RFootPos,3),R0);
    //x0.segment(kine::ori,3) << R0.transpose()*ori;

    x0.segment(kine::ori,3) << ori;

    sim.setState(x0,0);

    for (int k=0;k<10;++k)
    {
        u.setValue(uk,k);

        ///give the input to the simulator
        ///we only need to give one value and the
        ///simulator takes automatically the appropriate value
        sim.setInput(uk,k);
    }

    ///set the sampling perdiod to the functor
    imu.setSamplingPeriod(dt);


    sim.simulateDynamicsTo(kmax);

    x = sim.getStateArray(1,kmax);


    x.writeInFile("state.dat");
    u.writeInFile("input.dat");
    //std::cout <<x[kmax].norm ()<< " "<< x[kmax].transpose() << std::endl;

    // imu.test(x0); // test unitaires


}

int main()
{


    return test();

}





