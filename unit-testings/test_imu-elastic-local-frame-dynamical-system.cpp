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
double h, thetax, thetay, thetaz;
// Definition of the length and the unit vector between the two foots.
Vrl=vLFootPos-vRFootPos;
j=kine::rotationVectorToAngleAxis(Vrl);
h = j.angle(); // length between the ankles
axis = j.axis(); // unit vector between the ankles
//Definition of the transformation (rotation) between (x,y,z) and (perpendicular of j, j, z).
thetax=atan2(axis[1],axis[2]);
//theta=theta(0,0,atan(axis[0]/axis[1]));
thetay = atan2(axis[0],axis[2]);//theta[1];
thetaz = atan2(axis[0],axis[1]);//theta[2];
R << cos(thetay)*cos(thetaz), -cos(thetay)*sin(thetaz), sin(thetay),
sin(thetaz), cos(thetaz), 0,
-sin(thetay)*cos(thetaz), sin(thetay)*sin(thetaz), cos(thetay);
return h;
}

int test()
{

    int t, i;
    t = time(NULL);
    srand(t);

    Matrix3 R0;
    Vector3 Vrl0, axis0, ori, theta0;
    AngleAxis  j0;
    double h0;

    /// The number of samples
    const unsigned kmax=2000;

    ///sampling period
    const double dt=5e-3;

    ///Sizes of the states for the state, the measurement, and the input vector
    const unsigned stateSize=18;
    const unsigned measurementSize=6;
    const unsigned inputSize=48;


    ///The array containing all the states, the measurements and the inputs
    IndexedMatrixArray x;
    IndexedMatrixArray u;

    ///simulation of the signal
    /// the IMU dynamical system functor
    flexibilityEstimation::IMUElasticLocalFrameDynamicalSystem imu(dt);
    imu.setInputSize(inputSize);

    ///the simulator initalization
    DynamicalSystemSimulator sim;
    sim.setDynamicsFunctor(& imu); // choice of the dynamical system

    double h;
//std::cout << "hello" << std::endl;
    /// Input initialization

    Vector6 Inert;
    Vector uk=Vector::Zero(imu.getInputSize(),1);

//    // Inertia, the robot is considered as a full cylindar
//    Inert << 0.25*hrp2::m*hrp2::R*hrp2::R+0.33*hrp2::m*hrp2::H*hrp2::H, 0.25*hrp2::m*hrp2::R*hrp2::R+0.33*hrp2::m*hrp2::H*hrp2::H, 0.5*hrp2::m*hrp2::R*hrp2::R, 0.0, 0.0, 0.0;
//    uk.segment(input::Inertia,6)=Inert;
//
//    // Contacts
//    imu.setContactsNumber(1);
//    uk.segment(input::contacts,6) << 0.0094904630937003645, -0.095000000000000001, 1.9819700013135044e-07,0,0,0;
//    uk.segment(input::contacts+6,6) << 0.0094904630936998632, 0.095000000000000001, 1.9819700018686159e-07,0,0,0;
//
//    // Input init
//    uk.segment(input::posCom,3) <<  0,
//                                    0,
//                                    0.80771;//hrp2::H*0.5;
//
//    uk.segment(input::accCom,3) <<  0,
//                                    0,
//                                    -9.8;

    uk <<    0.0135672,
            0.001536,
            0.80771,
            -2.50425e-06,
            -1.03787e-08,
            5.4317e-08,
            -2.50434e-06,
            -1.03944e-08,
            5.45321e-08,
            48.1348,
            46.9498,
            1.76068,
            -0.0863332,
            -0.59487,
            -0.0402246,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -0.098,
            -1.21619e-10,
            1.1174,
            3.06752e-22,
            -1.06094e-20,
            7.75345e-22,
            -2.84609e-06,
            -1.18496e-08,
            -4.52691e-18,
            2.95535e-20,
            -1.0346e-18,
            7.58731e-20,
            -0.000284609,
            -1.18496e-06,
            -4.52691e-16,
            0.0094904630937003645,
            0,//-0.15000000000000001,
            1.9819700013135044e-07,
            0,
            0,
            0;
//            0.0094904630936998632,
//            0.15000000000000001,
//            1.9819700018686159e-07,
//            0,
//            0,
//            0;

//    // contacts position
//    double angleZ;
//    angleZ=-PI/4;
//
//    Vector6 v1, v2;
//    v1 <<   -0.095*tan(angleZ),//0.0094904630936998632,
//            0.095000000000000001,
//            1.9819700018686159e-07,
//            0,
//            0,
//            angleZ;
//
//    v2 <<    0.095*tan(angleZ),//0.0094904630936998632+
//            -0.095000000000000001,
//            1.9819700018686159e-07,
//            0,
//            0,
//            angleZ;
//
//    Matrix4 cont1, cont2;
//    cont1=kine::vector6ToHomogeneousMatrix(v1);
//    cont2=kine::vector6ToHomogeneousMatrix(v2);
//    uk.segment(input::contacts,6) << cont1(0,3), cont1(1,3), cont1(2,3), cont1(1,2), cont1(0,2), cont1(0,1);//0.0094904630937003645, -0.095000000000000001, 1.9819700013135044e-07,0,0,0;
//    uk.segment(input::contacts+6,6) << cont2(0,3), cont2(1,3), cont2(2,3), cont2(1,2), cont2(0,2), cont2(0,1);//0.0094904630936998632, 0.095000000000000001, 1.9819700018686159e-07,0,0,0;
//    std::cout << uk.transpose() << std::endl;


//    for (i=0;i<uk.size();++i)
//    {
//        if(uk(i)<0.001 && uk(i)>-0.001)
//        {
//            uk(i)=0.0;
//        }
//    }


    /// State initialization
    Vector x0=(Vector::Zero(stateSize,1));
    x0.segment(kine::pos,3) <<  0,
                                0,
                                0;//-9.8*hrp2::m / (2*hrp2::linKe); // flexibility excited by the weight

    // To put a well interpretable excitation despite of a rotated axe for the abnkles.


std::cout << "hello" << std::endl;
//    ori    <<   0;//-2*PI/180,
//                0;//2*PI/180,
//                0;

//    h=rotationMatrixFromContactsPositionKine(uk.segment(input::contacts,6),uk.segment(input::contacts+6,6),R0);
//    x0.segment(kine::ori,3) << R0.transpose()*ori;

//    x0.segment(kine::ori,3) << ori;

//    for (i=0;i<x0.size();++i)
//    {
//        if(x0(i)<0.1 && x0(i)>-0.1)
//        {
//            x0(i)=0.0;
//        }
//    }

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

    //imu.test(); // test unitaires
}

int main()
{


    return test();

}





