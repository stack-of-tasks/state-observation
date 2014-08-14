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

    int t, i;
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

    // Contacts
    imu.setContactsNumber(2);
    uk.segment(input::contacts,6) << 0.0094904630937003645, -0.095000000000000001, 1.9819700013135044e-07,0,0,0;
    uk.segment(input::contacts+6,6) << 0.0094904630936998632, 0.095000000000000001, 1.9819700018686159e-07,0,0,0;

    // Input init
    uk.segment(input::posCom,3) <<  0,
                                    0,
                                    0.80771;//hrp2::H*0.5;

    uk.segment(input::accCom,3) <<  0,
                                    0,
                                    -9.8;

//    uk <<   -0.0135671,
//            0.001536,
//            0.80771,
//            -1.04704e-06,
//            -4.33367e-09,
//            2.27103e-08,
//            -1.04708e-06,
//            -4.34023e-09,
//            -9.8,//2.28003e-08,
//            48,//9.59632,// //DENUT INERTIE
//            48,//8.39063,//
//            2.87,//1.73947,//
//            -0.00189113,
//            0.150994,
//            -0.0378544,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            -0.0980002,
//            -7.76242e-10,
//            1.1174,
//            1.71227e-21,
//            -6.05385e-20,
//            4.42689e-21,
//            -1.18997e-06,
//            -4.94785e-09,
//            -1.9513e-18,
//            8.67051e-21,
//            -3.1192e-19,
//            2.21891e-20,
//            1.28561e-05,
//            5.35406e-08,
//            2.03768e-17,
//            0.00949046,
//            -0.095,
//            1.98197e-07,
//            6.05496e-22,
//            -0,
//            5.3869e-21,
//            0.00949046,
//            0.095,
//            1.98197e-07,
//            5.54576e-22,
//            -1.01427e-16,
//            -4.87717e-21;

//    uk <<   0.0135673,
//            0.001536,
//            0.80771,
//            -2.63605e-06,
//            -1.09258e-08,
//            5.71759e-08,
//            -1.04708e-06, //2.71345,
//            -4.34023e-09, //0.3072,
//            -9.8,//161.542,
//            33.5272,
//            32.3215,
//            1.73947,
//            -0.00189113,
//            0.150993,
//            -0.0378544,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            -0.098,
//            -6.23712e-11,
//            1.1174,
//            1.56933e-22,
//            -5.40778e-21,
//            3.86235e-22,
//            -2.99589e-06,
//            -1.24742e-08,
//            -4.76329e-18,
//            3.13865e-20,
//            -1.08156e-18,
//            7.72471e-20,
//            -0.000299589,
//            -1.24742e-06,
//            -4.76329e-16,
//            0.00949046,
//            -0.095,
//            1.98197e-07,
//            8.75907e-23,
//            -0,
//            7.17001e-22,
//            0.00949046,
//            0.095,
//            1.98197e-07,
//            6.29302e-23,
//            -1.01427e-16,
//            -6.61701e-22;

//    uk <<   0.013567,
//            0.001536,
//            0.80771,
//            5.37764e-16,
//            -1.0842e-17,
//            -1.12797e-17,
//            3.46945e-16,
//            0,
//            0,
//            48.1348,
//            46.9498,
//            1.76068,
//            -0.0863332,
//            -0.594859,
//            -0.0402246,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            0,
//            -0.0980003,
//            -1.24554e-09,
//            1.1174,
//            -1.36435e-29,
//            -4.23861e-28,
//            2.80599e-29,
//            1.76985e-17,
//            -1.21597e-17,
//            -8.04986e-26,
//            3.86788e-29,
//            1.58758e-29,
//            -2.77334e-30,
//            1.97152e-15,
//            5.15681e-21,
//            3.22928e-27,
//            0.00949046,
//            -0.095,
//            -0.0014998,
//            0,
//            5.91866e-16,
//            1.56868e-24,
//            0.00949046,
//            0.095,
//            1.98197e-07,
//            0,
//            5.46736e-16,
//            1.52532e-24;

    uk <<   0.013067,
            0.00303608,
            0.807721,
            -0.100001,
            0.300014,
            0.00234917,
            -0.100004,
            0.300016,
            0.00226504,
            48.1365,
            46.9496,
            1.7603,
            -0.168322,
            -0.56512,
            -0.0493796,
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
            -0.0985685,
            0.00171151,
            1.1174,
            -4.15862e-15,
            -2.00491e-16,
            -2.2381e-16,
            -0.113649,
            0.342301,
            -1.79819e-13,
            -8.31724e-13,
            -4.00981e-14,
            -4.47619e-14,
            -11.3649,
            34.2301,
            -1.79819e-11,
            0.00949058,
            -0.0950002,
            -0.00149675,
            -2.16874e-15,
            -6.74997e-16,
            7.85277e-17,
            0.00949035,
            0.0950002,
            3.23101e-06,
            -2.1585e-15,
            8.01036e-16,
            1.7284e-17;


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
                                -9.8*hrp2::m / (2*hrp2::linKe); // flexibility excited by the weight

//    // To put a well interpretable excitation despite of a rotated axe for the abnkles.
//    ori    <<   -2*PI/180,
//                2*PI/180,
//                0;

    // h=rotationMatrixFromContactsPositionKine(uk.segment(input::LFootPos,3),uk.segment(input::RFootPos,3),R0);
    //x0.segment(kine::ori,3) << R0.transpose()*ori;

//    x0.segment(kine::ori,3) << ori;

//    x0 <<   -0.000986063,
//            0.000261362,
//            0.00808331,
//            0.00057878,
//            0.00228587,
//            -9.98593e-05,
//            -0.0232814,
//            0.00821096,
//            0.366858,
//            0.0246247,
//            0.104175,
//            -0.00356503,
//            -0.859072,
//            0.191626,
//            -30.3446,
//            0.906266,
//            3.486,
//            -0.124809;

    x0 <<   4.43601e-05,
            -8.16189e-07,
            -0.00653148,
            -0.0011089,
            0.0457922,
            0.000198551,
            -4.41254e-08,
            3.18629e-09,
            -1.74849e-07,
            9.83174e-08,
            -1.85835e-05,
            -6.63314e-08,
            -2.1157e-08,
            7.81833e-08,
            -2.09045e-07,
            -7.88384e-06,
            -2.18611e-05,
            2.30185e-07;

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





