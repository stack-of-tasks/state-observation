#include <iostream>
#include <fstream>

#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/gaussian-white-noise.hpp>
#include <state-observation/dynamical-system/imu-dynamical-system.hpp>
#include <state-observation/dynamical-system/dynamical-system-simulator.hpp>
#include <state-observation/observer/extended-kalman-filter.hpp>

using namespace stateObservation;

int test()
{
    const unsigned kmax=3000;
    const double dt=1e-3;

    IMUDynamicalSystem imu;

    DynamicalSystemSimulator sym;

    //GaussianWhiteNoise vk(imu.getMeasurementSize());

    sym.setDynamicsFunctor(&imu);

    Vector x0=Vector::Zero(imu.getStateSize(),1);

    sym.setState(x0,0);

    DiscreteTimeArray u;

    for (int i=0;i<kmax/10;++i)
    {
        Vector uk=Vector::Zero(imu.getInputSize(),1);

        uk[0]=0.4 * sin(M_PI/10*i);
        uk[1]=0.6 * sin(M_PI/12*i);
        uk[2]=0.2 * sin(M_PI/5*i);

        uk[3]=10  * sin(M_PI/12*i);
        uk[4]=0.07  * sin(M_PI/15*i);
        uk[5]=0.05 * sin(M_PI/5*i);

        for (int j=0;j<10;++j)
        {
            u.pushBack(uk,i*10+j);
        }
        sym.setInput(uk,10*i);

    }

    imu.setSamplingPeriod(dt);

    sym.simulateDynamicsTo(kmax+1);

    DiscreteTimeArray y = sym.getMeasurementArray(1,kmax);
    DiscreteTimeArray x = sym.getStateArray(1,kmax);

    std::ofstream f;

    f.open("trajectory.dat");

    IMUDynamicalSystem imuFunctor;

    imuFunctor.setSamplingPeriod(dt);

    imuFunctor=imu;

    ExtendedKalmanFilter filter(imuFunctor.getStateSize(), imuFunctor.getMeasurementSize(),
                                imuFunctor.getMeasurementSize(), false);

    filter.setFunctor(& imuFunctor);

    Vector xh0=filter.stateVectorRandom()*10;

    filter.setState(xh0,0);

    Matrix p=filter.getPmatrixZero();

    for (unsigned i=0;i<filter.getStateSize();++i)
    {
        p(i,i)=xh0[i];
    }
    p=p*p.transpose();

    filter.setStateCovariance(p);

    Matrix r1=filter.getRmatrixIdentity()*0.0001;

    Matrix q1=filter.getQmatrixZero();

    Matrix r=r1*r1.transpose();
    Matrix q=q1*q1.transpose();

    ///set the covariance matrices for the extended Kalman filter
    filter.setR(r);
    filter.setQ(q);

    ///set initial input
    filter.setInput(u[0],0);

    ///set the derivation step for the finite difference method
    Vector dx=filter.stateVectorConstant(1)*1e-8;


    DiscreteTimeArray xh;
    xh.pushBack(xh0,0);

    for (int i=y.getFirstTime();i<=y.getLastTime();++i)
    {
        filter.setMeasurement(y[i],i);
        if (i<y.getLastTime())
            filter.setInput(u[i],i);

        Matrix a=filter.getAMatrixFD(dx);
        Matrix c= filter.getCMatrixFD(dx);

        //std::cout << a << std::endl << c << std::endl;

        filter.setA(a);
        filter.setC(c);

        xh.pushBack(filter.getEstimateState(i),i);

        Vector3 g;
        {
            Matrix3 R;
            Vector3 orientationV=Vector(x[i]).segment(9,3);
            double angle=orientationV.norm();
            if (angle > cst::epsilonAngle)
                R = AngleAxis(angle, orientationV/angle).toRotationMatrix();
            else
                R = Matrix3::Identity();
            g=R.transpose()*Vector3::UnitZ();
        }

        Vector3 gh;
        {
            Matrix3 Rh;
            Vector3 orientationV=Vector(xh[i]).segment(9,3);
            double angle=orientationV.norm();
            if (angle > cst::epsilonAngle)
                Rh = AngleAxis(angle, orientationV/angle).toRotationMatrix();
            else
                Rh = Matrix3::Identity();
            gh=Rh.transpose()*Vector3::UnitZ();
        }


        f << i<< " \t "<< g.transpose()-gh.transpose()<< " \t\t\t "
        << g.transpose() << " \t\t\t " << gh.transpose() << std::endl;

        //std::getchar();
    }






    return 0;
}

int main()
{


    return test();

}
