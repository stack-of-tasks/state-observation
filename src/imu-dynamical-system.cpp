#include <state-observation/dynamical-system/imu-dynamical-system.hpp>


#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{

    IMUDynamicalSystem::IMUDynamicalSystem()
    :processNoise_(0x0),dt_(1),orientationVector_(Vector3::Zero()),
        quaternion_(Quaternion::Identity())
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
       std::cout<<std::endl<<"IMUFixedContactDynamicalSystem Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR
        //ctor
    }

    IMUDynamicalSystem::~IMUDynamicalSystem()
    {
        //dtor
    }

    Vector IMUDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& u, unsigned)
    {
        assertStateVector_(x);
        assertInputVector_(u);

        Vector3 position=x.segment(kine::pos,3);
        Vector3 velocity=x.segment(kine::linVel,3);
        Vector3 acceleration=x.segment(kine::linAcc,3);

        Vector3 orientationV=x.segment(kine::ori,3);
        Vector3 angularVelocity=x.segment(kine::angVel,3);
        Vector3 angularAcceleration=x.segment(kine::angAcc,3);

        Quaternion orientation=computeQuaternion_(orientationV);

        integrateKinematics
                (position, velocity, acceleration, orientation, angularVelocity,
                        angularAcceleration, dt_);

        //x_{k+1}
        Vector xk1=Vector::Zero(18,1);

        xk1.segment(kine::pos,3) = position;
        xk1.segment(kine::linVel,3) = velocity;

        AngleAxis orientationAA(orientation);

        orientationV=orientationAA.angle()*orientationAA.axis();

        xk1.segment(kine::ori,3) =  orientationV;
        xk1.segment(kine::angVel,3) = angularVelocity;


        //inputs
        Vector3 accelerationInput =u.head(3);
        Vector3 angularAccelerationInput =u.tail(3);

        xk1.segment(kine::linAcc,3)+=accelerationInput;
        xk1.segment(kine::angAcc,3)+=angularAccelerationInput;

        if (processNoise_!=0x0)
            return processNoise_->addNoise(xk1);
        else
            return xk1;

    }

    Quaternion IMUDynamicalSystem::computeQuaternion_(const Vector3 & x)
    {
        if (orientationVector_!=x)
        {
            quaternion_ = kine::rotationVectorToAngleAxis(x);
            orientationVector_=x;
        }

        return quaternion_;
    }

    Vector IMUDynamicalSystem::measureDynamics (const Vector& x, const Vector& , unsigned k)
    {
        assertStateVector_(x);

        Vector3 acceleration=x.segment(kine::linAcc,3);

        Vector3 orientationV=x.segment(kine::ori,3);
        Vector3 angularVelocity=x.segment(kine::angVel,3);

        Quaternion q=computeQuaternion_(orientationV);

        Vector v=Vector::Zero(10,1);

        v[0]=q.w();
        v[1]=q.x();
        v[2]=q.y();
        v[3]=q.z();

        v.segment(4,3)=acceleration;
        v.tail(3)=angularVelocity;

        sensor_.setState(v,k);

        return sensor_.getMeasurements();
    }

    void IMUDynamicalSystem::setProcessNoise( NoiseBase * n)
    {
        processNoise_=n;
    }

    void IMUDynamicalSystem::resetProcessNoise()
    {
        processNoise_=0x0;
    }

    void IMUDynamicalSystem::setMeasurementNoise( NoiseBase * n)
    {
        sensor_.setNoise(n);
    }
    void IMUDynamicalSystem::resetMeasurementNoise()
    {
        sensor_.resetNoise();
    }

    void IMUDynamicalSystem::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    unsigned IMUDynamicalSystem::getStateSize() const
    {
        return stateSize_;
    }

    unsigned IMUDynamicalSystem::getInputSize() const
    {
        return inputSize_;
    }

    unsigned IMUDynamicalSystem::getMeasurementSize() const
    {
        return measurementSize_;
    }

    NoiseBase * IMUDynamicalSystem::getProcessNoise() const
    {
        return processNoise_;
    }

    NoiseBase * IMUDynamicalSystem::getMeasurementNoise() const
    {
        return sensor_.getNoise();
    }
}
