#include <state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{
    using namespace stateObservation;

    IMUFixedContactDynamicalSystem::
                    IMUFixedContactDynamicalSystem():
        processNoise_(0x0), dt_(1),orientationVector_(Vector3::Zero()),
        quaternion_(Quaternion::Identity()),
        measurementSize_(measurementSizeBase_)
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
        std::cout<<std::endl<<"IMUFixedContactDynamicalSystem Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR
    }

    IMUFixedContactDynamicalSystem::
                    ~IMUFixedContactDynamicalSystem()
    {
        //dtor
    }


    Vector IMUFixedContactDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& , unsigned)
    {
        assertStateVector_(x);

        Vector3 positionFlex=x.head(3);
        Vector3 velocityFlex=x.segment(3,3);
        Vector3 accelerationFlex=x.segment(6,3);

        Vector3 orientationFlexV=x.segment(9,3);
        Vector3 angularVelocityFlex=x.segment(12,3);
        Vector3 angularAccelerationFlex=x.tail(3);

        Quaternion orientationFlex=computeQuaternion_(orientationFlexV);

        integrateKinematics
                (positionFlex, velocityFlex, accelerationFlex, orientationFlex,
                 angularVelocityFlex, angularAccelerationFlex, dt_);

        //x_{k+1}
        Vector xk1(x);

        xk1.head(3) = positionFlex;
        xk1.segment(3,3) = velocityFlex;

        AngleAxis orientationAA(orientationFlex);
        orientationFlexV=orientationAA.angle()*orientationAA.axis();

        xk1.segment(9,3) =  orientationFlexV;
        xk1.segment(12,3) = angularVelocityFlex;

        if (processNoise_!=0x0)
            return processNoise_->addNoise(xk1);
        else
            return xk1;
    }

    Quaternion IMUFixedContactDynamicalSystem::computeQuaternion_
                (const Vector3 & x)
    {
        if (orientationVector_!=x)
        {
            orientationVector_ = x;
            quaternion_ = tools::rotationVectorToAngleAxis(x);
        }

        return quaternion_;
    }

    Vector IMUFixedContactDynamicalSystem::measureDynamics
                (const Vector& x, const Vector& u, unsigned k)
    {
        assertStateVector_(x);

        Vector3 positionFlex=x.head(3);
        Vector3 velocityFlex=x.segment(3,3);
        Vector3 accelerationFlex=x.segment(6,3);

        Vector3 orientationFlexV=x.segment(9,3);
        Vector3 angularVelocityFlex=x.segment(12,3);
        Vector3 angularAccelerationFlex=x.tail(3);

        Quaternion qFlex = computeQuaternion_(orientationFlexV);
        Matrix3 rFlex = qFlex.toRotationMatrix();


        assertInputVector_(u);

        Vector3 positionControl=u.head(3);
        Vector3 velocityControl=u.segment(3,3);
        Vector3 accelerationControl=u.segment(6,3);

        Vector3 orientationControlV=u.segment(9,3);
        Vector3 angularVelocityControl=u.segment(12,3);

        Quaternion qControl=computeQuaternion_(orientationControlV);

        Quaternion q = qFlex * qControl;

        Vector3 acceleration =
         (tools::skewSymmetric(angularAccelerationFlex)
              + tools::square(tools::skewSymmetric(angularVelocityFlex)))
                  * rFlex * positionControl
         + 2*tools::skewSymmetric(angularVelocityFlex) * rFlex * velocityControl
         + accelerationFlex + rFlex * accelerationControl;

        Vector3 angularVelocity = angularVelocityFlex +
                                    rFlex * angularVelocityControl;

        Vector v=Vector::Zero(10,1);

        v[0]=q.w();
        v[1]=q.x();
        v[2]=q.y();
        v[3]=q.z();

        v.segment(4,3)=acceleration;
        v.tail(3)=angularVelocity;

        sensor_.setState(v,k);

        Vector y = Matrix::Zero(measurementSize_,1);

        y.head(sensor_.getMeasurementSize()) = sensor_.getMeasurements();

        for (unsigned i=0; i<contactPositions_.size();++i)
        {
            y.segment(sensor_.getMeasurementSize()+i*3,3)=
                rFlex * contactPositions_[i] + positionFlex;
        }

        return y;
    }

    void IMUFixedContactDynamicalSystem::setProcessNoise(NoiseBase * n)
    {
        processNoise_=n;
    }
    void IMUFixedContactDynamicalSystem::resetProcessNoise()
    {
        processNoise_=0x0;
    }

    void IMUFixedContactDynamicalSystem::setMeasurementNoise
                ( NoiseBase * n)
    {
        sensor_.setNoise(n);
    }

    void IMUFixedContactDynamicalSystem::resetMeasurementNoise()
    {
        sensor_.resetNoise();
    }

    void IMUFixedContactDynamicalSystem::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    unsigned IMUFixedContactDynamicalSystem::getStateSize()
    {
        return stateSize_;
    }

    unsigned IMUFixedContactDynamicalSystem::getInputSize()
    {
        return inputSize_;
    }

    unsigned IMUFixedContactDynamicalSystem::getMeasurementSize()
    {
        return measurementSize_;
    }

    NoiseBase * IMUFixedContactDynamicalSystem::getProcessNoise() const
    {
        return processNoise_;
    }

    NoiseBase * IMUFixedContactDynamicalSystem::
                                                    getMeasurementNoise() const
    {
        return sensor_.getNoise();
    }

    void IMUFixedContactDynamicalSystem::setContactsNumber(unsigned i)
    {
        measurementSize_ = measurementSizeBase_ + 3 * i;
        contactPositions_.resize(i, Vector3::Zero());
    }

    void IMUFixedContactDynamicalSystem::setContactPosition
                                        (unsigned i, const Vector3 & position)
    {
        BOOST_ASSERT( i< contactPositions_.size() &&
                    "ERROR: The index of contact is out of range.");

        contactPositions_[i] = position;
    }
}
}
