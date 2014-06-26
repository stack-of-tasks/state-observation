/*
 * dynamical-system.cpp
 *
 *  Created on: 19 mai 2014
 *      Author: alexis
 */

#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <iostream>


namespace stateObservation
{
namespace flexibilityEstimation
{
    using namespace stateObservation;



   IMUElasticLocalFrameDynamicalSystem::
    	IMUElasticLocalFrameDynamicalSystem(double dt):
        processNoise_(0x0), dt_(dt),orientationVector_(Vector3::Zero()),
        quaternion_(Quaternion::Identity()),
        measurementSize_(measurementSizeBase_)
    {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
        std::cout<<std::endl<<"IMUElasticLocalFrameDynamicalSystem Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR
    }

    IMUElasticLocalFrameDynamicalSystem::
                    ~IMUElasticLocalFrameDynamicalSystem()
    {
        //dtor
    }


    Vector3 IMUElasticLocalFrameDynamicalSystem::computeFc(const stateObservation::Vector& x){


        // Vector we want
        Vector3 Fc;
        double ke=hrp2::linKe;
        double kv=hrp2::linKv;

        // stifness ans viscosity
        Vector3 Ke(ke,ke,ke);
        Vector3 Kv(kv,kv,kv);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));

        Fc=-Ke.cwiseProduct(positionFlex)-Kv.cwiseProduct(velocityFlex);

        return Fc;

    }

    double IMUElasticLocalFrameDynamicalSystem::rotationMatrixFromCOntactsPositiontr(const Vector3 vLFootPos, const Vector3 vRFootPos, Matrix3& R )
    {

        Vector3 Vrl, axis, theta;
        AngleAxis j;
        double h, thetay, thetaz;

        // Definition of the length and the unit vector between the two foots.
        Vrl=vLFootPos-vRFootPos;
        j=kine::rotationVectorToAngleAxis(Vrl);
        h = j.angle(); // length between the ankles
        axis = j.axis(); // unit vector between the ankles

        // Definition of the transformation (rotation) between (x,y,z) and (perpendicular of j, j, z).
        theta=kine::unitVectorToRotationVector(axis);
        thetay = theta[1];
        thetaz = theta[2];


        R <<    cos(thetay)*cos(thetaz), -cos(thetay)*sin(thetaz), sin(thetay),
                sin(thetaz), cos(thetaz), 0,
                -sin(thetay)*cos(thetaz), sin(thetay)*sin(thetaz), cos(thetay);

        return h;

    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::computeTc(const stateObservation::Vector& x, const stateObservation::Vector& u)
    {

        // Isotropic stifness and viscosity for a simple contact case
        double kte=hrp2::angKe;
        double ktv=hrp2::angKv;
        double kfe=hrp2::linKe;
        double kfv=hrp2::linKv;

        // anisotropic stifness and vicosity
        Vector3 Ke;
        Vector3 Kv;

        // Vector we want
        Vector3 Tc;

        // Number of contacts
        unsigned nbContacts = getContactsNumber();

        // Flexibility state
        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlexV(x.segment(kine::angVel,3));

        Vector3 angularVelocityFlex;
        Vector3 orientationFlex;
        Matrix3 R;
        double h;
        Vector3 vLFootPos=u.segment(input::contact1Pos,3);
        Vector3 vRFootPos=u.segment(input::contact2Pos,3);

        if (nbContacts==1)
        {
            // Definittion of isotropic stiffnes and viscosity in (x,y,z).
            Ke << kte,kte,kte;
            Kv << ktv,ktv,ktv;

            // Expression of the Orientation/Velocity of the flexibility in the (x,y,z) frame.
            orientationFlex=orientationFlexV;
            angularVelocityFlex=angularVelocityFlexV;
        }
        else if (nbContacts==2)
        {
            h=rotationMatrixFromCOntactsPositiontr(vLFootPos,vRFootPos,R);

            // Definittion of anisotropic stiffnes and viscosity in (j, perpendicular of j, z) frame.
            Ke << h*h*0.5*kfe,kte,h*h*0.5*kfe; // (perpendicular of j, j, z)
            Kv << h*h*0.5*kfv,ktv,h*h*0.5*kfv; // (perpendicular of j, j, z)

            // Expression of the Orientation/Velocity of the flexibility in the (j, perpendicular of j, z).
            orientationFlex=R*orientationFlexV;
            angularVelocityFlex=R*angularVelocityFlexV;
        }

        // Determination of the Tc vector in either the (x,y,z) or (j, perpendicular of j, z) frame.
        Tc =-Ke.cwiseProduct(orientationFlex);
        Tc += -Kv.cwiseProduct(angularVelocityFlex);
        Tc += kine::skewSymmetric(positionFlex)*computeFc(x);

        // Determination of the Tc vector in the (x,y,z) frame.
        if (nbContacts==1)
        {
            return Tc;
        }
        else if (nbContacts==2)
        {
            Tc=R.transpose()*Tc;
            return Tc;
        }


    }


    void IMUElasticLocalFrameDynamicalSystem::test(const Vector& x)
    {
        Vector6 v;
        Vector3 u;

        v << 1,2,3,4,5,6;
        u << 1,2,3;

        Matrix4 m(kine::vector6ToHomogeneousMatrix(v));
        Matrix3 mat(kine::skewSymmetric(u));

        std::cout << m << std::endl;
        std::cout << mat << std::endl;
    }


    Vector3 IMUElasticLocalFrameDynamicalSystem::computeAccelerationAngular
    	(const Vector& x, const Vector& u, unsigned k)
    {
   // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW); // FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW,  FE_ALL_EXCEPT

    //fedisableexcept(FE_ALL_EXCEPT);
        // Vector we want
        Vector3 AccAngular;
    //fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);// || FE_INVALID || FE_OVERFLOW || FE_UNDERFLOW);
        // To simplify reading => to remplace by physics meaning of Mat and Vec.
        Matrix3 Mat;
        Vector3 Vec;

        // State vector
        Vector3 positionFlex(x.segment(kine::pos,3));               // t
        Vector3 velocityFlex(x.segment(kine::linVel,3));            // tdot
        Vector3 orientationFlexV(x.segment(kine::ori,3));           // Omega (position)
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));     // omega (velocity)
        // Input vector
        Vector3 positionCom(u.segment(input::posCom,3));
        Vector3 velocityCom(u.segment(input::velCom,3));
        Vector3 accelerationCom(u.segment(input::accCom,3));
        Vector3 AngMomentum(u.segment(input::AngMoment,3));
        Vector3 dotAngMomentum(u.segment(input::dotAngMoment,3));

//        std::cout << "angularVelocityFlex" << angularVelocityFlex << std::endl;

        // To be human readable
       // const Vector3 Fc(computeFc(x));

        const Vector3 Tc(computeTc(x,u));
        const Quaternion qFlex (computeQuaternion_(orientationFlexV));
        const Matrix3 R (qFlex.toRotationMatrix());
        const Matrix3 Inertia(kine::computeInertiaTensor(u.segment(input::Inertia,6)));
        const Matrix3 dotInertia(kine::computeInertiaTensor(u.segment(input::dotInertia,6)));

        Mat = R*Inertia*R.transpose();
        Mat += hrp2::m*kine::skewSymmetric(R*positionCom)*kine::skewSymmetric(R*positionCom);
        Mat = Mat.inverse().eval();

        Vec=Tc   -   (  (kine::skewSymmetric(angularVelocityFlex)*R*Inertia*R.transpose()+R*dotInertia*R.transpose())*angularVelocityFlex
                                +R*dotAngMomentum
                                +kine::skewSymmetric(angularVelocityFlex)*R*AngMomentum
                                +hrp2::m*kine::skewSymmetric(positionFlex)*
                                    (
                                        kine::skewSymmetric(angularVelocityFlex)*kine::skewSymmetric(angularVelocityFlex)*R*positionCom
                                        +2*kine::skewSymmetric(angularVelocityFlex)*R*velocityCom
                                        +R*accelerationCom
                                    )
                                +hrp2::m*kine::skewSymmetric(R*positionCom+positionFlex)*cst::gravity
                      )
                -   (kine::skewSymmetric(positionFlex)+kine::skewSymmetric(R*positionCom))*
                    (   computeFc(x)
                        -(  R*hrp2::m*accelerationCom
                            +2*kine::skewSymmetric(orientationFlexV)*R*velocityCom
                            +hrp2::m*kine::skewSymmetric(orientationFlexV)*kine::skewSymmetric(orientationFlexV)*R*positionCom
                            +hrp2::m*cst::gravity
                         )
                    );


        AccAngular = Mat*Vec;
//std::cout << "AccAng" << AccAngular << "\n\n" << std::endl;
        return AccAngular;

    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::computeAccelerationLinear
    	(const Vector& x, const Vector& u, unsigned k)
    {
        // Vector we want
        Vector3 AccLinear;

        // State vector
        Vector3 positionFlex(x.segment(kine::pos,3));               // t
        Vector3 velocityFlex(x.segment(kine::linVel,3));            // tdot
        Vector3 orientationFlexV(x.segment(kine::ori,3));           // \Omega (position)
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));     // \omega (velocity)

        // To be human readable
        const Vector3 Fc(computeFc(x));
        //const Vector3 Tc(computeTc(x,u));
        const Quaternion qFlex (computeQuaternion_(orientationFlexV));
        const Matrix3 R (qFlex.toRotationMatrix());
        const Matrix3 Inertia(kine::computeInertiaTensor(u.segment(input::Inertia,6)));
        const Matrix3 dotInertia(kine::computeInertiaTensor(u.segment(input::dotInertia,6)));

        // Input vector
        Vector3 positionCom(u.segment(input::posCom,3));
        Vector3 velocityCom(u.segment(input::velCom,3));
        Vector3 accelerationCom(u.segment(input::accCom,3));
        Vector3 AngMomentum(u.segment(input::AngMoment,3));
        Vector3 dotAngMomentum(u.segment(input::dotAngMoment,3));

        AccLinear   =   1/hrp2::m*
                            (   Fc
                                -(  R*hrp2::m*accelerationCom
                                    +2*kine::skewSymmetric(orientationFlexV)*R*velocityCom
                                    +hrp2::m*kine::skewSymmetric(orientationFlexV)*kine::skewSymmetric(orientationFlexV)*R*positionCom
                                    +hrp2::m*cst::gravity
                                 )
                            )
                        +kine::skewSymmetric(R*positionCom)*computeAccelerationAngular(x,u,k);

        return AccLinear;

    }

    Vector IMUElasticLocalFrameDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& u, unsigned k)
    {

        assertStateVector_(x);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        //Vector3 accelerationFlex(x.segment(kine::linAcc,3));

    //feenableexcept(FE_DIVBYZERO || FE_INVALID || FE_OVERFLOW || FE_UNDERFLOW); // FE_DIVBYZERO, FE_INEXACT, FE_INVALID, FE_OVERFLOW, FE_UNDERFLOW,  FE_ALL_EXCEPT
    //fedisableexcept(FE_ALL_EXCEPT);
        Vector3 accelerationFlex(computeAccelerationLinear(x, u, k));
    //fedisableexcept(FE_ALL_EXCEPT);

        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));
        //Vector3 angularAccelerationFlex(x.segment(kine::angAcc,3));
        Vector3 angularAccelerationFlex(computeAccelerationAngular(x, u, k));

        Quaternion orientationFlex(computeQuaternion_(orientationFlexV));

//        Vector3 positionControl(u.segment(kine::pos,3));
//        Vector3 velocityControl(u.segment(kine::linVel,3));
//        Vector3 accelerationControl(u.segment(kine::linAcc,3));
//
//        Vector3 orientationControlV(u.segment(kine::ori,3));
//        Vector3 angularVelocityControl(u.segment(kine::angVel,3));

        Vector3 positionControl(u.segment(input::posIMU,3));
        Vector3 velocityControl(u.segment(input::linVelIMU,3));
        Vector3 accelerationControl(u.segment(input::linAccIMU,3));

        Vector3 orientationControlV(u.segment(input::oriIMU,3));
        Vector3 angularVelocityControl(u.segment(input::angVelIMU,3));

        integrateKinematics
                (positionFlex, velocityFlex, accelerationFlex, orientationFlex,
                 angularVelocityFlex, angularAccelerationFlex, dt_);

        //x_{k+1}
        Vector xk1(x);


        xk1.segment(kine::pos,3) = positionFlex;
        xk1.segment(kine::linVel,3) = velocityFlex;
        xk1.segment(kine::linAcc,3) = accelerationFlex;

        AngleAxis orientationAA(orientationFlex);
        orientationFlexV=orientationAA.angle()*orientationAA.axis();

        xk1.segment(kine::ori,3) =  orientationFlexV;
        xk1.segment(kine::angVel,3) = angularVelocityFlex;
        xk1.segment(kine::angAcc,3) = angularAccelerationFlex;

         if (processNoise_!=0x0)
            return processNoise_->addNoise(xk1);
        else
            return xk1;



    }

    Quaternion IMUElasticLocalFrameDynamicalSystem::computeQuaternion_
                (const Vector3 & x)
    {
        if (orientationVector_!=x)
        {
            orientationVector_ = x;
            quaternion_ = kine::rotationVectorToAngleAxis(x);
        }

        return quaternion_;
    }

    Vector IMUElasticLocalFrameDynamicalSystem::measureDynamics
                (const Vector& x, const Vector& u, unsigned k)
    {
        assertStateVector_(x);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 accelerationFlex(x.segment(kine::linAcc,3));

        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));
        Vector3 angularAccelerationFlex(x.segment(kine::angAcc,3));

        Quaternion qFlex (computeQuaternion_(orientationFlexV));
        Matrix3 rFlex (qFlex.toRotationMatrix());


        assertInputVector_(u);

//        Vector3 positionControl(u.segment(kine::pos,3));
//        Vector3 velocityControl(u.segment(kine::linVel,3));
//        Vector3 accelerationControl(u.segment(kine::linAcc,3));
//
//        Vector3 orientationControlV(u.segment(kine::ori,3));
//        Vector3 angularVelocityControl(u.segment(kine::angVel,3));

        Vector3 positionControl(u.segment(input::posIMU,3));
        Vector3 velocityControl(u.segment(input::linVelIMU,3));
        Vector3 accelerationControl(u.segment(input::linAccIMU,3));

        Vector3 orientationControlV(u.segment(input::oriIMU,3));
        Vector3 angularVelocityControl(u.segment(input::angVelIMU,3));

        Quaternion qControl(computeQuaternion_(orientationControlV));

        Quaternion q = qFlex * qControl;

        Vector3 acceleration (
         (kine::skewSymmetric(angularAccelerationFlex)
              + tools::square(kine::skewSymmetric(angularVelocityFlex)))
                  * rFlex * positionControl
         + 2*kine::skewSymmetric(angularVelocityFlex) * rFlex * velocityControl
         + accelerationFlex + rFlex * accelerationControl);

        Vector3 angularVelocity( angularVelocityFlex +
                                    rFlex * angularVelocityControl);

        Vector v(Vector::Zero(10,1));

        v[0]=q.w();
        v[1]=q.x();
        v[2]=q.y();
        v[3]=q.z();

        v.segment(4,3)=acceleration;
        v.tail(3)=angularVelocity;

        sensor_.setState(v,k);

        Vector y (Matrix::Zero(measurementSize_,1));

        y.head(sensor_.getMeasurementSize()) = sensor_.getMeasurements();

        for (unsigned i=0; i<contactPositions_.size();++i)
        {
            y.segment(sensor_.getMeasurementSize()+i*3,3)=
                rFlex * contactPositions_[i] + positionFlex - contactPositions_[i];
        }

        return y;
    }

    void IMUElasticLocalFrameDynamicalSystem::setProcessNoise(NoiseBase * n)
    {
        processNoise_=n;
    }

    void IMUElasticLocalFrameDynamicalSystem::resetProcessNoise()
    {
        processNoise_=0x0;
    }

    void IMUElasticLocalFrameDynamicalSystem::setMeasurementNoise
                ( NoiseBase * n)
    {
        sensor_.setNoise(n);
    }

    void IMUElasticLocalFrameDynamicalSystem::resetMeasurementNoise()
    {
        sensor_.resetNoise();
    }

    void IMUElasticLocalFrameDynamicalSystem::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getStateSize()
    {
        return stateSize_;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getInputSize()
    {
        return inputSize_;
    }

    unsigned IMUElasticLocalFrameDynamicalSystem::getMeasurementSize()
    {
        return measurementSize_;
    }

    NoiseBase * IMUElasticLocalFrameDynamicalSystem::getProcessNoise() const
    {
        return processNoise_;
    }

    NoiseBase * IMUElasticLocalFrameDynamicalSystem::
                                                    getMeasurementNoise() const
    {
        return sensor_.getNoise();
    }

    void IMUElasticLocalFrameDynamicalSystem::setContactsNumber(unsigned i)
    {
        measurementSize_ = measurementSizeBase_ + 3 * i;
        contactPositions_.resize(i, Vector3::Zero());
        nbContacts_=i;
    }

    void IMUElasticLocalFrameDynamicalSystem::setContactPosition
                                        (unsigned i, const Vector3 & position)
    {
        BOOST_ASSERT( i< contactPositions_.size() &&
                    "ERROR: The index of contact is out of range.");

        contactPositions_[i] = position;
    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::getContactPosition(unsigned i)
    {
        return contactPositions_[i];
    }


    unsigned IMUElasticLocalFrameDynamicalSystem::getContactsNumber(void)
    {
        return  nbContacts_;
    }
}
}


