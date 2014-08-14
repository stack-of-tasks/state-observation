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
       // std::cout<<std::endl<<"IMUElasticLocalFrameDynamicalSystem Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR

       calculationState <<  -1,
                           -1,
                           -1;

    }

    IMUElasticLocalFrameDynamicalSystem::
                    ~IMUElasticLocalFrameDynamicalSystem()
    {
        //dtor
    }


    Vector3 IMUElasticLocalFrameDynamicalSystem::computeFc(const stateObservation::Vector& x, const stateObservation::Vector& u)
    {

        unsigned nbContacts(getContactsNumber());
        Fci.resize(3,nbContacts);
        Fci.setZero();
        //std::cout << nbContacts << std::endl;
        unsigned i;

        Vector3 Fc, Fcu;
        Fc << 0,0,0;

        // Isotropic stifness and viscosity for a simple contact case
        double kfe=hrp2::linKe;
        double kfv=hrp2::linKv;

        Matrix3 Kte, Ktv, Kfe, Kfv;

        Matrix4 homoi;
        Matrix3 Rci;
        Vector3 tci;

        // Flexibility state
        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlexV(x.segment(kine::angVel,3));

        Vector6 posFlex;
        Matrix3 Rflex;
        Vector3 tflex;

        posFlex << positionFlex,orientationFlexV;
        Matrix4 homoFlex(kine::vector6ToHomogeneousMatrix(posFlex));
        Rflex = homoFlex.block(0,0,3,3);
        tflex = homoFlex.block(0,3,3,1);


        Kfe <<  kfe,0,0,
                0,kfe,0,
                0,0,kfe;

        Kfv <<  kfv,0,0,
                0,kfv,0,
                0,0,kfv;

        for(i=0;i<nbContacts;i++)
        {
            homoi = kine::vector6ToHomogeneousMatrix(u.segment(input::contacts+6*i,6));
                                            //std::cout << u.segment(input::contacts+6*i,6) << std::endl;
            Rci = homoi.block(0,0,3,3);
            tci = homoi.block(0,3,3,1);

            Fcu.noalias() = - Rci*Kfe*Rci.transpose()*(Rflex*tci+tflex-tci);
//            std::cout << "Fc" << i << "0 " << Fc.transpose() << std::endl;
//            if(sqrt(Fc.squaredNorm())>100)
//            {
//                std::cout << "Rci*Kfe*Rci.transpose()" << Rci*Kfe*Rci.transpose() << std::endl;
//                std::cout << "Rflex*tci+tflex-tci" << Rflex*tci+tflex-tci << std::endl;
//            }
            Fcu.noalias() += - Rci*Kfv*Rci.transpose()*(kine::skewSymmetric(angularVelocityFlexV)*Rflex*tci+velocityFlex);
//            std::cout << "Fc" << i << "1 " << Fc.transpose() << std::endl;
//            if(sqrt(Fc.squaredNorm())>100)
//            {
//                std::cout << "Rci*Kfv*Rci.transpose()" << Rci*Kfv*Rci.transpose() << std::endl;
//                std::cout << "kine::skewSymmetric(angularVelocityFlexV)*Rflex*tci+velocityFlex" << kine::skewSymmetric(angularVelocityFlexV)*Rflex*tci+velocityFlex << std::endl;
//            }

//            std::cout << "Fc" << i << " " << Fc.transpose() << std::endl;
//            if(sqrt(Fc.squaredNorm())>100)
//            {
//                std::cout << "Rci" << Rci << std::endl;
//                std::cout << "tci" << tci.transpose() << std::endl;
//                std::cout << "Rflex" << Rflex << std::endl;
//                std::cout << "tflex" << tflex.transpose() << std::endl;
//                std::cout << "angularVelocityFlexV" << angularVelocityFlexV.transpose() << std::endl;
//                std::cout << "velocityFlex" << velocityFlex.transpose() << std::endl;
//            }


            Fci.block(0,i,3,1)=Fcu;
            Fc += Fcu;

//            if(i>0)
//            {
//                Fci.block(0,i,3,1)=Fc-Fci.block(0,i-1,3,1);
//            }

        }

        return Fc;
    }

    double IMUElasticLocalFrameDynamicalSystem::rotationMatrixFromContactsPositiont(const Vector3 contact1, const Vector3 contact2, Matrix3& R )
    {

        Vector3 Vrl, axis, theta;
        AngleAxis j;
        double h, thetay, thetaz;

        // Definition of the length and the unit vector between the two feet.
        Vrl=contact1-contact2;
        if(Vrl(1)<0) Vrl=-Vrl ; // make sure the axis will be in the right side

        j=kine::rotationVectorToAngleAxis(Vrl);
        h = j.angle(); // length between the ankles
        axis = j.axis(); // unit vector between the ankles

        // Taking care for the first iteration were there is no still contacts position
        if(axis(2) == 1 && axis(1) == 0 && axis(0) == 0)
        {
           axis <<  0,
                    1,
                    0;
           h=0.19;
        }

        // Definition of the transformation (rotation) between (x,y,z) and (perpendicular of j, j, z).
        /// In this case there is onlys z rotation
        theta=kine::unitVectorToRotationVector(axis);
        thetaz = theta[2];
        R <<    cos(thetaz),-sin(thetaz),0,
                sin(thetaz),cos(thetaz),0,
                0,0,1;

        return h;
    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::computeTc(const stateObservation::Vector& x, const stateObservation::Vector& u)
    {
        unsigned nbContacts(getContactsNumber());
        unsigned i;

        Vector3 Tc;
        //Vector3 Fc(computeFc(x,u));
        Tc << 0,0,0;

        // Isotropic stifness and viscosity for a simple contact case
        double kte=hrp2::angKe;
        double ktv=hrp2::angKv;

        Matrix3 Kte, Ktv, Kfe, Kfv;

        Matrix4 homoi;
        Matrix3 Rci;
        Vector3 tci;

        // Flexibility state
        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlexV(x.segment(kine::angVel,3));

        Vector6 posFlex;
        Matrix3 Rflex;
        Vector3 tflex;

        posFlex << positionFlex,orientationFlexV;
        Matrix4 homoFlex(kine::vector6ToHomogeneousMatrix(posFlex));
        Rflex = homoFlex.block(0,0,3,3);
        tflex = homoFlex.block(0,3,3,1);

        Kte <<  kte,0,0,
                0,kte,0,
                0,0,kte;

        Ktv <<  ktv,0,0,
                0,ktv,0,
                0,0,ktv;

        for(i=0;i<nbContacts;i++)
        {
            homoi = kine::vector6ToHomogeneousMatrix(u.segment(input::contacts+6*i,6));
            Rci = homoi.block(0,0,3,3);
            tci = homoi.block(0,3,3,1);

//            std::cout << "conttact" << i+1 << "\t Rci: " << Rci << "\n tci: " << tci.transpose() << std::endl;

            Tc.noalias() += -Rci*Kte*Rci.transpose()*orientationFlexV;
            Tc.noalias() += -Rci*Ktv*Rci.transpose()*angularVelocityFlexV;
            Tc.noalias() += kine::skewSymmetric(Rflex*tci+tflex)*Fci.block(0,i,3,1);

//            std::cout << "Tc" << Tc << std::endl;

        }

        return Tc;

    }


    void IMUElasticLocalFrameDynamicalSystem::test()
    {
        Vector x, u;
        x.resize(18);
        u.resize(54);
//        x << 1.95675e-08,-2.22298e-10,-0.0109575,1.54011e-12,1.35566e-10,3.81192e-23,1.41111e-06,-1.60311e-08,-0.0489989,6.56636e-12,-3.14893e-09,-9.62102e-11,-5.76107e-05,6.54491e-07,9.47828,-8.92673e-08,-8.41481e-06,-7.82905e-09;
//        u << 0.0135673,0.001536,0.80771,-2.63605e-06,-1.09258e-08,5.71759e-08,2.71345,0.3072,161.542,1.55342e-318,-6.53681e-66,-6.39612e-66,-6.54978e-66,-6.53171e-66,-6.55957e-66,0,0,0,0,0,0,0,0,0,0,0,0,-0.098,-6.23712e-11,1.1174,1.56933e-22,-5.40778e-21,3.86235e-22,-2.99589e-06,-1.24742e-08,-4.76329e-18,3.13865e-20,-1.08156e-18,7.72471e-20,-0.000299589,-1.24742e-06,-4.76329e-16,0.00949046,-0.095,1.98197e-07,8.75907e-23,-0,7.17001e-22,0.00949046,0.095,1.98197e-07,6.29302e-23,1.01427e-16,-6.61701e-22;

        x << 0,0,-0.0109575,0,0,0,0,0,-0.0489989,0,0,0,0,0,0,0,0,0;
        u << 0.0135673,0.001536,0.80771,0,0,0,2.71345,0.3072,161.542,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.098,0,1.1174,0,0,0,0,0,0,0,0,0,-0.000299589,0,0,0.00949046,-0.095,0,0,-0,0,0.00949046,0.095,0,0,0,0;

//        std::cout << "AccAng" << computeAccelerationAngular(x,u,1) << std::endl;
//        std::cout << "AccLin" << computeAccelerationLinear(x,u,1) << std::endl;
    }


    Vector3 IMUElasticLocalFrameDynamicalSystem::computeAccelerationAngular
    	(const Vector& x, const Vector& u, unsigned k)
    {

        int i,j;
        bool q;

        // Vector we want
        Vector3 AccAngular;

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

        // To be human readable
        const Vector3 Fc(computeFc(x,u)); // first because Tc neef Fci where i is the indice of the contact.
        const Vector3 Tc(computeTc(x,u));
        const Quaternion qFlex (computeQuaternion_(orientationFlexV));
        const Matrix3 R (qFlex.toRotationMatrix());
        const Matrix3 Inertia(kine::computeInertiaTensor(u.segment(input::Inertia,6)));
        const Matrix3 dotInertia(kine::computeInertiaTensor(u.segment(input::dotInertia,6)));

       // std::cout << "Input again" << u.transpose() << std::endl;

        Mat = R*Inertia*R.transpose();
     //   std::cout << "Mat 0 \n"  << Mat << std::endl;
        Mat += hrp2::m*kine::skewSymmetric(R*positionCom)*kine::skewSymmetric(R*positionCom);
     //   std::cout << "R\n" << R << "\npositionCom\n" << positionCom << std::endl;
      //  std::cout << "kine::skewSymmetric(R*positionCom)\n" << kine::skewSymmetric(R*positionCom) << std::endl;
      //  std::cout << "kine::skewSymmetric(R*positionCom)*kine::skewSymmetric(R*positionCom)\n" << kine::skewSymmetric(R*positionCom)*kine::skewSymmetric(R*positionCom) << std::endl;
      //  std::cout << "Mat 1 \n"  << Mat << std::endl;

        q=false;
        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                if(Mat(i,j)<=0)
                {
                    q=false;//true;
                }
            }
        }

        if(q==true)
        {
            Mat <<  0,0,0,
                    0,0,0,
                    0,0,0;
        }
        else
        {
            Mat = Mat.inverse().eval();
        }


      //  std::cout << "Mat 2 \n"  << Mat << std::endl;

//        if(Mat!=Mat)// || isinf(sqrt(Mat.squaredNorm())))
//        {
//            Mat <<  0,0,0,
//                    0,0,0,
//                    0,0,0;
//        }
//
//        std::cout << "Mat 3 "  << Mat << std::endl;

        Vec = -(    (kine::skewSymmetric(angularVelocityFlex)*R*Inertia*R.transpose()+R*dotInertia*R.transpose())*angularVelocityFlex
                                +R*dotAngMomentum
                                +kine::skewSymmetric(angularVelocityFlex)*R*AngMomentum
                                +hrp2::m*kine::skewSymmetric(positionFlex)*
                                    (
                                        kine::skewSymmetric(angularVelocityFlex)*kine::skewSymmetric(angularVelocityFlex)*R*positionCom
                                        +2*kine::skewSymmetric(angularVelocityFlex)*R*velocityCom
                                        +R*accelerationCom
                                    )
                                +hrp2::m*kine::skewSymmetric(R*positionCom+positionFlex)*cst::gravity
                );
//        std::cout << "Vec 0 "  << Vec.transpose() << std::endl;
        Vec -= (kine::skewSymmetric(positionFlex)+kine::skewSymmetric(R*positionCom))*
                    (   Fc
                        -(  R*hrp2::m*accelerationCom
                            +2*kine::skewSymmetric(orientationFlexV)*R*velocityCom
                            +hrp2::m*kine::skewSymmetric(orientationFlexV)*kine::skewSymmetric(orientationFlexV)*R*positionCom
                            +hrp2::m*cst::gravity
                         )
                    );
//        std::cout << "Vec 1 "  << Vec.transpose() << std::endl;
        Vec += Tc;
//        std::cout << "Vec 2 "  << Vec.transpose() << std::endl;

        AccAngular = Mat*Vec;
//        std::cout << "Acc angular" << Mat*Vec << std::endl;
        return AccAngular;

    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::computeAccelerationLinear
    	(const Vector& x, const Vector& u, unsigned k)
    {
        // Vector we want
        Vector3 AccLinear;

        // State vector
        Vector3 orientationFlexV(x.segment(kine::ori,3));           // \Omega (position)

        // To be human readable
        const Vector3 Fc(computeFc(x,u));
        //std::cout << "Fc: " << Fc.transpose() << std::endl;
        const Quaternion qFlex (computeQuaternion_(orientationFlexV));
        const Matrix3 R (qFlex.toRotationMatrix());

        // Input vector
        Vector3 positionCom(u.segment(input::posCom,3));
        Vector3 velocityCom(u.segment(input::velCom,3));
        Vector3 accelerationCom(u.segment(input::accCom,3));

       //         std::cout << "x " << x.transpose() << std::endl;
       // std::cout << "u " << u.transpose() << std::endl;

        AccLinear = Fc;
//        std::cout << "accLinear Fc " << AccLinear.transpose() << std::endl;
        AccLinear -=  R*hrp2::m*accelerationCom
                                    +2*kine::skewSymmetric(orientationFlexV)*R*velocityCom
                                    +hrp2::m*kine::skewSymmetric(orientationFlexV)*kine::skewSymmetric(orientationFlexV)*R*positionCom
                                    +hrp2::m*cst::gravity;
//               std::cout << "accLinear 1 " << AccLinear.transpose() << std::endl;
        AccLinear /= hrp2::m;
//                std::cout << "accLinear 2 " << AccLinear.transpose() << std::endl;
        AccLinear += kine::skewSymmetric(R*positionCom)*computeAccelerationAngular(x,u,k);


//        std::cout << "accLinear fin" << AccLinear.transpose() << std::endl;

        return AccLinear;

    }

    Vector IMUElasticLocalFrameDynamicalSystem::stateDynamics
        (const Vector& x, const Vector& u, unsigned k)
    {

        assertStateVector_(x);

        Vector3 positionFlex(x.segment(kine::pos,3));
        Vector3 velocityFlex(x.segment(kine::linVel,3));
        Vector3 accelerationFlex(x.segment(kine::linAcc,3));
        Vector3 orientationFlexV(x.segment(kine::ori,3));
        Vector3 angularVelocityFlex(x.segment(kine::angVel,3));
        Vector3 angularAccelerationFlex(x.segment(kine::angAcc,3));


        Quaternion orientationFlex(computeQuaternion_(orientationFlexV));

        integrateKinematics(positionFlex, velocityFlex, accelerationFlex, orientationFlex,angularVelocityFlex, angularAccelerationFlex, dt_);

        accelerationFlex = computeAccelerationLinear(x, u, k);
        angularAccelerationFlex = computeAccelerationAngular(x, u, k);
//        std::cout << "accelerationFlex" << accelerationFlex.transpose() << std::endl;
//        std::cout << "angularAccelerationFlec" << angularAccelerationFlex.transpose() << std::endl;

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

//       std::cout << "===> x  " << x.transpose() << std::endl;
//       std::cout << "===> u  " << u.transpose() << std::endl;
//       std::cout << "===> xk1  " << xk1.transpose() << std::endl;
        //std::cout << "accelerationFlex  " << accelerationFlex.transpose() << std::endl;
        //std::cout << "angularAccelerationFlex  " << angularAccelerationFlex.transpose() << std::endl;
       // std::cout << "xhk1 " << xk1.transpose() << std::endl;

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

        Vector3 positionControl(u.segment(input::posIMU,3));
        Vector3 velocityControl(u.segment(input::linVelIMU,3));
        Vector3 accelerationControl(u.segment(input::linAccIMU,3));
        Vector3 orientationControlV(u.segment(input::oriIMU,3));
        Vector3 angularVelocityControl(u.segment(input::angVelIMU,3));

        Quaternion qControl(computeQuaternion_(orientationControlV));
        Quaternion q = qFlex * qControl;

//        std::cout << "Mesure x " << x.transpose() << std::endl;
//        std::cout << "Mesure u " << u.transpose() << std::endl;

        // Translation sensor dynamic
        Vector3 acceleration;

        acceleration << 0,0,0;

        acceleration += 2*kine::skewSymmetric(angularVelocityFlex) * rFlex * velocityControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        acceleration += accelerationFlex + rFlex * accelerationControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        acceleration += (
                    kine::skewSymmetric(angularAccelerationFlex)
                    + tools::square(kine::skewSymmetric(angularVelocityFlex))
                )*rFlex * positionControl;

//        std::cout << "acceleration " << acceleration.transpose() << std::endl;

        // Rotation sensor dynamic
        Vector3 angularVelocity( angularVelocityFlex + rFlex * angularVelocityControl);

        // Set sensor state before measurement
        Vector v(Vector::Zero(10,1));
        v[0]=q.w();
        v[1]=q.x();
        v[2]=q.y();
        v[3]=q.z();
        v.segment(4,3)=acceleration;
        v.tail(3)=angularVelocity;
        sensor_.setState(v,k);

        // Measurement
        Vector y (Matrix::Zero(measurementSize_,1));
        y.head(sensor_.getMeasurementSize()) = sensor_.getMeasurements();

//        std::cout << "===> v state sensor  " << v.transpose() << std::endl;
//        std::cout << "===> y  " << y.transpose() << std::endl;

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

    void IMUElasticLocalFrameDynamicalSystem::setInputSize(unsigned i)
    {
        inputSize_=i;
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

    Vector3 IMUElasticLocalFrameDynamicalSystem::getFc(unsigned k)
    {
        Vector3 Fc;
        unsigned nbContacts(getContactsNumber());
        unsigned i;
        Fc  <<  0,
                0,
                0;

        for(i=0;i<nbContacts;++i)
        {
            Fc += getFc(i, k);
        }

        return Fc;
    }

    Vector3 IMUElasticLocalFrameDynamicalSystem::getFc(unsigned i, unsigned k)
    {
        if(calculationState(i+1)==k)
        {
            return Fci.block(0,i,3,1);
        }
        else
        {
            // on calcul tous les Fci, on met à jour calculationState et on retourne le Fci voulu;
            // compute Fci
            calculationState(i+1)=k;
            return Fci.block(0,i,3,1);
        }
    }
}
}


