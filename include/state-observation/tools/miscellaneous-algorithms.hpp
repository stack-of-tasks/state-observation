/**
 * \file      miscellaneous-algorithms.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */


#ifndef STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
#define STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS

#include <state-observation/tools/definitions.hpp>


namespace stateObservation
{
    namespace tools
    {

        /// Puts the orientation vector norm between 0 and Pi if its
        /// get close to 2pi
        inline Vector regulateOrientationVector(const Vector3 & v )
        {
            double n2=v.squaredNorm();
            if (n2 > (3./2.) * M_PI * (3./2.) * M_PI )
            {
                double n=sqrt(n2);
                unsigned k =  ceil((n - M_PI) / (2*M_PI)) ;
                return (v / n) * ( n - k*2*M_PI );
            }
            else
                return v;
        }

        /// Transform the rotation vector into angle axis
        inline AngleAxis rotationVectorToAngleAxis(const Vector3 & v)
        {
            double angle=v.squaredNorm();
            if (angle > cst::epsilonAngle * cst::epsilonAngle)
            {
                angle=sqrt(angle);
                return AngleAxis(angle, v/angle);
            }
            else
                return AngleAxis(0.0 , Vector3::UnitZ());
        }

        inline Matrix3 skewSymmetric(const Vector3 & v)
        {
            Matrix3 R ;
            R << 0,     -v[2],    v[1],
                 v[2],   0,      -v[0],
                -v[1],   v[0],    0;

            return R;
        }

        template <class T>
        inline T square (const T & x)
        {
            return x*x;
        }

        inline Vector6 homogeneousMatrixToVector6(const Matrix4 & M)
        {
            Vector6 v;
            AngleAxis a = AngleAxis(Matrix3(M.block(0,0,3,3)));

            v.head(3) = M.block(0,3,3,1);
            v.tail(3) = a.angle() * a.axis();

            return v;
        }

        inline Matrix4 vector6ToHomogeneousMatrix(const Vector6 & v)
        {
            Matrix4 M;
            M.block(0,0,3,3) = rotationVectorToAngleAxis(Vector3(v.tail(3)))
                                    .toRotationMatrix();
            M.block(0,3,3,1) = v.head(3);
            M(3,3) = 1.0;

            return M;
        }

        inline void fixedPointRotationToTranslation
            (const Matrix3 & R, const Vector3 & rotationVelocity,
                const Vector3 & rotationAcceleration, const Vector3 & fixedPoint,
                Vector3 & outputTranslation, Vector3 & outputLinearVelocity,
                Vector3 & outputLinearAcceleration)
        {
            Matrix3 omega_x = skewSymmetric(rotationVelocity);
            outputTranslation = fixedPoint - R * fixedPoint;
            outputLinearVelocity = omega_x * R * fixedPoint,
            outputLinearAcceleration =
                    (skewSymmetric(rotationAcceleration) + square(omega_x))
                            * R * fixedPoint;
        }

        template <class T>
        inline T derivate(const T & o1 , const T & o2 , double dt)
        {
            return (o2-o1)/dt;
        }

        inline Vector3 derivateRotationFD
            (const Quaternion & q1, const Quaternion & q2, double dt)
        {
            AngleAxis aa (q2 * q1.inverse());

            double a=aa.angle();

            Vector3 v =  aa.axis();

            return (aa.angle()/dt)*aa.axis();
        }

//        inline Vector3 position(const Vector & v)
//        {
//            return v.head(3);
//        }
//
//        inline Vector3 velocity(const Vector & v)
//        {
//            return v.segment(3,3);
//        }
//
//        inline Vector3 acceleration(const Vector & v)
//        {
//            return v.segment(6,3);
//        }
//
//        inline Vector3 orientation(const Vector & v)
//        {
//            return v.segment(9,3);
//        }
//
//        inline Vector3 angularVelocity(const Vector & v)
//        {
//            return v.segment(12,3);
//        }
//
//        inline Vector3 angularAcceleration(const Vector & v)
//        {
//            return v.tail(3);
//        }

        inline DiscreteTimeArray reconstructStateTrajectory
            (const DiscreteTimeArray & positionOrientation,
             double dt)
        {
            Vector r(Vector::Zero(18,1));

            const DiscreteTimeArray & po= positionOrientation;

            unsigned i0=positionOrientation.getFirstTime();
            unsigned i1=positionOrientation.getLastTime()+1;

            DiscreteTimeArray a;
            a.setValue(r,i0);
            a.resize(po.size(),r);

            for (unsigned i=i0; i<i1; ++i)
            {
                Vector poi = po[i];

                r.head(3)       = poi.head(3);
                r.segment(9,3) = poi.tail(3);
                a.setValue(r,i);
            }

            for (unsigned i=i0; i<i1-1; ++i)
            {
                r = a[i];

                Vector poi = po[i];
                Vector poi1 = po[i+1];

                r.segment(3,3)  = derivate(Vector3(poi.head(3)),
                                            Vector3(poi1.head(3)), dt);
                r.segment(12,3) = derivateRotationFD(
                          Quaternion(rotationVectorToAngleAxis (poi.tail(3))),
                          Quaternion(rotationVectorToAngleAxis (poi1.tail(3))),
                          dt);

                a.setValue(r,i);

            }

            for (unsigned i=i0; i<i1-2; ++i)
            {
                r = a[i];
                Vector r2 = a[i+1];

                r.segment(6,3) = derivate(Vector3(r.segment(6,3)),
                                        Vector3(r2.segment(6,3)), dt);
                r.tail(3) = derivate(Vector3(r.segment(12,3)),
                            Vector3(r2.segment(12,3)), dt);

                a.setValue(r,i);

            }

            return a;
        }
    }
}



#endif //STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
