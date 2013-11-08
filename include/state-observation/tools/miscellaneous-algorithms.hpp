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

            if (v.squaredNorm() > (3./2.) * M_PI * (3./2.) * M_PI )
            {
                double n=v.norm();
                return v*( n - 2*M_PI )/n;
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

        inline Vector3 derivateRotationFD
            (const Quaternion & q1, const Quaternion & q2, double dt)
        {
            AngleAxis aa (q2 * q2.inverse());

            return (aa.angle()/dt)*aa.axis();
        }

    }
}



#endif //STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
