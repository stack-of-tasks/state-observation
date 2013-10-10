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
       /// puts the orientation vector norm between 0 and Pi if its get close to 2pi
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
    }
}



#endif //STATEOBSERVATIONTOOLSMISCELANEOUSALGORITHMS
