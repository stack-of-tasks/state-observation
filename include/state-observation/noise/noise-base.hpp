/**
 * \file      noise-base.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 *
 */



#ifndef SENSORSIMULATIONNOISEBASEHPP
#define SENSORSIMULATIONNOISEBASEHPP

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{

    /**
     * \class  NoiseBase
     * \brief
     *
     * \details
     *
     */

    class NoiseBase
    {
    public:
        virtual ~NoiseBase(){}

        virtual Vector addNoise(const Vector & state, const Vector & measurement, unsigned timeIndex)=0;


    protected:

    };

}

#endif //SENSORSIMULATIONNOISEBASEHPP
