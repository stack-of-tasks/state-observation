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
        /// Virtual destructor
        virtual ~NoiseBase(){}

        /// The method to overload to add a noise to a given vector
        virtual Vector addNoise(const Vector &)=0;


    protected:

    };

}

#endif //SENSORSIMULATIONNOISEBASEHPP
