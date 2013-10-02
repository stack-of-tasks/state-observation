/**
 * \file      sensor-base.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 */



#ifndef SIMULATIONSENSORBASEHPP
#define SIMULATIONSENSORBASEHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/noise/noise-base.hpp>

namespace stateObservation
{
    /**
     * \class  ObserverBase
     * \brief  The base class for observers.
     *
     * \details
     *
     */

    class SensorBase
    {
    public:
        SensorBase();

        virtual ~SensorBase(){}

        virtual Vector getMeasurements(bool noisy=true)=0;

        virtual void setState(const Vector & state, unsigned k)=0;

        virtual void setNoise(NoiseBase *);

        virtual void resetNoise();

        virtual unsigned getTime() const=0;

        virtual unsigned getStateSize() const=0;

        virtual unsigned getMeasurementSize() const=0;

        virtual Vector stateVectorZero() const;

        virtual bool checkStateVector(const Vector &) const;


    protected:

        NoiseBase *  noise_;

    };



}

#endif //SIMULATIONSENSORBASEHPP
