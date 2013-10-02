/**
 * \file      algebraic-sensor.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 */



#ifndef SIMULATIONALGEBRAICSENSORHPP
#define SIMULATIONALGEBRAICSENSORHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/sensors-simulation/sensor-base.hpp>

namespace stateObservation
{
    /**
     * \class  ObserverBase
     * \brief  The base class for observers.
     *
     * \details
     *
     */

    class AlgebraicSensor: public SensorBase
    {
    public:
        AlgebraicSensor();

        virtual ~AlgebraicSensor(){}

        virtual Vector getMeasurements(bool noisy=true);

        virtual void setState(const Vector & state, unsigned k);

        virtual unsigned getTime() const;

        virtual unsigned getStateSize() const=0;

        virtual unsigned getMeasurementSize() const=0;

    protected:
        virtual Vector computeNoiselessMeasurement_()=0;

        virtual Vector computeNoisyMeasurement_();

        virtual void checkState_(const Vector &);

        unsigned time_;

        Vector state_;

        bool storedNoisyMeasurement_;

        Vector noisyMeasurement_;

        bool storedNoiselessMeasurement_;

        Vector noiselessMeasurement_;

    };

}

#endif //SIMULATIONALGEBRAICSENSORHPP
