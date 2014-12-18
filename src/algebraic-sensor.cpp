#include <state-observation/sensors-simulation/algebraic-sensor.hpp>

namespace stateObservation
{
    AlgebraicSensor::AlgebraicSensor()
        :time_(0),storedNoisyMeasurement_(),storedNoiselessMeasurement_(false)
    {
    }

    void AlgebraicSensor::setState(const Vector & state, unsigned k)
    {
        if (state.size()!=state_.size() || state!=state_ || time_!=k)
        {
            state_=state;
            time_=k;
            storedNoisyMeasurement_=false;
            storedNoiselessMeasurement_=false;
        }
    }

    Vector AlgebraicSensor::getMeasurements(bool noisy)
    {


        if (noisy && noise_!=0x0)
        {
            if (!storedNoisyMeasurement_)
            {
                noisyMeasurement_ =  computeNoisyMeasurement_();
                storedNoisyMeasurement_=true;
            }

            return noisyMeasurement_;
        }
        else
        {
            if (!storedNoiselessMeasurement_)
            {
                noiselessMeasurement_ = computeNoiselessMeasurement_();
                storedNoiselessMeasurement_=true;
            }
            return noiselessMeasurement_;
        }

    }

    void AlgebraicSensor::checkState_(const Vector & v)
    {
        (void)v;//avoid warning
        BOOST_ASSERT( checkStateVector(v) && "ERROR: The state vector is incorrectly set.");
    }

    unsigned AlgebraicSensor::getTime() const
    {
        return time_;
    }

    Vector AlgebraicSensor::computeNoisyMeasurement_()
    {
        if (!storedNoiselessMeasurement_)
        {
            BOOST_ASSERT(checkStateVector(state_)&& "The state is not set or incorrectly set");

            noiselessMeasurement_ = computeNoiselessMeasurement_();
            storedNoiselessMeasurement_=true;
        }

        return noise_->addNoise( noiselessMeasurement_);
    }


}
