#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>

namespace stateObservation
{
    unsigned AccelerometerGyrometer::getStateSize() const
    {
        return stateSize_;
    }

    unsigned AccelerometerGyrometer::getMeasurementSize() const
    {
        return measurementSize_;
    }

    Vector AccelerometerGyrometer::computeNoiselessMeasurement_()
    {
        Quaternion q(state_[0],state_[1],state_[2],state_[3]);
        Matrix3 r(q.toRotationMatrix());

        Vector3 acceleration = state_.segment(4,3);
        Vector3 rotationVector = state_.tail(3);

        Vector output=Vector::Zero(measurementSize_,1);

        output.head(3)=accelerationMeasure(acceleration,r);
        output.tail(3)=rotationVelocityMeasure(rotationVector, r);

        return output;
    }

}
