#ifndef FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H
#define FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H

#include <state-observation/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp>
#include <vector>

namespace stateObservation
{
    namespace examples
    {
        stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const stateObservation::DiscreteTimeArray & u,
            const Matrix & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3> & contactsPositions,
            double dt);

        stateObservation::DiscreteTimeArray offlineEKFFlexibilityEstimation(
            const stateObservation::DiscreteTimeArray & y,
            const Matrix & xh0,
            unsigned numberOfContacts,
            const std::vector<Vector3> & contactsPositions,
            double dt);



#include <state-observation/examples/offline-ekf-flexibility-estimation.hxx>

    }

}

#endif // FLEXIBILITYESTIMATION_OFFLINEEKFFLEXIBILITYESTIMATION_H
