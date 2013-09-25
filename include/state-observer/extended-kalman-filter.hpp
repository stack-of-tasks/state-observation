/**
 * \file      extended-kalman-filter-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defined the class to intanciate to use an extended Kalman filter.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *             y_k=h(x_k,u_k)+w_k
 *
 * \details
 *
 *
 */


#ifndef STATEOBSERVER_EXTENDEDKALMANFILTERHPP
#define STATEOBSERVER_EXTENDEDKALMANFILTERHPP

#include <state-observer/kalman-filter-base.hpp>

namespace stateObserver
{
/**
     * \class  ExtendedKalmanFilter
     * \brief
     *
     *        The class to intanciate to use an extended Kalman filter.
     *        To use this class, one needs to provide a pointer on a functor
     *        that describes the state dynamics and the measurement dynamics.
     *        The functor type needs to be a derived class from the class
     *        DynamicsFunctorBase.
     *
     *        x_{k+1}=f(x_k,u_k)+v_k
     *
     *        y_k=h(x_k,u_k)+w_k
     *
     *
     *
     */

    class ExtendedKalmanFilter: public KalmanFilterBase
    {

        /**
     * \class  DynamicsFunctorBase
     * \brief
     *        This is the base class of any functor that describes the dynamics
     *        of the state and the measurement.
     *        This class is to be derived in order to be given
     *        to the Extended Kalman Filter.
     *
     */
    public:
        class DynamicsFunctorBase
        {
        public:
            ///The function to oberload to describe the dynamics of the state
            virtual ObserverBase::StateVector stateDynamics
            (const ObserverBase::StateVector& x,
             const ObserverBase::InputVector& u,
             unsigned k)=0;

            ///The function to overloas to describe the dynamics of the sensor (measurements)
            virtual ObserverBase::MeasureVector measureDynamics
            (const ObserverBase::StateVector& x,
             const ObserverBase::InputVector& u,
             unsigned k)=0;


            ///The method to overload if the functor needs to be reset when the
            ///Exteded Kalman filter is reset itself
            virtual void reset(){}

        };

        /// The constructor.
        ///  \li n : size of the state vector
        ///  \li m : size of the measurements vector
        ///  \li p : size of the input vector
        ///  \li The parameter directInputOutputFeedthrough defines whether (true) or not (false) the measurement y_k requires the input u_k
        ExtendedKalmanFilter(unsigned n,unsigned m,unsigned p=0,bool directInputOutputFeedthrough=true)
            :KalmanFilterBase(n,m,p)
        {
            if (p==0)
                directInputOutputFeedthrough=false;
        }


        /// Set a pointer to the functor that defines the dynamics of the states
        ///and the measurement the user is responsible for the validity of the
        ///pointer during the execution of the kalman filter
        void setFunctor(DynamicsFunctorBase* f);

        /// Clear the value of the functor
        ///Does not destroy the pointed object
        void clearFunctor();

        /// Precise whether (true) or not (false) the measurement y_k requires
        ///the input u_k
        void setDirectInputOutputFeedthrough(bool b=true);

        /// A function that gives the prediction (this is NOT the estimation of the state,
        /// for the estimation call getEstimateState method
        /// it is only an execution of the state synamics with the current state
        /// estimation and the current input value
        virtual StateVector getPrediction(unsigned k);

        ///Give an estimation of A matrix using
        ///finite difference method (the forward difference method)
        ///the parameter dx is the step vector for derivation
        virtual Amatrix getAMatrixFD(const StateVector &dx);

        ///Give an estimation of C matrix using
        ///finite difference method (the forward difference method)
        ///the parameter dx is the step vector for derivation
        virtual Cmatrix getCMatrixFD(const StateVector &dx);

        /// Reset the extended kalman filter (call also the reset function of the dynamics functor)
        virtual void reset();

    protected:
        /// simulate the dynamics of the state using the functor
        virtual StateVector prediction_(unsigned k);

        /// simulate the dynamic of the measurement using the functor
        virtual MeasureVector simulateSensor_(const StateVector& x, unsigned k);

        /// container for the prediction
        DiscreteTimeMatrix xbar_;

        /// boolean that provides if theris a need of not for input for the masurement
        bool directInputOutputFeedthrough_;

        /// pointer on the dynamics functor
        DynamicsFunctorBase* f_;
    };

}

#endif //EXTENDEDKALMANFILTERHPP
