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

namespace observation
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
     *         \li n : size of the state vector
     *         \li m : size of the measurements vector
     *         \li p : size of the input vector
     *
     * \details
     *
     */

    template <unsigned n,unsigned m, unsigned p=0>
    class ExtendedKalmanFilter: public KalmanFilterBase<n,m,p>
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
            virtual typename ObserverBase<n,m,p>::StateVector stateDynamics
            (const typename ObserverBase<n,m,p>::StateVector& x,
             const typename ObserverBase<n,m,p>::InputVector& u,
             unsigned k)=0;

            ///The function to overloas to describe the dynamics of the sensor (measurements)
            virtual typename ObserverBase<n,m,p>::MeasureVector measureDynamics
            (const typename ObserverBase<n,m,p>::StateVector& x,
             const typename ObserverBase<n,m,p>::InputVector& u,
             unsigned k)=0;


            ///The method to overload if the functor needs to be reset when the
            ///Exteded Kalman filter is reset itself
            virtual void reset(){}

        };

        /// The constructor. The parameter directInputOutputFeedthrough defines
        /// whether (true) or not (false) the measurement y_k requires the input u_k
        ExtendedKalmanFilter(bool directInputOutputFeedthrough=true):
                directInputOutputFeedthrough_(directInputOutputFeedthrough),f_(0x0)
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
        virtual typename ObserverBase<n,m,p>::StateVector getPrediction(unsigned k);

        ///Give an estimation of A matrix using
        ///finite difference method (the forward difference method)
        ///the parameter dx is the step vector for derivation
        virtual typename KalmanFilterBase<n,m,p>::Amatrix getAMatrixFD(const typename ObserverBase<n,m,p>::StateVector &dx);

        ///Give an estimation of C matrix using
        ///finite difference method (the forward difference method)
        ///the parameter dx is the step vector for derivation
        virtual typename KalmanFilterBase<n,m,p>::Cmatrix getCMatrixFD(const typename ObserverBase<n,m,p>::StateVector &dx);

        /// Reset the extended kalman filter (call also the reset function of the dynamics functor)
        virtual void reset();

    protected:
        /// simulate the dynamics of the state using the functor
        virtual typename ObserverBase<n,m,p>::StateVector prediction_(unsigned k);

        /// simulate the dynamic of the measurement using the functor
        virtual typename ObserverBase<n,m,p>::MeasureVector simulateSensor_
        (const typename ObserverBase<n,m,p>::StateVector& x, unsigned k);

        /// container for the prediction
        typename ObserverBase<n,m,p>::State xbar_;

        /// boolean that provides if theris a need of not for input for the masurement
        bool directInputOutputFeedthrough_;

        /// pointer on the dynamics functor
        DynamicsFunctorBase* f_;
    };



#include <state-observer/extended-kalman-filter.hxx>

}

#endif //EXTENDEDKALMANFILTERHPP
