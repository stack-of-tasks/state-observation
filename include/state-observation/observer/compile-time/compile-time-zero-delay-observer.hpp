/**
 * \file      zero-delay-observer.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of online zero delay observers.
 *             Zero delay observers are the classical state observers where
 *             input and state values at instant k and the measurement value
 *             at instant k+1 are enough to provide the estimation of the state
 *             at instant k+1.
 *
 * \details
 *
 *
 */


#ifndef TEMPLATEZERODELAYOBSERVER_H
#define TEMPLATEZERODELAYOBSERVER_H

#include <deque>

#include <state-observation/observer/compile-time/compile-time-observer-base.hpp>


namespace stateObservation
{
    namespace compileTime
    {
        /**
         * \class  ZeroDelayObserver
         * \brief  Defines the base class of online zero delay observers.
         *         Zero delay observers are the classical state observers where
         *         input and state values at instant k and the measurement value
         *         at instant k+1 are enough to provide the estimation of the state
         *         at instant k+1.
         *         This class mostly defines the data structures for storing the
         *         vectors, it describes the set routines and the observation
         *         loop mechanism. It requires to be derviated to implement the
         *         new oneStepEstimation_() method
         *
         *         \li n : size of the state vector
         *         \li m : size of the measurements vector
         *         \li p : size of the input vector
         *
         * \details
         *
         */
        template <unsigned n,unsigned m, unsigned p=0>
        class ZeroDelayObserver: public ObserverBase<n,m,p>
        {
        public:
            virtual ~ZeroDelayObserver(){};

            ///Set the value of the state vector at time index k. Only the value
            ///with the highest time-index is kept and others are deleted, the
            ///highest index is called the current time k_0
            virtual void setState(
                const typename ObserverBase<n,m,p>::StateVector& x_k,unsigned k);

            ///Remove all the given past values of the state
            virtual void clearState();

            ///Set the value of the measurements vector at time index k. The
            ///measurements have to be inserted in chronological order without gaps.
            virtual void setMeasurement(
                const typename ObserverBase<n,m,p>::MeasureVector& y_k,unsigned k);

            ///Remove all the given past values of the measurements
            virtual void clearMeasurements();

            ///Set the value of the input vector at time index k. The
            ///inputs have to be inserted in chronological order without gaps.
            ///If there is no input in the system (p==0), this instruction has no effect
            virtual void setInput(
                const typename ObserverBase<n,m,p>::InputVector& u_k,unsigned k);

            ///Remove all the given past values of the inputs
            ///If there is no input, this instruction has no effect
            virtual void clearInputs();

            ///Run the observer loop and gets the state estimation of the state at
            ///instant k.
            ///In order to estimate the state k, two conditions hqve to be met:
            /// \li the time index k must be superior to the current time k_0, the
            ///     does *not* record past values of the state and cannot observe
            ///     past states.
            /// \li the observer has to be able to reconstruct all the state
            ///     values from k_0 to k. Thqt means all the measurements or input
            ///     values reauired have to be provided before.
            ///
            ///That means generally (for most zero delay observers) that when
            ///current time is k_0 (we know an estimation of x_{k_0}) and we want
            ///to reconstruct the state at time k>k_0 we need to have the values of
            ///y_{k_0+1} to y_{k} and u_{k_0} to u_{k-1}
            ///
            /// This method sets the current time to k
            virtual typename ObserverBase<n,m,p>::StateVector
            getEstimatedState(unsigned k);

            ///Get the value of the current time index
            virtual unsigned getCurrentTime()const;

        protected:

            ///This method describes one loop of the observer (from k_0 to k_0+1)
            /// it has to be implemented in derived classes.
            virtual typename ObserverBase<n,m,p>::StateVector
            oneStepEstimation_()=0;

            ///while the measurements and iputs are put in lists

            ///The state estimation of the observer (only one state is recorded)
            typename ObserverBase<n,m,p>::State x_;

            ///Container for the measurements.
            std::deque<typename ObserverBase<n,m,p>::Measure,Eigen::aligned_allocator<typename ObserverBase<n,m,p>::Measure> > y_;

            ///Container for the actual measurements.
            std::deque<typename ObserverBase<n,m,p>::Input,Eigen::aligned_allocator<typename ObserverBase<n,m,p>::Input>  > u_;

        };
#include <state-observation/observer/compile-time/compile-time-zero-delay-observer.hxx>
    }
}



#endif //TEMPLATEZERODELAYOBSERVER
