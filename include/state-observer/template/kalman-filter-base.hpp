/**
 * \file      kalman-filter-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a Kalman filter
 *
 *             It mostly implements the equations of Kalman filtering
 *             It is suitablle by derivation to be used incases of Linear,
 *             linearized and extended Kalman filtering. It may be
 *             derived to unscented Kalman filtering, but non-straighforwardly
 *             because the state vector is modified.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *
 *             y_k=h(x_k,u_k)+w_k
 *
 * \details
 *
 *
 */



#ifndef TEMPLATEKALMANFILTERBASEHPP
#define TEMPLATEKALMANFILTERBASEHPP

#include <state-observer/template/zero-delay-observer.hpp>

namespace stateObserver
{
    namespace compileTime
    {

        /**
             * \class  KalmanFilterBase
             * \brief
             *        It mostly implements the equations of Kalman filtering
             *        It is suitablle by derivation to be used incases of Linear,
             *        linearized and extended Kalman filtering. It may be
             *        derived to unscented Kalman filtering, but non-straighforwardly
             *        because the state vector is modified. This class requires
             *        to be derived to overload the update routine and the measurements
             *        simulation routine.
             *
             *             x_{k+1}=f(x_k,u_k)+v_k
             *
             *             y_k=h(x_k,u_k)+w_k
             *
             *         \li n : size of the state vector
             *         \li m : size of the measurements vector
             *         \li p : size of the input vector
             *
             * \details
             *
             */
        template <unsigned n,unsigned m, unsigned p=0>
        class KalmanFilterBase: public ZeroDelayObserver<n,m,p>
        {
        public:

            /// The type of the jacobian df/dx
            typedef Eigen::Matrix<double, n,n> Amatrix;

            /// The type of the jacobian dh/dx
            typedef Eigen::Matrix<double, m,n> Cmatrix;

            /// The type of the covariance matrix of the process noise v
            typedef Eigen::Matrix<double, n,n> Qmatrix;

            /// The type of the covariance matrix of the measurement noise w
            typedef Eigen::Matrix<double, m,m> Rmatrix;

            /// The type of the covariance matrix of the state estimation error.
            typedef Eigen::Matrix<double, n,n> Pmatrix;



            /// Set the value of the jacobian df/dx
            virtual void setA(const Amatrix& A);

            /// Clear the jacobian df/dx
            virtual void clearA();


            /// Set the value of the Jacobian dh/dx
            virtual void setC(const Cmatrix& C);

            /// Clear the jacobian dh/dx
            virtual void clearC();


            /// Set the measurement noise covariance matrix
            virtual void setR(const Rmatrix& R);

            /// Clear the measurement noise covariance matrix
            virtual void clearR();


            /// Set the process noise covariance matrix
            virtual void setQ(const Qmatrix& Q);

            /// Clear the process noise covariance matrix
            virtual void clearQ();


            /// Set the covariance matrix of the current time state estimation error
            virtual void setStateCovariance(const Pmatrix& P);

            /// Clear the covariace matrix of the current time state estimation
            /// error
            virtual void clearStateCovariance();

            /// Get the covariance matrix of the current time state estimation
            virtual Pmatrix getStateCovariance(unsigned k);

            /// Resets all the observer
            virtual void reset();

        protected:

            /// The type of Kalman gain matrix
            typedef Eigen::Matrix<double, n,m> Kmatrix;

            /// The Kalman filter loop
            virtual typename ObserverBase<n,m,p>::StateVector oneStepEstimation_();

            /// The abstract method to overload to implement f(x,u)
            virtual typename ObserverBase<n,m,p>::StateVector prediction_(unsigned k)=0;

            /// The abstract method to overload to implement h(x,u)
            virtual typename ObserverBase<n,m,p>::MeasureVector simulateSensor_
            (const typename ObserverBase<n,m,p>::StateVector& x, unsigned k)=0;

            /// The internal type for storing the jacobian matrix of the process
            typedef DiscreteTimeMatrix<n,n> Amat_;

            /// The internal type for storing the jacobian matrix of the measurement
            typedef DiscreteTimeMatrix<m,n> Cmat_;

            /// The internal type for storing the process noise covariance matrices
            typedef DiscreteTimeMatrix<n,n> Qmat_;

            /// The internal type for storing the measurement noise covariance matrices
            typedef DiscreteTimeMatrix<m,m> Rmat_;

            /// The internal type for storing the covariance matrix of the
            /// estimation error
            typedef DiscreteTimeMatrix<n,n> Pmat_;

            /// Containers for the jacobian matrix of the process
            Amat_ a_;

            /// Containers for the jacobian matrix of the measurement
            Cmat_ c_;

            /// Container for the process noise covariance matrice
            Qmat_ q_;

            /// Container for the measurement noise covariance matrice
            Rmat_ r_;

            /// Container for the covariance matrix of the estimation error
            Pmat_ p_;

        };


#include <state-observer/template/kalman-filter-base.hxx>
    }
}

#endif //TEMPLATEKALMANFILTERBASEHPP
