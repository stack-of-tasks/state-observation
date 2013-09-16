/**
 * \file      observer-base.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the base class of a Kalman filter
 *
 *             It mostly implements the equations of Kalman filtering
 *             It is suitablle by derivation to be used incases of Linear,
 *             linearized and extended Kalman filtering. It can may be
 *             derived to unscented Kalman filtering, but non-straighforwardly
 *             because the state vector is modified.
 *
 *             x_{k+1}=f(x_k,u_k)+v_k
 *             y_k=h(x_k,u_k)
 *
 * \details
 *
 *
 */



#ifndef KALMANFILTERBASEHPP
#define KALMANFILTERBASEHPP

#include <state-observer/zero-delay-observer.hpp>

namespace observation
{


    /**
     * \class  KalmanFilterBase
     * \brief  The base class of linear /extended Kalman filters.
     *         It implements one loop of the Kalman filter. This class requires
     *         to be derived to overload the update routine and the measurements
     *         simulation routine.
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



        /// The type of the jacobian df/dx to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> AmatrixMember;

        /// The type of the jacobian dh/dx to be used to set as class members
        typedef Eigen::Matrix<double, m,n,Eigen::DontAlign> CmatrixMember;

        /// The type of the covariance matrix of the process noise v to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> QmatrixMember;

        /// The type of the covariance matrix of the measurement noise w to be used to set as class members
        typedef Eigen::Matrix<double, m,m,Eigen::DontAlign> RmatrixMember;

        /// The type of the covariance matrix of the state estimation error. to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> PmatrixMember;



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

        /// The internal type for storing the jacobian matrices
        typedef DiscreteTimeMatrix<n,n> Amat_;
        typedef DiscreteTimeMatrix<m,n> Cmat_;

        /// The internal type for storing the noises covariance matrices
        typedef DiscreteTimeMatrix<n,n> Qmat_;
        typedef DiscreteTimeMatrix<m,m> Rmat_;

        /// The internal type for storing the covariance matrix of the
        /// estimation error
        typedef DiscreteTimeMatrix<n,n> Pmat_;

        /// Containers for the jacobian matrices
        Amat_ a_;
        Cmat_ c_;

        /// Containers for the noises covariance matrices
        Qmat_ q_;
        Rmat_ r_;

        /// Container for the covariance matrix of the estimation error
        Pmat_ p_;

    };

    /**
     * \class  KalmanFilterBase
     * \brief  The base class of linear /extended Kalman filters, specialized for
     *         the case p=0 (there is no input for the system)
     *         There is no need to fill the input vector in that case
     *
     *         \li n : size of the state vector
     *         \li m : size of the measurements vector
     *
     * \details
     *
     */
    template <unsigned n,unsigned m>
    class KalmanFilterBase<n,m,0>: public ZeroDelayObserver<n,m,0>
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




        /// The type of the jacobian df/dx to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> AmatrixMember;

        /// The type of the jacobian dh/dx to be used to set as class members
        typedef Eigen::Matrix<double, m,n,Eigen::DontAlign> CmatrixMember;

        /// The type of the covariance matrix of the process noise v to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> QmatrixMember;

        /// The type of the covariance matrix of the measurement noise w to be used to set as class members
        typedef Eigen::Matrix<double, m,m,Eigen::DontAlign> RmatrixMember;


        /// The type of the covariance matrix of the state estimation error. to be used to set as class members
        typedef Eigen::Matrix<double, n,n,Eigen::DontAlign> PmatrixMember;




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
        virtual typename ObserverBase<n,m,0>::StateVector oneStepEstimation_();

        /// The abstract method to overload to implement f(x,u)
        virtual typename ObserverBase<n,m,0>::StateVector prediction_(unsigned k)=0;

        /// The abstract method to overload to implement h(x,u)
        virtual typename ObserverBase<n,m,0>::MeasureVector
        simulateSensor_(const typename ObserverBase<n,m,0>::StateVector& x,unsigned k)=0;


        /// Containers for the jacobian matrices
        typedef DiscreteTimeMatrix<n,n>  Amat_;
        typedef DiscreteTimeMatrix<m,n>  Cmat_;

        /// The internal type for storing the noises covariance matrices
        typedef DiscreteTimeMatrix<n,n> Qmat_;
        typedef DiscreteTimeMatrix<m,m> Rmat_;

        /// The internal type for storing the covariance matrix of the
        /// estimation error
        typedef DiscreteTimeMatrix<n,n> Pmat_;


        /// Containers for the jacobian matrices
        Amat_ a_;
        Cmat_ c_;

        /// Containers for the noises covariance matrices
        Qmat_ q_;
        Rmat_ r_;

        /// Container for the covariance matrix of the estimation error
        Pmat_ p_;

    };


#include <state-observer/kalman-filter-base.hxx>

}

#endif //KALMANFILTERBASEHPP
