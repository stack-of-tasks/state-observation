#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/random.hpp>

#include <state-observer/kalman-filter.hpp>
#include <state-observer/extended-kalman-filter.hpp>


typedef observation::KalmanFilter<4,3,2> filter;
typedef observation::ExtendedKalmanFilter<4,3,1> ekf;
filter f;
ekf f2;

class KalmanFunctor:
            public ekf::DynamicsFunctorBase
{
public:
    KalmanFunctor()
    {
        s_=ekf::StateVector::Random()*0.1;
        m_=ekf::MeasureVector::Random()*0.1;
        t_=ekf::StateVector::Random();
        n_=ekf::MeasureVector::Random();
    }

    virtual ekf::StateVector stateDynamics(const ekf::StateVector& x, const ekf::InputVector& u, unsigned k)
    {
        (void)k;//unused
        (void)u;//unused

        std::cout<<"here";
        ekf::StateVector xk1;
        xk1=a_*x+(x.transpose()*x)[0]*s_+t_+(u*u.transpose())[0]*s_;
        return xk1;
    }

    virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
    {
        (void)k;//unused
        (void)u;//unused


        ekf::MeasureVector yk;
        yk=c_*x+(x.transpose()*x)[0]*m_+n_;
        return yk;
    }

    void setA(const ekf::Amatrix& a)
    {
        a_=a;
    }

    void setC(const ekf::Cmatrix& c)
    {
        c_=c;
    }

    void doNothing () {}

private:
    ekf::StateVectorMember s_;
    ekf::StateVectorMember t_;
    ekf::MeasureVectorMember m_;
    ekf::MeasureVectorMember n_;

    ekf::AmatrixMember a_;
    ekf::CmatrixMember c_;

};


void testExtendedKalmanFilter()
{
    ekf::Amatrix a;
    KalmanFunctor func;
    KalmanFunctor* debug=&func;
    f2.setFunctor(&func);
    a<<	-0.6785714 ,0.1156463 ,	0.4392517,	0.2863946  ,
    0.0865306 ,	-0.0273810,	0.3355102  ,0.0184150 ,
    - 0.4172789,-0.2036735,	-0.4434014 ,-0.2666667,
    0.4200680  ,0.5387075 ,	0.4883673 , -0.6598639 ;

    ekf::Cmatrix c=filter::Cmatrix::Random();

    const unsigned kmax=1000;

    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector uk[kmax];

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;

    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=filter::Rmatrix::Random()*0.001;

    ekf::Qmatrix q1=filter::Qmatrix::Zero()*0.0001;

    ekf::Rmatrix r=r1*r1.transpose();
    ekf::Qmatrix q=q1*q1.transpose();

    func.setA(a);
    func.setC(c);

    uk[0]=ekf::InputVector::Random();

    for (unsigned k=1; k<=kmax; ++k)
    {
        ekf::StateVector v;
        for (unsigned i=0;i<ekf::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        ekf::MeasureVector w;
        for (unsigned i=0;i<ekf::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        uk[k]=ekf::InputVector::Random();

        xk[k]=x=func.stateDynamics(x,uk[k-1],k-1)+v;
        yk[k-1]=func.measureDynamics(x,uk[k-1],k-1)+w;
    }

    ekf::StateVector xh=filter::StateVector::Random();

    f2.setState(xh,0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();

    for (unsigned i=0;i<filter::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    f2.setStateCovariance(p);

    f2.setA(a);
    f2.setC(c);

    f2.setR(r);
    f2.setQ(q);

    f2.setInput(uk[0],0);

    for (unsigned i=1;i<=kmax;++i)
    {
        f2.setMeasurement(yk[i-1],i);
        f2.setInput(uk[i],i);

        func.reset();
        std::cout << "Success in doNothing." << std::endl;
        //f2.setInput(u,i);
        std::cout<<xk[i].transpose()<<std::endl;
        std::cout<<yk[i-1].transpose()<<std::endl;
        std::cout<<uk[i-1].transpose()<<std::endl;
        std::cout<<f2.getEstimateState(i).transpose()<<std::endl;
        std::cout<<f2.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }
    system("pause");
}


void testKalmanFilter()
{
    filter::Amatrix a;

    a<<	-0.6785714 ,0.1156463 ,	0.4392517 ,	0.2863946 ,
    0.0865306  ,-0.0273810,	0.3355102 ,	0.0184150 ,
    - 0.4172789,-0.2036735,	-0.4434014,	-0.2666667,
    0.4200680  ,0.5387075 ,	0.4883673 , -0.6598639;

    filter::Bmatrix b=filter::Bmatrix::Random();
    filter::Cmatrix c=filter::Cmatrix::Random();
    filter::Dmatrix d=filter::Dmatrix::Random();

    const unsigned kmax=1000;

    filter::StateVector xk[kmax+1];
    filter::MeasureVector yk[kmax];
    filter::InputVector uk[kmax];

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;

    boost::lagged_fibonacci1279 gen_;

    filter::Rmatrix r1=filter::Rmatrix::Random();

    filter::Qmatrix q1=filter::Qmatrix::Random();

    filter::Rmatrix r=r1*r1.transpose();
    filter::Qmatrix q=q1*q1.transpose();

    uk[0]=filter::InputVector::Random();

    for (unsigned k=1; k<=kmax; ++k)
    {

        filter::StateVector v;
        for (unsigned i=0;i<filter::stateSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            v[i]=g(gen_);
        }
        v=q1*v;

        filter::MeasureVector w;
        for (unsigned i=0;i<filter::measureSize;++i)
        {
            boost::normal_distribution<> g(0, 1);
            w[i]=g(gen_);
        }
        w=r1*w;

        uk[k]=filter::InputVector::Random();



        xk[k]=x=a*x+v+b*uk[k-1];
        yk[k-1]=c*x+w+d*uk[k];
    }

    filter::StateVector xh=filter::StateVector::Random();

    f.setState(xh,0);

    filter::Pmatrix p=filter::Pmatrix::Zero();

    for (unsigned i=0;i<filter::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    f.setStateCovariance(p);

    f.setA(a);
    f.setB(b);
    f.setC(c);
    f.setD(d);

    f.setR(r);
    f.setQ(q);

    f.setInput(uk[0],0);

    for (unsigned i=1;i<=kmax;++i)
    {
        //filter::InputVector u=filter::InputVector::Random();
        f.setMeasurement(yk[i-1],i);
        f.setInput(uk[i],i);
        //f.setInput(u,i);
        std::cout<<xk[i].transpose()<<std::endl;
        std::cout<<yk[i-1].transpose()<<std::endl;
        std::cout<<uk[i-1].transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }
    system("pause");
}

int main()
{
    std::cout<<"Starting"<<std::endl;

    testExtendedKalmanFilter();

    return 0;

}
