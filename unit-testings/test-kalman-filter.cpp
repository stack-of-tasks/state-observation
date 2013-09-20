#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <boost/random.hpp>

#include <algorithm>

#include <state-observer/kalman-filter.hpp>
#include <state-observer/extended-kalman-filter.hpp>





void testExtendedKalmanFilter()
{
    const static unsigned kmax=1000;


    typedef observation::ExtendedKalmanFilter<4,3,1> ekf;


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

            ekf::StateVector xk1;
            xk1=a_*x+cos(10*(x.transpose()*x)[0])*s_+t_+(u*u.transpose())[0]*s_;
            return xk1;
        }

        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            yk=c_*x+cos(10*(x.transpose()*x)[0])*m_+n_;
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

private:
        ///Special instructions to have a static-sized eigen vector as a member
        enum { NeedsToAlign = ((sizeof(ekf::StateVector)%16)==0)||
                              ((sizeof(ekf::MeasureVector)%16)==0)||
                              ((sizeof(ekf::Amatrix)%16)==0)||
                              ((sizeof(ekf::Cmatrix)%16)==0)
             };
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        ekf::StateVector s_;
        ekf::StateVector t_;
        ekf::MeasureVector m_;
        ekf::MeasureVector n_;

        ekf::Amatrix a_;
        ekf::Cmatrix c_;

    };


    ekf f2;

    ekf::Amatrix a;
    KalmanFunctor func;

    f2.setFunctor(&func);


    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector uk[kmax+1];

    ekf::StateVector x=ekf::StateVector::Zero();

    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=ekf::Rmatrix::Random()*0.01;

    ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.01;

    ekf::Rmatrix r=r1*r1.transpose();
    ekf::Qmatrix q=q1*q1.transpose();

    xk[0]=x;
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

        x=func.stateDynamics(x,uk[k-1],k-1)+v;

        xk[k]=x;
        yk[k-1]=func.measureDynamics(x,uk[k],k)+w;

    }

    ekf::StateVector xh=ekf::StateVector::Random();

    f2.setState(xk[0],0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();



    f2.setStateCovariance(p);

    f2.setR(r);
    f2.setQ(q);

    f2.setInput(uk[0],0);

    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;

    std::cout<<"i "<<0<<std::endl;
    std::cout<<"x "<<xk[0].transpose()<<std::endl;
    std::cout<<"uk "<<uk[0].transpose()<<std::endl;
    std::cout<<"xh "<<f2.getEstimateState(0).transpose()<<std::endl;
    std::cout<<"er "<<f2.getEstimateState(0).transpose()-xk[0].transpose()<<std::endl<<std::endl;

    for (unsigned i=1;i<=kmax;++i)
    {
        f2.setMeasurement(yk[i-1],i);
        f2.setInput(uk[i],i);


        ////Debug code (to be removed)
        //KalmanFunctor* debug=&func;
        //func.reset();
        //debug->reset();

        std::cout<<"i "<<i<<std::endl;
        std::cout<<"yk "<<yk[i-1].transpose()<<std::endl;
        std::cout<<"uk "<<uk[i-1].transpose()<<std::endl;

        ekf::Amatrix a=f2.getAMatrixFD(dx);
        ekf::Cmatrix c= f2.getCMatrixFD(dx);

        std::cout<<"aFD "<<std::endl<<a<<std::endl;
        std::cout<<"cFD "<<std::endl<<c<<std::endl;

        f2.setA(a);
        f2.setC(c);

        std::cout<<"x "<<xk[i].transpose()<<std::endl;
        std::cout<<"xh "<<f2.getEstimateState(i).transpose()<<std::endl;
        std::cout<<"er "<<f2.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }

}

void testExtendedKalmanFilterLTV()
{
    const static unsigned kmax=1000;

    typedef observation::ExtendedKalmanFilter<4,3,1> ekf;


    struct KalmanFunctorLTV:
                public ekf::DynamicsFunctorBase
    {

    public:
        KalmanFunctorLTV()
        {
            s_=ekf::StateVector::Random();
            n_=ekf::MeasureVector::Random();

            for(unsigned i=0;i<=kmax;++i)
            {
                a.push_back( ekf::Amatrix::Random()*0.5 );
                c.push_back( ekf::Cmatrix::Random() );
            }

        }

        virtual ekf::StateVector stateDynamics(const ekf::StateVector& x, const ekf::InputVector& u, unsigned k)
        {
            ekf::StateVector xk1;
            unsigned kk=std::min(k,kmax);
            xk1=a[kk]*x+(u*u.transpose())[0]*s_;

            return xk1;
        }

        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            unsigned kk=std::min(k,kmax);
            yk=c[kk]*x+(u*u.transpose())[0]*n_;
            return yk;
        }

        std::vector<ekf::Amatrix,Eigen::aligned_allocator<ekf::Amatrix> > a;
        std::vector<ekf::Cmatrix,Eigen::aligned_allocator<ekf::Cmatrix> > c;
    private:
        ekf::StateVector s_;
        ekf::MeasureVector n_;
    };


    ekf f2;

    ekf::Amatrix a;
    KalmanFunctorLTV func;

    f2.setFunctor(&func);



    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector uk[kmax+1];

    ekf::StateVector x=ekf::StateVector::Zero();



    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=ekf::Rmatrix::Random()*1;

    ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.1;

    ekf::Rmatrix r=r1*r1.transpose();
    ekf::Qmatrix q=q1*q1.transpose();

    xk[0]=x;
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

        x=func.stateDynamics(x,uk[k-1],k-1)+v;

        xk[k]=x;
        yk[k-1]=func.measureDynamics(x,uk[k],k)+w;

    }

    ekf::StateVector xh=ekf::StateVector::Random();

    f2.setState(xh,0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();



    f2.setStateCovariance(p);

    f2.setR(r);
    f2.setQ(q);

    f2.setInput(uk[0],0);

    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;

    std::cout<<"i "<<0<<std::endl;
    std::cout<<"x "<<xk[0].transpose()<<std::endl;
    std::cout<<"uk "<<uk[0].transpose()<<std::endl;
    std::cout<<"xh "<<f2.getEstimateState(0).transpose()<<std::endl;
    std::cout<<"er "<<f2.getEstimateState(0).transpose()-xk[0].transpose()<<std::endl<<std::endl;

    for (unsigned i=1;i<=kmax;++i)
    {
        f2.setMeasurement(yk[i-1],i);
        f2.setInput(uk[i],i);


        ////Debug code (to be removed)
        //KalmanFunctor* debug=&func;
        //func.reset();
        //debug->reset();

        std::cout<<"i "<<i<<std::endl;
        std::cout<<"yk "<<yk[i-1].transpose()<<std::endl;
        std::cout<<"uk "<<uk[i-1].transpose()<<std::endl;

        ekf::Amatrix a=f2.getAMatrixFD(dx);
        ekf::Cmatrix c= f2.getCMatrixFD(dx);

        f2.setA(a);
        f2.setC(c);


        std::cout<<"aFD "<<std::endl<<a<<std::endl;
        std::cout<<"cFD "<<std::endl<<c<<std::endl;
        std::cout<<"a "<<std::endl<<func.a[i-1]<<std::endl;
        std::cout<<"c "<<std::endl<<func.c[i-1]<<std::endl;
        std::cout<<"x "<<xk[i].transpose()<<std::endl;
        std::cout<<"xh "<<f2.getEstimateState(i).transpose()<<std::endl;
        std::cout<<"er "<<f2.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }

}

void testExtendedKalmanFilterZeroInput()
{
    const static unsigned kmax=1000;


    typedef observation::ExtendedKalmanFilter<4,3> ekf;


    class KalmanFunctor:
                public ekf::DynamicsFunctorBase
    {

    public:
        KalmanFunctor()
        {
            m_=ekf::MeasureVector::Random()*0.1;
            t_=ekf::StateVector::Random();
            n_=ekf::MeasureVector::Random();
        }

        virtual ekf::StateVector stateDynamics(const ekf::StateVector& x, const ekf::InputVector& u, unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::StateVector xk1;
            xk1=a_*x+cos(10*(x.transpose()*x)[0])*s_+t_;
            return xk1;
        }

        virtual ekf::MeasureVector measureDynamics(const ekf::StateVector& x, const ekf::InputVector &u,unsigned k)
        {
            (void)k;//unused
            (void)u;//unused

            ekf::MeasureVector yk;
            yk=c_*x+cos(10*(x.transpose()*x)[0])*m_+n_;
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

private:
        ///Special instructions to have a static-sized eigen vector as a member
        enum { NeedsToAlign = ((sizeof(ekf::StateVector)%16)==0)||
                              ((sizeof(ekf::MeasureVector)%16)==0)||
                              ((sizeof(ekf::Amatrix)%16)==0)||
                              ((sizeof(ekf::Cmatrix)%16)==0)
             };
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)

        ekf::StateVector s_;
        ekf::StateVector t_;
        ekf::MeasureVector m_;
        ekf::MeasureVector n_;

        ekf::Amatrix a_;
        ekf::Cmatrix c_;

    };


    ekf f2;

    ekf::Amatrix a;
    KalmanFunctor func;

    f2.setFunctor(&func);


    ekf::StateVector xk[kmax+1];
    ekf::MeasureVector yk[kmax];
    ekf::InputVector u= ekf::InputVector::Zero();

    ekf::StateVector x=ekf::StateVector::Zero();

    boost::lagged_fibonacci1279 gen_;

    ekf::Rmatrix r1=ekf::Rmatrix::Random()*0.01;

    ekf::Qmatrix q1=ekf::Qmatrix::Random()*0.01;

    ekf::Rmatrix r=r1*r1.transpose();
    ekf::Qmatrix q=q1*q1.transpose();

    xk[0]=x;


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

        x=func.stateDynamics(x,u,k-1)+v;

        xk[k]=x;
        yk[k-1]=func.measureDynamics(x,u,k)+w;

    }

    ekf::StateVector xh=ekf::StateVector::Random();

    f2.setState(xk[0],0);

    ekf::Pmatrix p=ekf::Pmatrix::Zero();



    f2.setStateCovariance(p);

    f2.setR(r);
    f2.setQ(q);

    for (unsigned i=0;i<ekf::stateSize;++i)
    {
        p(i,i)=xh[i];
    }
    p=p*p.transpose();

    ekf::StateVector dx=ekf::StateVector::Constant(1)*1e-8;

    std::cout<<"i "<<0<<std::endl;
    std::cout<<"x "<<xk[0].transpose()<<std::endl;
    std::cout<<"xh "<<f2.getEstimateState(0).transpose()<<std::endl;
    std::cout<<"er "<<f2.getEstimateState(0).transpose()-xk[0].transpose()<<std::endl<<std::endl;

    for (unsigned i=1;i<=kmax;++i)
    {
        f2.setMeasurement(yk[i-1],i);


        ////Debug code (to be removed)
        //KalmanFunctor* debug=&func;
        //func.reset();
        //debug->reset();

        std::cout<<"i "<<i<<std::endl;
        std::cout<<"yk "<<yk[i-1].transpose()<<std::endl;

        ekf::Amatrix a=f2.getAMatrixFD(dx);
        ekf::Cmatrix c= f2.getCMatrixFD(dx);

        std::cout<<"aFD "<<std::endl<<a<<std::endl;
        std::cout<<"cFD "<<std::endl<<c<<std::endl;

        f2.setA(a);
        f2.setC(c);

        std::cout<<"x "<<xk[i].transpose()<<std::endl;
        std::cout<<"xh "<<f2.getEstimateState(i).transpose()<<std::endl;
        std::cout<<"er "<<f2.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }

}


void testKalmanFilter()
{

    typedef observation::KalmanFilter<4,3,2> filter;


    filter f;
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
    filter::InputVector uk[kmax+1];

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;

    boost::lagged_fibonacci1279 gen_;

    filter::Rmatrix r1=filter::Rmatrix::Random()*0.1;

    filter::Qmatrix q1=filter::Qmatrix::Random()*0.1;

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
        f.setMeasurement(yk[i-1],i);
        f.setInput(uk[i],i);

        std::cout<<xk[i].transpose()<<std::endl;
        std::cout<<yk[i-1].transpose()<<std::endl;
        std::cout<<uk[i-1].transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }
}

void testKalmanFilterZeroInput()
{
    typedef observation::KalmanFilter<4,3> filter;

    filter f;
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

    filter::StateVector x=filter::StateVector::Zero();

    xk[0]=x;

    boost::lagged_fibonacci1279 gen_;

    filter::Rmatrix r1=filter::Rmatrix::Random()*0.1;

    filter::Qmatrix q1=filter::Qmatrix::Random()*0.1;

    filter::Rmatrix r=r1*r1.transpose();
    filter::Qmatrix q=q1*q1.transpose();

    filter::InputVector u=filter::InputVector::Zero();

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

        xk[k]=x=a*x+v;
        yk[k-1]=c*x+w;
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

    for (unsigned i=1;i<=kmax;++i)
    {
        f.setMeasurement(yk[i-1],i);

        std::cout<<xk[i].transpose()<<std::endl;
        std::cout<<yk[i-1].transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()<<std::endl;
        std::cout<<f.getEstimateState(i).transpose()-xk[i].transpose()<<std::endl<<std::endl;
    }
}

int main()
{
    std::cout<<"Starting"<<std::endl;

    //testKalmanFilter();
    testKalmanFilterZeroInput();
    //testExtendedKalmanFilter();
    //testExtendedKalmanFilterLTV();
    //testExtendedKalmanFilterZeroInput();

    return 0;

}
