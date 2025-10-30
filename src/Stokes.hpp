/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef STOKESVECTOR_H
#define STOKESVECTOR_H

#include "Matrix2D.hpp"
#include "Typedefs.hpp"

class StokesVector
{
public:
    StokesVector()
    {
        sI = 0;
        sQ = 0;
        sU = 0;
        sV = 0;
        sT = 0;
        sSp1 = 0;
        sSp2 = 0;
        sSp3 = 0;
        sSp4 = 0;
    }

    StokesVector(double val)
    {
        sI = val;
        sQ = val;
        sU = val;
        sV = val;
        sT = val;
        sSp1 = val;
        sSp2 = val;
        sSp3 = val;
        sSp4 = val;
    }

    StokesVector(double I, double Q, double U, double V)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = 0;
        sSp1 = 0;
        sSp2 = 0;
        sSp3 = 0;
        sSp4 = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp1 = 0;
        sSp2 = 0;
        sSp3 = 0;
        sSp4 = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T, double Sp)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp1 = Sp;
        sSp2 = 0;
        sSp3 = 0;
        sSp4 = 0;
    }
    
    StokesVector(double I, double Q, double U, double V, double T, 
            double Sp1, double Sp2, double Sp3, double Sp4)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp1 = Sp1;
        sSp2 = Sp2;
        sSp3 = Sp3;
        sSp4 = Sp4;
    }

    StokesVector(const StokesVector & st)
    {
        sI = st.I();
        sQ = st.Q();
        sU = st.U();
        sV = st.V();
        sT = st.T();
        sSp1 = st.Sp1();
        sSp2 = st.Sp2();
        sSp3 = st.Sp3();
        sSp4 = st.Sp4();
    }

    ~StokesVector(void)
    {}

    // linearly polarized intensity
    double iPol() const;

    // totally polarized intensity
    double tPol() const;

    // degree of linear polarization
    double linPol() const;

    // degree of circular polarization
    double circPol() const;

    // polarization angle
    double getAngle() const;

    void setI(double _I);

    void setQ(double _Q);

    void setU(double _U);

    void setV(double _V);

    void setT(double _T);

    void setSp1(double _Sp1);
    void setSp2(double _Sp2);
    void setSp3(double _Sp3);
    void setSp4(double _Sp4);

    void set(double _I, double _Q, double _U, double _V, double _T);

    void set(double _I, double _Q, double _U, double _V);

    void set(const StokesVector & _S);

    void addI(double _I);

    void addQ(double _Q);

    void addU(double _U);

    void addV(double _V);

    void addT(double _T);

    void addSp1(double _Sp1);
    void addSp2(double _Sp2);
    void addSp3(double _Sp3);
    void addSp4(double _Sp4);

    void addS(StokesVector _S);

    void multI(double _I);

    void multQ(double _Q);

    void multU(double _U);

    void multV(double _V);

    void multT(double _T);

    void multSp1(double _Sp1);
    void multSp2(double _Sp2);
    void multSp3(double _Sp3);
    void multSp4(double _Sp4);

    void multStokesParam(double _S);

    double I() const;

    double Q() const;

    double U() const;

    double V() const;

    double T() const;

    double Sp1() const;
    double Sp2() const;
    double Sp3() const;
    double Sp4() const;

    void rot(double phi);

    void rot(double sin_phi, double cos_phi);

    bool isConsistent();

    void normalize();

    void clear();

    void resetIntensity();

    void depolarize();

    StokesVector & operator=(const StokesVector & ex);

    StokesVector & operator=(double v);

    StokesVector operator+(const StokesVector & ex) const;

    StokesVector operator-(const StokesVector & ex) const;

    StokesVector & operator+=(const StokesVector & ex);

    StokesVector & operator-=(const StokesVector & ex);

    StokesVector & operator*=(double val);

    StokesVector & operator*=(const Matrix2D & dM);

    StokesVector & operator/=(double val);

    friend ostream & operator<<(ostream & out, const StokesVector & ex);

    friend StokesVector operator*(const Matrix2D & dM, const StokesVector & v);

    friend StokesVector operator*(const StokesVector & v, const StokesVector & u);

    friend StokesVector operator*(const StokesVector & v, double val);

    friend StokesVector operator/(const StokesVector & v, double val);

private:
    //stokes parameters
    double sI;
    double sQ;
    double sU;
    double sV; 
    
    //optical depth
    double sT;
    
    //special parameters e.g. column density
    double sSp1;
    double sSp2;
    double sSp3;
    double sSp4;
};

#endif /* STOKESVECTOR_H */
