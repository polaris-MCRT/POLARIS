/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CVECTOR3D_H
#define CVECTOR3D_H

#include "Matrix2D.hpp"
#include "Typedefs.hpp"

class Vector3D
{
public:
    Vector3D(void)
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector3D(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    ~Vector3D(void)
    {}

    void printValues() const;

    void set(double _x, double _y, double _z);

    void get(double & _x, double & _y, double & _z);

    void set(double val, uint i);

    Vector3D crossProduct(const Vector3D & vec);

    void rot(const Vector3D & n, double cos_a, double sin_a);

    void spher2cart();

    void cart2spher();

    void cart2cyl();

    void cyl2cart();

    Vector3D getSphericalCoord() const;

    Vector3D getCylindricalCoord() const;

    double getPhiCoord() const;

    void rot(const Vector3D & n, double ang);

    void cross(const Vector3D & vec);

    Vector3D cross(const Vector3D & vec1, const Vector3D & vec2);

    double length() const;

    double sq_length() const;

    double sqrel_length(Vector3D & c) const;

    // random direction according to http://mathworld.wolfram.com/SpherePointPicking.html
    void rndDir(double r1, double r2, uint exponentThetaBias = 1);

    // for TRUST benchmark restrict photon directions to increase SNR
    void rndDirTRUST(double r1, double r2);

    void normalize();

    Vector3D normalized();

    void setX(double _x);

    void setY(double _y);

    void setZ(double _z);

    void setR(double _R);

    void setPhi(double _phi);

    void setTheta(double _theta);

    double X() const;

    double Y() const;

    double Z() const;

    double R() const;

    double Phi() const;

    double Theta() const;

    void rotX(double phi);

    void rotY(double phi);

    void rotZ(double phi);

    static double sign(double x);

    bool ctob(char c);

    static double angle(double x, double y);

    static double atan3(double x, double y);

    static double getAngleTheta(Vector3D lhs, Vector3D rhs);
    
    static double getAngleThetaOff(Vector3D lhs, Vector3D rhs);

    Vector3D projection(const Vector3D & v, const Vector3D & w);

    static double getAnglePhi(const Vector3D & ex, const Vector3D & ey, const Vector3D & proj);

    double rot_sgn(Vector3D n, Vector3D o_d, Vector3D n_d);

    Vector3D & operator=(const Vector3D & rhs);

    bool operator==(const Vector3D & rhs);

    Vector3D & operator=(double val);

    double operator()(uint i);

    Vector3D operator+(const Vector3D & rhs) const;

    Vector3D operator-(const Vector3D & rhs) const;

    Vector3D operator-();

    Vector3D & operator+=(const Vector3D & rhs);

    Vector3D & operator-=(const Vector3D & rhs);

    bool operator!=(const int val);

    Vector3D operator*(double val) const;

    Vector3D operator/(double val) const;

    void operator*=(double v);

    void operator*=(float v);

    void operator*=(int v);

    void operator*=(uint v);

    void operator/=(double v);

    void operator/=(float v);

    void operator/=(int v);

    void operator/=(uint v);

    friend const Vector3D operator*(double val, const Vector3D & rhs);

    friend Vector3D operator*(const Matrix2D & dM, const Vector3D & v);

    friend double operator*(const Vector3D & lhs, const Vector3D & rhs);

    friend ostream & operator<<(ostream & out, const Vector3D & rhs);

private:
    double x, y, z;
};

#endif /* CVECTOR3D_H */
