#include "Matrix2D.h"
#include "Typedefs.h"

#ifndef CVECTOR3D
#define CVECTOR3D

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

    void printValues()
    {
        cout << "(" << float(x) << ":" << float(y) << ":" << float(z) << ")";
    }

    void set(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }

    void get(double & _x, double & _y, double & _z)
    {
        _x = x;
        _y = y;
        _z = z;
    }

    // Operators

    Vector3D & operator=(const Vector3D & rhs)
    {
        x = rhs.X();
        y = rhs.Y();
        z = rhs.Z();
        return *this;
    }

    bool operator==(const Vector3D & rhs)
    {
        if(x != rhs.X())
            return false;
        if(y != rhs.Y())
            return false;
        if(z != rhs.Z())
            return false;

        return true;
    }

    Vector3D & operator=(double val)
    {
        x = val;
        y = val;
        z = val;
        return *this;
    }

    double operator()(uint i)
    {
        if(i == 1)
            return y;

        if(i == 2)
            return z;

        return x;
    }

    void set(double val, uint i)
    {
        if(i == 1)
            y = val;

        if(i == 2)
            z = val;

        x = val;
    }

    Vector3D crossProduct(const Vector3D & vec)
    {
        Vector3D tmp;
        double ax = x, bx = vec.X();
        double ay = y, by = vec.Y();
        double az = z, bz = vec.Z();

        tmp.setX(ay * bz - az * by);
        tmp.setY(az * bx - ax * bz);
        tmp.setZ(ax * by - ay * bx);

        return tmp;
    }

    void rot(Vector3D & n, double cos_a, double sin_a)
    {
        Vector3D cx(x, y, z);
        Vector3D cr = cross(n, cx);
        double dot = x * n.X() + y * n.Y() + z * n.Z();

        Vector3D res = n * dot + cross(cr, n) * cos_a + cr * sin_a;

        x = res.X();
        y = res.Y();
        z = res.Z();
    }

    void spher2cart()
    {
        double r = x;
        double ph = y;
        double th = z;

        double sin_ph = sin(ph);
        double cos_ph = cos(ph);

        double sin_th = sin(th);
        double cos_th = cos(th);

        x = r * sin_th * cos_ph;
        y = r * sin_th * sin_ph;
        z = r * cos_th;
    }

    void cart2spher()
    {
        double r = sqrt(x * x + y * y + z * z);

        double arg;
        if(r != 0)
            arg = z / r;
        else
            arg = 1;

        if(arg > 1)
            arg = 1;

        if(arg < -1)
            arg = -1;

        double ph = atan2(y, x), th = acos(arg);

        if(ph < 0)
            ph += PIx2;

        x = r;
        y = ph;
        z = th;
    }

    void cart2cyl()
    {
        double r = sqrt(x * x + y * y);
        double ph = atan2(y, x);

        if(ph < 0.0)
            ph += PIx2;

        x = r;
        y = ph;
    }

    void cyl2cart()
    {
        double r = x;
        double ph = y;

        x = r * cos(ph);
        y = r * sin(ph);
    }

    Vector3D getSphericalCoord() const
    {
        if(x * x + y * y + z * z == 0)
            return Vector3D(0, 0, 0);

        double r = sqrt(x * x + y * y + z * z);

        double arg;
        if(r != 0)
            arg = z / r;
        else
            arg = 1;

        if(arg > 1)
            arg = 1;

        if(arg < -1)
            arg = -1;

        double ph = atan2(y, x), th = acos(arg);

        if(ph < 0)
            ph += PIx2;

        return Vector3D(r, ph, th);
    }

    Vector3D getCylindricalCoord()const
    {
        if(x * x + y * y + z * z == 0)
            return Vector3D(0, 0, 0);

        double r = sqrt(x * x + y * y);
        double ph = atan2(y, x);

        if(ph < 0)
            ph += PIx2;

        return Vector3D(r, ph, z);
    }

    double getPhiCoord() const
    {
        if(x * x + y * y + z * z == 0)
            return 0;

        double ph = atan2(y, x);

        if(ph < 0)
            ph += PIx2;

        return ph;
    }

    void rot(const Vector3D & n, double ang)
    {
        Vector3D cx(x, y, z);
        Vector3D cr = cross(n, cx);
        double dot = x * n.X() + y * n.Y() + z * n.Z();

        double cos_a = cos(ang);
        double sin_a = sin(ang);

        Vector3D res = n * dot + cross(cr, n) * cos_a + cr * sin_a;

        x = res.X();
        y = res.Y();
        z = res.Z();
    }

    void cross(const Vector3D & vec)
    {
        double ax = x, bx = vec.X();
        double ay = y, by = vec.Y();
        double az = x, bz = vec.Z();

        x = ay * bz - az * by;
        y = az * bx - ax * bz;
        z = ax * by - ay * bx;
    }

    Vector3D cross(const Vector3D & vec1, const Vector3D & vec2)
    {
        double ax = vec1.X(), bx = vec2.X();
        double ay = vec1.Y(), by = vec2.Y();
        double az = vec1.Z(), bz = vec2.Z();

        return Vector3D(ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx);
    }

    Vector3D operator+(const Vector3D & rhs) const
    {
        return Vector3D(x + rhs.X(), y + rhs.Y(), z + rhs.Z());
    }

    Vector3D operator-(const Vector3D & rhs) const
    {
        return Vector3D(x - rhs.X(), y - rhs.Y(), z - rhs.Z());
    }

    Vector3D operator-()
    {
        return Vector3D(-x, -y, -z);
    }

    Vector3D & operator+=(const Vector3D & rhs)
    {
        x += rhs.X();
        y += rhs.Y();
        z += rhs.Z();
        return *this;
    }

    Vector3D & operator-=(const Vector3D & rhs)
    {
        x -= rhs.X();
        y -= rhs.Y();
        z -= rhs.Z();
        return *this;
    }

    bool operator!=(const int val)
    {
        if(x == val && y == val && z == val)
        {
            return false;
        }
        return true;
    }

    Vector3D operator*(double val) const
    {
        Vector3D tmp(x * val, y * val, z * val);
        return tmp;
    }

    Vector3D operator/(double val)
    {
        Vector3D tmp(x / val, y / val, z / val);
        return tmp;
    }

    void operator*=(double v)
    {
        x *= v;
        y *= v;
        z *= v;
    }

    void operator*=(float v)
    {
        x *= (double)v;
        y *= (double)v;
        z *= (double)v;
    }

    void operator*=(int v)
    {
        x *= (double)v;
        y *= (double)v;
        z *= (double)v;
    }

    void operator*=(uint v)
    {
        x *= (double)v;
        y *= (double)v;
        z *= (double)v;
    }

    void operator/=(double v)
    {
        x /= v;
        y /= v;
        z /= v;
    }

    void operator/=(float v)
    {
        x /= (double)v;
        y /= (double)v;
        z /= (double)v;
    }

    void operator/=(int v)
    {
        x /= (double)v;
        y /= (double)v;
        z /= (double)v;
    }

    void operator/=(uint v)
    {
        x /= (double)v;
        y /= (double)v;
        z /= (double)v;
    }

    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    double sq_length() const
    {
        return (x * x + y * y + z * z);
    }

    double sqrel_length(Vector3D & c) const
    {
        return (x - c.X()) * (x - c.X()) + (y - c.Y()) * (y - c.Y()) + (z - c.Z()) * (z - c.Z());
    }

    // random direction according to http://mathworld.wolfram.com/SpherePointPicking.html

    void rndDir(double r1, double r2)
    {
        double u = 1.0 - 2.0 * r1;
        double th = PIx2 * r2;
        double sqr = sqrt(1.0 - u * u);
        x = sqr * cos(th);
        y = sqr * sin(th);
        z = u;
        normalize();
    }

    void normalize()
    {
        double len = length();
        if(len != 0)
        {
            x /= len;
            y /= len;
            z /= len;
        }
    }

    Vector3D normalized()
    {
        Vector3D tmp;
        double len = length();

        if(len != 0)
            tmp.set(x / len, y / len, z / len);
        return tmp;
    }

    void setX(double _x)
    {
        x = _x;
    }

    void setY(double _y)
    {
        y = _y;
    }

    void setZ(double _z)
    {
        z = _z;
    }

    void setR(double _R)
    {
        x = _R;
    }

    void setPhi(double _phi)
    {
        y = _phi;
    }

    void setTheta(double _theta)
    {
        z = _theta;
    }

    double X() const
    {
        return x;
    }

    double Y() const
    {
        return y;
    }

    double Z() const
    {
        return z;
    }

    double R() const
    {
        return x;
    }

    double Phi() const
    {
        return y;
    }

    double Theta() const
    {
        return z;
    }

    void rotX(double phi)
    {
        double ty = y, tz = z;
        double s = sin(phi), c = cos(phi);
        y = ty * c - tz * s;
        z = ty * s + tz * c;
    }

    void rotY(double phi)
    {
        double tx = x, tz = z;
        double s = sin(phi), c = cos(phi);
        x = tx * c - tz * s;
        z = tx * s + tz * c;
    }

    void rotZ(double phi)
    {
        double tx = x, ty = y;
        double s = sin(phi), c = cos(phi);
        x = tx * c - ty * s;
        y = tx * s + ty * c;
    }

  private:
    double x, y, z;
};

namespace
{

double sign(double x)
{
    if(x < 0)
        return -1;
    if(x > 0)
        return 1;

    return 0;
}

bool ctob(char c)
{
    if(c == '1')
        return true;

    return false;
}

double angle(double x, double y)
{
    if(x == 0 && y == 0)
        return 0;

    if(x == 0)
    {
        if(y >= 0)
            return PI / 2.0;
        else
            return 3.0 * PI / 2.0;
    }

    double p = atan(abs(y) / abs(x));

    if(y >= 0 && x < 0)
        return PI - p;
    if(y <= 0 && x < 0)
        return PI + p;
    if(y < 0 && x > 0)
        return 2.0 * PI - p;

    return p;
}

double atan3(double x, double y)
{
    double angle = atan2(y, x);
    if(angle > 0)
        return angle;
    else
        return angle + PIx2;
}

// inline Vector3D operator-(const Vector3D & lhs,const Vector3D & rhs) { return
// Vector3D(lhs.X()-rhs.X(),lhs.X()-rhs.Y(),lhs.X()-rhs.Z());};

inline Vector3D operator*(double val, const Vector3D & rhs)
{
    return Vector3D(rhs.X() * val, rhs.Y() * val, rhs.Z() * val);
}

inline Vector3D operator*(const Matrix2D & dM, const Vector3D & v)
{
    Vector3D tmp;

#ifdef DEBUG
    if(dM.col() != 3 || dM.row() != 3)
        return tmp;
#endif

    tmp.setX(v.X() * dM(0, 0) + v.Y() * dM(0, 1) + v.Z() * dM(0, 2));
    tmp.setY(v.X() * dM(1, 0) + v.Y() * dM(1, 1) + v.Z() * dM(1, 2));
    tmp.setZ(v.X() * dM(2, 0) + v.Y() * dM(2, 1) + v.Z() * dM(2, 2));

    return tmp;
}

inline double operator*(const Vector3D & lhs, const Vector3D & rhs)
{
    return lhs.X() * rhs.X() + lhs.Y() * rhs.Y() + lhs.Z() * rhs.Z();
}

inline ostream & operator<<(ostream & out, const Vector3D & rhs)
{
    out << " x: " << rhs.X() << " y: " << rhs.Y() << " z: " << rhs.Z() << " length: " << rhs.length() << " ";

    return out;
}

double getAngleTheta(Vector3D lhs, Vector3D rhs)
{
    if(rhs * lhs == 0.0)
        return PI2;

    double arg = (lhs.X() * rhs.X() + lhs.Y() * rhs.Y() + lhs.Z() * rhs.Z()) / (lhs.length() * rhs.length());

    if(arg < -1.0)
        arg = -1.0;

    if(arg > 1.0)
        arg = 1.0;

    return abs(acos(arg));
    // return abs(0.5*PI-acos((lhs.X()*rhs.X()+lhs.Y()*rhs.Y()+lhs.Z()*rhs.Z())
    //	/(lhs.length()*rhs.length()) ));
}

Vector3D projection(Vector3D v, Vector3D w)
{
    double tmp = w.sq_length();
    tmp = (v * w) / tmp;

    Vector3D res = tmp * w;
    return res;
}

double getAnglePhi(Vector3D ex, Vector3D ey, Vector3D proj)
{
    Vector3D px, py;
    double lx, ly, sx, sy;

    px = (proj.X() * ex.X() + proj.Y() * ex.Y() + proj.Z() * ex.Z()) / ex.length() * ex;
    py = (proj.X() * ey.X() + proj.Y() * ey.Y() + proj.Z() * ey.Z()) / ey.length() * ey;

    if(0 == px.X() && 0 == ex.X())
    {
        if(0 == px.Y() && 0 == ex.Y())
            sx = sign(px.Z() / ex.Z());
        else
            sx = sign(px.Y() / ex.Y());
    }
    else
        sx = sign(px.X() / ex.X());

    if(0 == py.X() && 0 == ey.X())
    {
        if(0 == py.Y() && 0 == ey.Y())
            sy = sign(py.Z() / ey.Z());
        else
            sy = sign(py.Y() / ey.Y());
    }
    else
        sy = sign(py.X() / ey.X());

    lx = sx * px.length();
    ly = sy * py.length();

    return angle(lx, ly);
}

Vector3D cross(const Vector3D & vec1, const Vector3D & vec2)
{
    double ax = vec1.X(), bx = vec2.X();
    double ay = vec1.Y(), by = vec2.Y();
    double az = vec1.Z(), bz = vec2.Z();

    return Vector3D(ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx);
}

double rot_sgn(Vector3D n, Vector3D o_d, Vector3D n_d)
{
    Vector3D o = n.crossProduct(o_d);
    double tmp = o.sq_length();

    if(tmp == 0)
        return 0;

    tmp = (n_d.X() * o.X() + n_d.Y() * o.Y() + n_d.Z() * o.Z()) / tmp;
    Vector3D op(tmp * o.X(), tmp * o.Y(), tmp * o.Z());

    if(op.X() != 0)
    {
        return sign(o.X() / op.X());
    }
    else
    {
        if(op.Y() != 0)
        {
            return sign(o.Y() / op.Y());
        }
        else
        {
            if(op.Z() != 0)
            {
                return sign(o.Z() / op.Z());
            }
        }
    }

    return 0;
}
};

#endif
