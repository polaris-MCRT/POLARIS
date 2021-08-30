#pragma once
#include "Typedefs.h"
// m x n
// row    -> m (zeile)
// column |  n (spalte)

typedef vector<vector<vector<double> > > Matrix3D;

#ifndef MATRIX2D
#define MATRIX2D

class Matrix2D
{
  public:
    Matrix2D()
    {
        m_n = 0;
        m_m = 0;
        m_size = 0;
        m_data = 0;
        m_pos = 0;
    }

    Matrix2D(uint m, uint n)
    {
        m_n = n;
        m_m = m;
        m_size = n * m;
        m_pos = 0;
        m_data = new double[m_size];
        for(uint i = 0; i < m_size; i++)
            m_data[i] = 0;
    }

    Matrix2D(uint m, uint n, double * data)
    {
        m_n = n;
        m_m = m;
        m_size = n * m;
        m_pos = 0;
        m_data = new double[m_size];
        for(uint i = 0; i < m_size; i++)
            m_data[i] = data[i];
    }

    Matrix2D(uint m, uint n, vector<double> data)
    {
        m_n = n;
        m_m = m;
        m_size = n * m;
        m_pos = 0;
        m_data = new double[m_size];
        for(uint i = 0; i < m_size; i++)
            m_data[i] = data[i];
    }

    Matrix2D(const Matrix2D & rhs)
    {
        m_n = rhs.col();
        m_m = rhs.row();
        m_size = m_n * m_m;
        m_pos = 0;
        m_data = new double[m_size];

        for(uint i = 0; i < m_size; i++)
            m_data[i] = rhs(i);
    }

    uint get_n()
    {
        return m_n;
    }

    uint get_m()
    {
        return m_m;
    }

    void resize(uint m, uint n)
    {
        if(m_data != 0)
            delete[] m_data;
        m_n = n;
        m_m = m;
        m_size = n * m;
        m_data = new double[m_size];
        for(uint i = 0; i < m_size; i++)
            m_data[i] = 0;
    }

    void replace(double o_val, double n_val)
    {
        for(uint i = 0; i < m_size; i++)
            if(m_data[i] == o_val)
                m_data[i] = n_val;
    }

    void resize(uint m, uint n, double val)
    {
        if(m_data != 0)
            delete[] m_data;
        m_n = n;
        m_m = m;
        m_size = n * m;
        m_data = new double[m_size];
        for(uint i = 0; i < m_size; i++)
            m_data[i] = val;
    }

    ~Matrix2D(void)
    {
        if(m_data != 0)
            delete[] m_data;
    }

    void clear()
    {
        if(m_data != 0)
            delete m_data;

        m_n = 0;
        m_m = 0;
        m_size = 0;
        m_data = 0;
    }

    double & operator()(uint i, uint j)
    {
#ifdef DEBUG
        if(i * m_n + j >= m_size)
        {
            cout << "Matrix overflow" << endl;
            i = j = 0;
        }
#endif
        return m_data[i * m_n + j];
    }

    double & operator()(uint i)
    {
        return m_data[i];
    }

    void operator*=(double val)
    {
        for(uint i = 0; i < m_size; i++)
            m_data[i] *= val;
    }

    void operator+=(double val)
    {
        for(uint i = 0; i < m_size; i++)
            m_data[i] += val;
    }

    void operator+=(const Matrix2D & mat)
    {
        if(m_size == mat.size())
            for(uint i = 0; i < m_size; i++)
                m_data[i] += mat(i);
    }

    void operator-=(const Matrix2D & mat)
    {
        if(m_size == mat.size())
            for(uint i = 0; i < m_size; i++)
                m_data[i] -= mat(i);
    }

    void printMatrix()
    {
        for(uint i = 0; i < m_m; i++)
        {
            for(uint j = 0; j < m_n; j++)
                cout << m_data[i * m_n + j] << " ";

            cout << endl;
        }
    }

    double minElement()
    {
        if(m_size == 0)
            return 0;

        if(m_size == 1)
            return m_data[0];

        double min = m_data[0];

        for(uint i = 0; i < m_size; i++)
            if(m_data[i] < min)
                min = m_data[i];

        return min;
    }

    double maxElement()
    {
        if(m_size == 0)
            return 0;

        if(m_size == 1)
            return m_data[0];

        double max = m_data[0];

        for(uint i = 0; i < m_size; i++)
            if(m_data[i] > max)
                max = m_data[i];

        return max;
    }

    void setValue(uint i, uint j, double val)
    {
#pragma omp atomic write
        m_data[i * m_n + j] = val;
    }

    void addValue(uint i, uint j, double val)
    {
#pragma omp atomic update
        m_data[i * m_n + j] += val;
    }

    void addValue(uint k, double val)
    {
#pragma omp atomic update
        m_data[k] += val;
    }

    double operator()(uint i, uint j) const
    {
        return m_data[i * m_n + j];
    }

    double operator()(uint k) const
    {
#ifdef DEBUG
        if(k >= m_size)
        {
            cout << "Matrix overflow" << endl;
            k = 0;
        }

        if(m_data == 0)
            return 0;
#endif
        return m_data[k];
    }

    void set(uint k, double val)
    {
        m_data[k] = val;
    }

    void set(double * val, uint size)
    {
        if(m_size != size)
            return;

        for(uint k = 0; k < size; k++)
            m_data[k] = val[k];
    }

    void fill(double val)
    {
        for(uint k = 0; k < m_size; k++)
            m_data[k] = val;
    }

    void transpose(bool invert_values)
    {
        if(m_m != m_n)
            cout << "Matrix is not symmetric for transpose!" << endl;

        double * tmp_data = new double[m_size];
        for(uint i = 0; i < m_n; i++)
            for(uint j = 0; j < m_n; j++)
                tmp_data[i * m_n + j] = m_data[j * m_n + i];

        for(uint k = 0; k < m_size; k++)
            if(invert_values)
                m_data[k] = -tmp_data[k];
            else
                m_data[k] = tmp_data[k];

        delete[] tmp_data;
    }

    void transpose()
    {
        double * tmp_data = new double[m_size];
        for(uint i = 0; i < m_n; i++)
            for(uint j = 0; j < m_n; j++)
                tmp_data[i * m_n + j] = m_data[j * m_n + i];

        for(uint k = 0; k < m_size; k++)
            m_data[k] = tmp_data[k];

        delete[] tmp_data;
    }

    dlist get_n_list(uint j)
    {
        dlist res;
        res.resize(m_n);
        for(uint i = 0; i < m_n; i++)
            res[i] = m_data[i * m_n + j];
        return res;
    }

    dlist get_m_list(uint i)
    {
        dlist res;
        res.resize(m_n);
        for(uint j = 0; j < m_m; j++)
            res[j] = m_data[i * m_n + j];
        return res;
    }

    Matrix2D & operator=(const Matrix2D & rhs)
    {
        if(m_size != rhs.size())
            resize(rhs.row(), rhs.col());

        for(uint i = 0; i < m_size; i++)
            m_data[i] = rhs(i);

        return *this;
    }

    Matrix2D & operator=(Matrix2D * rhs)
    {
        cout << "Matrix2D move" << endl;
        if(m_size != rhs->size())
            resize(rhs->row(), rhs->col());

        for(uint i = 0; i < m_size; i++)
            m_data[i] = rhs->operator()(i);

        delete rhs;
        rhs = 0;

        return *this;
    }

    bool operator==(const Matrix2D & rhs)
    {
        return (m_n == rhs.col() && m_m == rhs.row());
    }

    uint size() const
    {
        return m_size;
    }

    uint row() const
    {
        return m_m;
    }

    uint col() const
    {
        return m_n;
    }

    void unityMatrix()
    {
        if(m_m == m_n)
            for(uint i = 0; i < m_size; i++)
            {
                if(i % (m_n + 1) == 0)
                    m_data[i] = 1;
                else
                    m_data[i] = 0;
            }
    }

    void unityMatrix(uint m, uint n)
    {
        if(m_data != 0)
            delete m_data;

        m_n = n;
        m_m = m;
        m_size = n * m;
        m_data = new double[m_size];

        if(m_m == m_n)
            for(uint i = 0; i < m_size; i++)
            {
                if(i % (m_n + 1) == 0)
                    m_data[i] = 1;
                else
                    m_data[i] = 0;
            }
    }

    double sum()
    {
        double sum = 0;
        for(uint i = 0; i < m_size; i++)
            sum += m_data[i];
        return sum;
    }

  private:
    uint m_n; // column
    uint m_m; // row
    uint m_size;
    uint m_pos;
    double * m_data;
};

namespace
{

inline Matrix2D operator*(const Matrix2D & lhs, const Matrix2D & rhs)
{
    Matrix2D res = Matrix2D(lhs.row(), rhs.col());
    if(lhs.col() != rhs.row())
        return res;

    for(uint i = 0; i < res.row(); i++)
        for(uint j = 0; j < res.col(); j++)
        {
            double tmp = 0;
            for(uint k = 0; k < lhs.col(); k++)
                tmp += lhs(i, k) * rhs(k, j);
            res.setValue(i,j,tmp);
        }

    return res;
}

inline Matrix2D operator*(double val, const Matrix2D & mat)
{
    Matrix2D res(mat);
    res *= val;
    return res;
}

inline Matrix2D operator*(const Matrix2D & mat, double val)
{
    Matrix2D res(mat);
    res *= val;
    return res;
}

inline Matrix2D operator+(const Matrix2D & mat, const Matrix2D & rhs)
{
    Matrix2D res(mat.row(), mat.col());

    for(uint i = 0; i < res.size(); i++)
        res.set(i, mat(i) + rhs(i));
    return res;
}

inline Matrix2D operator+(const Matrix2D & mat, double val)
{
    Matrix2D res(mat);
    res += val;
    return res;
}

inline ostream & operator<<(ostream & out, const Matrix2D & mat)
{

    for(uint i = 0; i < mat.row(); i++)
    {
        for(uint j = 0; j < mat.col(); j++)
            out << mat(i, j) << " ";

        out << endl;
    }

    return out;
}

}
#endif
