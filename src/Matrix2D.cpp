/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Matrix2D.hpp"

void Matrix2D::clear()
{
    if(m_data != 0)
        delete m_data;

    m_n = 0;
    m_m = 0;
    m_size = 0;
    m_data = 0;
}

uint Matrix2D::get_n()
{
    return m_n;
}

uint Matrix2D::get_m()
{
    return m_m;
}

void Matrix2D::resize(uint m, uint n)
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

void Matrix2D::replace(double o_val, double n_val)
{
    for(uint i = 0; i < m_size; i++)
        if(m_data[i] == o_val)
            m_data[i] = n_val;
}

void Matrix2D::resize(uint m, uint n, double val)
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

void Matrix2D::printMatrix()
{
    for(uint i = 0; i < m_m; i++)
    {
        for(uint j = 0; j < m_n; j++)
            cout << m_data[i * m_n + j] << " ";

        cout << endl;
    }
}

double Matrix2D::minElement()
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

double Matrix2D::maxElement()
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

void Matrix2D::setValue(uint i, uint j, double val)
{
    #pragma omp atomic write
    m_data[i * m_n + j] = val;
}

void Matrix2D::addValue(uint i, uint j, double val)
{
    #pragma omp atomic update
    m_data[i * m_n + j] += val;
}

void Matrix2D::addValue(uint k, double val)
{
    #pragma omp atomic update
    m_data[k] += val;
}

void Matrix2D::set(uint k, double val)
{
    m_data[k] = val;
}

void Matrix2D::set(double * val, uint size)
{
    if(m_size != size)
        return;

    for(uint k = 0; k < size; k++)
        m_data[k] = val[k];
}

void Matrix2D::fill(double val)
{
    for(uint k = 0; k < m_size; k++)
        m_data[k] = val;
}

void Matrix2D::transpose(bool invert_values)
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

void Matrix2D::transpose()
{
    double * tmp_data = new double[m_size];
    for(uint i = 0; i < m_n; i++)
        for(uint j = 0; j < m_n; j++)
            tmp_data[i * m_n + j] = m_data[j * m_n + i];

    for(uint k = 0; k < m_size; k++)
        m_data[k] = tmp_data[k];

    delete[] tmp_data;
}

dlist Matrix2D::get_n_list(uint j)
{
    dlist res;
    res.resize(m_n);
    for(uint i = 0; i < m_n; i++)
        res[i] = m_data[i * m_n + j];
    return res;
}

dlist Matrix2D::get_m_list(uint i)
{
    dlist res;
    res.resize(m_n);
    for(uint j = 0; j < m_m; j++)
        res[j] = m_data[i * m_n + j];
    return res;
}

uint Matrix2D::size() const
{
    return m_size;
}

uint Matrix2D::row() const
{
    return m_m;
}

uint Matrix2D::col() const
{
    return m_n;
}

void Matrix2D::unityMatrix()
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

void Matrix2D::unityMatrix(uint m, uint n)
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

double Matrix2D::sum()
{
    double sum = 0;
    for(uint i = 0; i < m_size; i++)
        sum += m_data[i];
    return sum;
}

double & Matrix2D::operator()(uint i, uint j)
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

double & Matrix2D::operator()(uint i)
{
    return m_data[i];
}

void Matrix2D::operator*=(double val)
{
    for(uint i = 0; i < m_size; i++)
        m_data[i] *= val;
}

void Matrix2D::operator+=(double val)
{
    for(uint i = 0; i < m_size; i++)
        m_data[i] += val;
}

void Matrix2D::operator+=(const Matrix2D & mat)
{
    if(m_size == mat.size())
        for(uint i = 0; i < m_size; i++)
            m_data[i] += mat(i);
}

void Matrix2D::operator-=(const Matrix2D & mat)
{
    if(m_size == mat.size())
        for(uint i = 0; i < m_size; i++)
            m_data[i] -= mat(i);
}

double Matrix2D::operator()(uint i, uint j) const
{
    return m_data[i * m_n + j];
}

double Matrix2D::operator()(uint k) const
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

Matrix2D & Matrix2D::operator=(const Matrix2D & rhs)
{
    if(m_size != rhs.size())
        resize(rhs.row(), rhs.col());

    for(uint i = 0; i < m_size; i++)
        m_data[i] = rhs(i);

    return *this;
}

Matrix2D & Matrix2D::operator=(Matrix2D * rhs)
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

bool Matrix2D::operator==(const Matrix2D & rhs)
{
    return (m_n == rhs.col() && m_m == rhs.row());
}

Matrix2D operator*(const Matrix2D & lhs, const Matrix2D & rhs)
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

Matrix2D operator*(double val, const Matrix2D & mat)
{
    Matrix2D res(mat);
    res *= val;
    return res;
}

Matrix2D operator*(const Matrix2D & mat, double val)
{
    Matrix2D res(mat);
    res *= val;
    return res;
}

Matrix2D operator+(const Matrix2D & mat, const Matrix2D & rhs)
{
    Matrix2D res(mat.row(), mat.col());

    for(uint i = 0; i < res.size(); i++)
        res.set(i, mat(i) + rhs(i));
    return res;
}

Matrix2D operator+(const Matrix2D & mat, double val)
{
    Matrix2D res(mat);
    res += val;
    return res;
}

ostream & operator<<(ostream & out, const Matrix2D & mat)
{

    for(uint i = 0; i < mat.row(); i++)
    {
        for(uint j = 0; j < mat.col(); j++)
            out << mat(i, j) << " ";

        out << endl;
    }

    return out;
}
