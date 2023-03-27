#pragma once

#include <iostream>
#include <functional>

template<typename T = double>
class ABMVector
{
    T* data;
    const uint32_t ndata;

public:
    ABMVector(uint32_t ndata, T value = T{}) : ndata(ndata)
    {
        data = new T[ndata];

        uint32_t i;
        for (i = 0; i < ndata; i++)
            data[i] = value;
    }

    ABMVector(std::initializer_list<T> values) : ndata(values.size())
    {
        data = new T[values.size()];
        std::copy(values.begin(), values.end(), data);
    }

    ABMVector(std::vector<T>& values) : ndata(values.size())
    {
        data = new T[values.size()];
        std::copy(values.begin(), values.end(), data);
    }

    ABMVector(const T* values, uint32_t ndata) : ndata(ndata)
    {
        data = new T[ndata];
        std::copy(&values[0], &values[0] + ndata, data);
    }

    ABMVector(ABMVector& o) : ndata(o.ndata)
    {
        data = new T[o.ndata];
        std::copy(&o.data[0], &o.data[0] + o.ndata, data);
    }

    ABMVector(const ABMVector& o) : ndata(o.ndata)
    {
        data = new T[o.ndata];
        std::copy(&o.data[0], &o.data[0] + o.ndata, data);
    }

    ABMVector(ABMVector&& o) : ndata(o.ndata)
    {
        data = new T[o.ndata];
        std::copy(&o.data[0], &o.data[0] + o.ndata, data);
    }

    ABMVector() = delete; // No default constructor                   

    ~ABMVector()
    {
        delete[] data;
    }

    uint32_t length() const
    {
        return ndata;
    }

    T& operator [](uint32_t i)
    {
        return data[i];
    }

    T operator [](uint32_t i) const
    {
        return data[i];
    }

    void add(T value)
    {
        uint32_t i;
        for (i = 0; i < ndata; i++)
            data[i] += value;
    }
    
    void substract(T value)
    {
        uint32_t i;
        for (i = 0; i < ndata; i++)
            data[i] -= value;
    }

    void multiply(T value)
    {
        uint32_t i;
        for (i = 0; i < ndata; i++)
            data[i] *= value;
    }

    void divide(T value)
    {
        uint32_t i;
        for (i = 0; i < ndata; i++)
            data[i] /= value;
    }

    ABMVector<T> operator +(T value) const
    {
        ABMVector<T> result(*this);
        result.add(value);
        return result;
    }

    ABMVector<T> operator -(T value) const
    {
        ABMVector<T> result(*this);
        result.substract(value);
        return result;
    }

    ABMVector<T> operator *(T value) const
    {
        ABMVector<T> result(*this);
        result.multiply(value);
        return result;
    }

    ABMVector<T> operator /(T value) const
    {
        ABMVector<T> result(*this);
        result.divide(value);
        return result;
    }        
};

// Adams-Bashforth-Moulton order 4 solver with order 5 correction + RK4 initialization
class ABMSolver
{
public:

    ABMSolver(double stepSize, const ABMVector<std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)>>& eqs) : eqs(eqs)
    {
        h = stepSize;
        nVariables = eqs.length();
    }

    ABMVector<double> RK4(double t, const ABMVector<double>& y0, const ABMVector<double>& x) const
    {
        ABMVector<double> result(y0);

        double* k = new double[4];
        constexpr double sixth = static_cast<double>(1.0/6.0);

        uint32_t i;       
        for (i = 0; i < nVariables; i++)
        {           
            k[0] = eqs[i](t,         y0             , x);
            k[1] = eqs[i](t + 0.5*h, y0 + 0.5*k[0]*h, x);
            k[2] = eqs[i](t + 0.5*h, y0 + 0.5*k[1]*h, x);
            k[3] = eqs[i](t +     h, y0 +     k[2]*h, x);

            result[i] += sixth * h*(k[0] + 2*k[1] + 2*k[2] + k[3]);           
        }

        delete[] k;
        return result;
    }

private:

    uint32_t nVariables;
    double h;
    const ABMVector <std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)>> eqs;
  
};