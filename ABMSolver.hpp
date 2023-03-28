#pragma once

#include <iostream>
#include <functional>
#include <deque>
#include <execution>

#define AB_CONST1  2.291666666667
#define AB_CONST2 -2.458333333333
#define AB_CONST3  1.541666666667
#define AB_CONST4 -0.375000000000
#define AM_CONST1  0.348611111111
#define AM_CONST2  0.897222222222
#define AM_CONST3 -0.366666666667
#define AM_CONST4  0.147222222222
#define AM_CONST5 -0.026388888889
#define SIXTH      0.166666666667

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

// Adams-Bashforth 4-step solver with 5-step Adams-Moulton correction + RK4 initialization [fixed-step PEC ode45 solver]
class ABMSolver
{
public:

    ABMSolver(double stepSize, const ABMVector<std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)>>& eqs, const ABMVector<double>& y0, const ABMVector<double>& x0) : eqs(eqs)
    {
        h = stepSize;
        nVariables = eqs.length();

        x = new std::deque<ABMVector<double>>(1, x0);
        y = new std::deque<ABMVector<double>>(1, y0);
    }

    ~ABMSolver()
    {
        delete x;
        delete y;
    }

    ABMVector<double> Step(const ABMVector<double> x)
    {
        if (stepCount < 3) // Initialize with RK4 method
        {
            this->x->push_back(x);
            ABMVector res = RK4(stepCount*h, this->y->back(), this->x->back());
            
            this->y->push_back(res);
            stepCount++;
            return res;
        }
        else // Use Adams-Bashforth-Adams method
        {
            this->x->pop_front();
            this->x->push_back(x);
            ABMVector res = ABM45(stepCount*h, *this->y, *this->x);
            
            this->y->pop_front();
            this->y->push_back(res);

            stepCount++;
            return res;
        }
    }

    void ChangeStepSize(double stepSize)
    {
        h = stepSize;
    }

private:

    uint32_t nVariables {1};
    double h = 1.0E-6;
    const ABMVector <std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)>> eqs;
    std::deque<ABMVector<double>> *x, *y;
    uint64_t stepCount {0};

    ABMVector<double> RK4(double t, const ABMVector<double>& y0, const ABMVector<double>& x) const
    {
        ABMVector<double> result(y0);
        double* k = new double[4];

        uint16_t i{};       
        for (i = 0; i < nVariables; i++)
        {           
            k[0] = eqs[i](t,         y0             , x);
            k[1] = eqs[i](t + 0.5*h, y0 + 0.5*k[0]*h, x);
            k[2] = eqs[i](t + 0.5*h, y0 + 0.5*k[1]*h, x);
            k[3] = eqs[i](t +     h, y0 +     k[2]*h, x);

            result[i] += SIXTH * h*(k[0] + 2*k[1] + 2*k[2] + k[3]);           
        }

        delete[] k;
        return result;
    }

    ABMVector<double> ABM45(double t, const std::deque<ABMVector<double>>& y, const std::deque<ABMVector<double>>& x) const
    {
        ABMVector<double> ABResult(y[3]);
        ABMVector<double> ABMResult(y[3]);
        double* k = new double[5];

        uint16_t i{};
        for (i = 0; i < nVariables; i++)
        {
            k[0] = AB_CONST1 * eqs[i](t      , y[3-0], x[3-0]); // t
            k[1] = AB_CONST2 * eqs[i](t -   h, y[3-1], x[3-1]); // t-1
            k[2] = AB_CONST3 * eqs[i](t - 2*h, y[3-2], x[3-2]); // t-2
            k[3] = AB_CONST4 * eqs[i](t - 3*h, y[3-3], x[3-3]); // t-3

            ABResult[i] += h*(k[0] + k[1] + k[2] + k[3]);
        }

        for (i = 0; i < nVariables; i++)
        {
            k[0] = AM_CONST1 * eqs[i](t +   h, ABResult[i], x[3-0]); // t+1
            k[1] = AM_CONST2 * eqs[i](t      ,      y[3-0], x[3-0]); // t
            k[2] = AM_CONST3 * eqs[i](t -   h,      y[3-1], x[3-1]); // t-1
            k[3] = AM_CONST4 * eqs[i](t - 2*h,      y[3-2], x[3-2]); // t-2
            k[4] = AM_CONST5 * eqs[i](t - 3*h,      y[3-3], x[3-3]); // t-3

            ABMResult[i] += h*(k[0] + k[1] + k[2] + k[3] + k[4]);
        }

        delete[] k;
        return ABMResult;
    }

};