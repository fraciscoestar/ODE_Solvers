#include <iostream>
#include <functional>
#include <fstream>
#include <thread>
#include <chrono>
#include "ABMSolver.hpp"

#define Cc_inv 9.442870633
#define Rc_inv 0.001
#define alpha 0.0386
#define R0 230.0
#define Tclk (1./2.4E6)
#define Tamb 21.0

// static constexpr double Rxx[6][6] = {{1.0E3, 1.0E3, 1.5E3, 1.5E3, 1.5E3, 1.5E3},
//                                              {1.0E3, 1.0E3, 1.5E3, 1.5E3, 1.5E3, 1.5E3},
//                                              {3.0E3, 3.0E3, 1.0E3, 3.5E3, 9.0E3, 3.5E3},
//                                              {3.0E3, 3.0E3, 3.5E3, 1.0E3, 3.5E3, 9.0E3},
//                                              {3.0E3, 3.0E3, 9.0E3, 3.5E3, 1.0E3, 3.5E3},
//                                              {3.0E3, 3.0E3, 3.5E3, 9.0E3, 3.5E3, 1.0E3}};

double T1(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a1 = Rc_inv + 1./1.0E3 + 1./1.5E3 + 1./1.5E3 + 1./1.5E3 + 1./1.5E3;
    const double ih_sqr_R0 = x[0]*x[0]*R0;
    const double b1 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a1 + alpha*ih_sqr_R0)*T[0] + (1./1.0E3)*T[1] + (1./1.5E3)*T[2] + (1./1.5E3)*T[3] + (1./1.5E3)*T[4] + (1./1.5E3)*T[5] + b1);
}

double T2(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a2 = Rc_inv + 1./1.0E3 + 1./1.5E3 + 1./1.5E3 + 1./1.5E3 + 1./1.5E3;
    const double ih_sqr_R0 = x[1]*x[1]*R0;
    const double b2 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a2 + alpha*ih_sqr_R0)*T[1] + (1./1.0E3)*T[0] + (1./1.5E3)*T[2] + (1./1.5E3)*T[3] + (1./1.5E3)*T[4] + (1./1.5E3)*T[5] + b2);
}

double T3(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a3 = Rc_inv + 1./1.5E3 + 1./1.5E3 + 1./2.0E3 + 1.0/2.0E3 + 1./2.5E3;
    const double ih_sqr_R0 = x[2]*x[2]*R0;
    const double b3 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a3 + alpha*ih_sqr_R0)*T[2] + (1./1.5E3)*T[0] + (1./1.5E3)*T[1] + (1./2.0E3)*T[3] + (1.0/2.5E3)*T[4] + (1./2.0E3)*T[5] + b3);
}

double T4(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a4 = Rc_inv + 1./1.5E3 + 1./1.5E3 + 1./2.0E3 + 1.0/2.0E3 + 1./2.5E3;
    const double ih_sqr_R0 = x[3]*x[3]*R0;
    const double b4 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a4 + alpha*ih_sqr_R0)*T[3] + (1./1.5E3)*T[0] + (1./1.5E3)*T[1] + (1./2.0E3)*T[2] + (1./2.0E3)*T[4] + (1./2.5E3)*T[5] + b4);
}

double T5(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a5 = Rc_inv + 1./1.5E3 + 1./1.5E3 + 1./2.0E3 + 1.0/2.0E3 + 1./2.5E3;
    const double ih_sqr_R0 = x[4]*x[4]*R0;
    const double b5 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a5 + alpha*ih_sqr_R0)*T[4] + (1./1.5E3)*T[0] + (1./1.5E3)*T[1] + (1./2.5E3)*T[2] + (1./2.0E3)*T[3] + (1./2.0E3)*T[5] + b5);
}

double T6(double t, const ABMVector<double>& T, const ABMVector<double>& x)
{
    constexpr double a6 = Rc_inv + 1./1.5E3 + 1./1.5E3 + 1./2.0E3 + 1.0/2.0E3 + 1./2.5E3;
    const double ih_sqr_R0 = x[5]*x[5]*R0;
    const double b6 = ih_sqr_R0 + Rc_inv*Tamb;

    return Cc_inv*((-a6 + alpha*ih_sqr_R0)*T[5] + (1./1.5E3)*T[0] + (1./1.5E3)*T[1] + (1./2.0E3)*T[2] + (1./2.5E3)*T[3] + (1./2.0E3)*T[4] + b6);
}

int main (int argc, char* argv[])
{
    std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)> eq1 = [](double t, const ABMVector<double>& y, const ABMVector<double>& x) -> double
    {
        const double ih_sqr_R0 = x[0]*x[0]*R0;
        return Cc_inv * ((ih_sqr_R0*alpha - Rc_inv)*y[0] + ih_sqr_R0);
    };

    ABMSolver solver(Tclk, ABMVector<std::function<double(double, const ABMVector<double>&, const ABMVector<double>&)>>({T1, T2, T3, T4, T5, T6}),
                     ABMVector<double>({Tamb, Tamb, Tamb, Tamb, Tamb, Tamb}), ABMVector<double>({9.0E-3, 9.0E-3, 9.0E-3, 9.0E-3, 9.0E-3, 9.0E-3}));

    std::ofstream csv(std::string("output.csv"), std::ios_base::out);
    csv << "t,T\n0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << std::endl;

    // ABMVector<double> T(1, 0.0);
    uint64_t i{};
    uint64_t niter = (uint64_t)(0.1 / Tclk);

    auto startTime = std::chrono::high_resolution_clock::now();
    for (i = 1; i <= niter; i++)
    {
        double current{};
        if (i%4000 < 1500)
            current = 9.0E-3;
        else if (i%4000 < 3200)
            current = 0.0;
        else
            current = 5.E-3;

        ABMVector<double> sol = solver.Step(ABMVector<double>({current, current, current, current, current, current}));
        csv << i*Tclk << ", " << sol[0] << ", " << sol[1] << ", " << sol[2] << ", " << sol[3] << ", " << sol[4] << ", " << sol[5] << "\n";
    }

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - startTime);
    std::cout << "Finished. Elapsed time: " << duration.count() << "ms" << std::endl;
    csv.close();

    return 0;
}