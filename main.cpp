#include <iostream>
#include <functional>
#include <fstream>
#include "ABMSolver.hpp"

#define Cc_inv 1./0.2152
#define Rc_inv 1000
#define alpha 0.0386
#define R0 230
#define Tclk (1./2.4E6)

int main (int argc, char* argv[])
{
    std::function<double(double t, const ABMVector<double>& y, const ABMVector<double>& x)> eq1 = [](double t, const ABMVector<double>& y, const ABMVector<double>& x) -> double
    {
        const double ih_sqr_R0 = x[0]*x[0]*R0;
        return Cc_inv * ((ih_sqr_R0*alpha - Rc_inv)*y[0] + ih_sqr_R0);
    };

    std::ofstream csv(std::string("output.csv"), std::ios_base::app | std::ios_base::out);
    csv << "t,T\n";

    ABMSolver solver(Tclk, ABMVector({eq1}));

    ABMVector<double> T(1, 0.0);
    uint32_t i{};
    uint32_t niter = (uint32_t)(1. / Tclk);

    std::cout << "Niter: " << niter << std::endl;
    for (i = 0; i < niter; i++)
    {
        T[0] = solver.RK4(i*Tclk, T, ABMVector<double>({9.0E-3}))[0];
        csv << i*Tclk << "," << T[0] << "\n";
    }

    std::cout << "Finished." << std::endl;
    csv.close();

    return 0;
}