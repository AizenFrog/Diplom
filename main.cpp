#include "probability.h"
#include <iostream>

int main(int argc, char** argv) {
    std::cout << "----------------First mode--------------" << std::endl;
    // первый режим обслуживания
    probability::Flow f11(0);
    probability::Flow f12(1);
    probability::Flow::changeState(first);

    f11.PFormation();
    f11.QFormation();
    f12.PFormation();
    f12.QFormation();

    double* trm = new double[statesCount * statesCount]{ 0.0 };
    probability::transitionMatrix(f11, f12, trm);
    double* mv = new double[statesCount] { 0.0 };
    probability::VectorOfMarginalProbability(trm, mv);;
    double* mv1 = new double[maxCarsInFlows[0] + 1];
    double* mv2 = new double[maxCarsInFlows[1] + 1];
    f11.getMarginalProbability(mv, mv1);
    f12.getMarginalProbability(mv, mv2);
    std::cout << "Expected Value of the number of customers in 1-st flow " << f11.expectedValue(mv1) << std::endl;
    std::cout << "Expected Value of the number of customers in 2-nd flow " << f12.expectedValue(mv2) << std::endl;
    std::cout << "Variance of the number of customers in 1-st flow " << f11.variance(mv1) << std::endl;
    std::cout << "Variance of the number of customers in 2-nd flow " << f12.variance(mv2) << std::endl;
    std::cout << "Covariance of the number of customers " << probability::сovariance(f11, f12, mv) << std::endl;
    delete[] mv2;
    delete[] mv1;
    delete[] mv;
    delete[] trm;

    std::cout << std::endl << std::endl << "----------------Second mode--------------" << std::endl;
    // второй режим обслуживания
    probability::Flow f21(0);
    probability::Flow f22(1);
    probability::Flow::changeState(second);

    f21.PFormation();
    f21.QFormation();
    f22.PFormation();
    f22.QFormation();

    trm = new double[statesCount * statesCount]{ 0.0 };
    probability::transitionMatrix(f21, f22, trm);
    mv = new double[statesCount] { 0.0 };
    probability::VectorOfMarginalProbability(trm, mv);;
    mv1 = new double[maxCarsInFlows[0] + 1];
    mv2 = new double[maxCarsInFlows[1] + 1];
    f21.getMarginalProbability(mv, mv1);
    f22.getMarginalProbability(mv, mv2);
    std::cout << "Expected Value of the number of customers in 1-st flow " << f21.expectedValue(mv1) << std::endl;
    std::cout << "Expected Value of the number of customers in 2-nd flow " << f22.expectedValue(mv2) << std::endl;
    std::cout << "Variance of the number of customers in 1-st flow " << f21.variance(mv1) << std::endl;
    std::cout << "Variance of the number of customers in 2-nd flow " << f22.variance(mv2) << std::endl;
    std::cout << "Covariance of the number of customers " << probability::сovariance(f21, f22, mv) << std::endl;
    delete[] mv2;
    delete[] mv1;
    delete[] mv;
    delete[] trm;
    return 0;
}