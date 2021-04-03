#include "probability.h"
#include <iostream>
//#include <boost/numeric/ublas/matrix.hpp>

int main(int argc, char** argv) {
    std::cout << "----------------First mode--------------" << std::endl;
    // первый режим обслуживания
    probability::Flow f11(0);
    probability::Flow f12(1);
    probability::Flow::changeState(first);

    f11.PFormation();
    f11.QFormation();
    f11.GFormation();
    f11.HFormation();
    f12.PFormation();
    f12.QFormation();
    f12.GFormation();
    f12.HFormation();

    /*std::cout << "flow 1, matrix G:" << std::endl;
    f11.printG();
    std::cout << "flow 1, matrix H:" << std::endl;
    f11.printH();
    std::cout << "flow 2, matrix G:" << std::endl;
    f12.printG();
    std::cout << "flow 2, matrix H:" << std::endl;
    f12.printH();*/

    double* trm = new double[statesCount * statesCount]{ 0.0 };
    probability::transitionMatrix(f11, f12, trm);
    double* mv = new double[statesCount] { 0.0 };
    probability::VectorOfMarginalProbability(trm, mv);
    double* mv1 = new double[maxCarsInFlows[0] + 1];
    double* mv2 = new double[maxCarsInFlows[1] + 1];
    f11.getMarginalProbability(mv, mv1);
    f12.getMarginalProbability(mv, mv2);
    double* vectorH11 = new double[maxCarsInFlows[0] + 1] {};
    double* vectorH12 = new double[maxCarsInFlows[1] + 1] {};
    double* Z1 = new double[statesCount];
    f11.VectorH(vectorH11);
    printVector(vectorH11, maxCarsInFlows[0] + 1);
    f12.VectorH(vectorH12);
    printVector(vectorH12, maxCarsInFlows[0] + 1);
    probability::ExpectedValueOfRequestTime(vectorH11, vectorH12, Z1);
    printMatrix(Z1, size_t(maxCarsInFlows[0] + 1));
    std::cout << "Expected Value of the number of customers in 1-st flow " << f11.expectedValue(mv1) << std::endl;
    std::cout << "Expected Value of the number of customers in 2-nd flow " << f12.expectedValue(mv2) << std::endl;
    std::cout << "Variance of the number of customers in 1-st flow " << f11.variance(mv1) << std::endl;
    std::cout << "Variance of the number of customers in 2-nd flow " << f12.variance(mv2) << std::endl;
    std::cout << "Covariance of the number of customers " << probability::covariance(f11, f12, mv) << std::endl;
    delete[] mv2;
    delete[] mv1;
    delete[] mv;
    delete[] trm;
    delete[] vectorH11;
    delete[] vectorH12;
    delete[] Z1;

    std::cout << std::endl << std::endl << "----------------Second mode--------------" << std::endl;
    // второй режим обслуживания
    probability::Flow f21(0);
    probability::Flow f22(1);
    probability::Flow::changeState(second);

    f21.PFormation();
    f21.QFormation();
    f21.GFormation();
    f21.HFormation();
    f22.PFormation();
    f22.QFormation();
    f22.GFormation();
    f22.HFormation();

    trm = new double[statesCount * statesCount]{ 0.0 };
    probability::transitionMatrix(f21, f22, trm);
    mv = new double[statesCount] { 0.0 };
    probability::VectorOfMarginalProbability(trm, mv);;
    mv1 = new double[maxCarsInFlows[0] + 1];
    mv2 = new double[maxCarsInFlows[1] + 1];
    f21.getMarginalProbability(mv, mv1);
    f22.getMarginalProbability(mv, mv2);
    double* vectorH21 = new double[maxCarsInFlows[0] + 1] {};
    double* vectorH22 = new double[maxCarsInFlows[1] + 1] {};
    double* Z2 = new double[statesCount];
    f21.VectorH(vectorH21);
    printVector(vectorH21, maxCarsInFlows[0] + 1);
    f22.VectorH(vectorH22);
    printVector(vectorH22, maxCarsInFlows[1] + 1);
    probability::ExpectedValueOfRequestTime(vectorH21, vectorH22, Z2);
    printMatrix(Z2, maxCarsInFlows[0] + 1);
    std::cout << "Expected Value of the number of customers in 1-st flow " << f21.expectedValue(mv1) << std::endl;
    std::cout << "Expected Value of the number of customers in 2-nd flow " << f22.expectedValue(mv2) << std::endl;
    std::cout << "Variance of the number of customers in 1-st flow " << f21.variance(mv1) << std::endl;
    std::cout << "Variance of the number of customers in 2-nd flow " << f22.variance(mv2) << std::endl;
    std::cout << "Covariance of the number of customers " << probability::covariance(f21, f22, mv) << std::endl;
    delete[] mv2;
    delete[] mv1;
    delete[] mv;
    delete[] trm;
    delete[] vectorH21;
    delete[] vectorH22;
    delete[] Z2;

    return 0;
}