#ifndef PROBABILITY
#define PROBABILITY

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

#include <iostream>

typedef unsigned char uint8;

constexpr uint8 numberOfFlows = 2;
constexpr size_t maxCarsInFlows[numberOfFlows] = { 10, 10 }; // ������� ��������
constexpr uint8 nomberOfStates = 2 * numberOfFlows; // ���������� ��������� �������������� ����������
constexpr size_t statesCount = (static_cast<size_t>(maxCarsInFlows[0]) + 1) * (maxCarsInFlows[1] + 1);
constexpr uint8 modeCount = 3;

enum mode { // ������ ������ �������������� ����������
    first,
    second,
    third
};

class State {
public:
    State();
    uint8 getP1();
    uint8 getQ1();
    uint8 getP2();
    uint8 getQ2();
    virtual ~State();
protected:
    State(const uint8 _p1, const uint8 _q1, const uint8 _p2, const uint8 _q2);
    uint8 p1, q1;
    uint8 p2, q2;
};

class FirstS : public State {
public:
    FirstS();
};

class SecondS : public State {
public:
    SecondS();
};

class ThirdS : public State {
public:
    ThirdS();
};

void LinearSolver(const double* Acopy, const double* b, double* x, const size_t size = statesCount);

void matrixMultiplication(const double* mat1, const double* mat2, double* const res, const size_t size);

template <class T>
static void printMatrix(const T* mat, const size_t size);

void printVector(const double* vec, const size_t size);

namespace probability {
    const double alpha[numberOfFlows] = { 0.1, 0.05 }; // ����������� �������� ������� � ����� i
    const double c[numberOfFlows] = { 1.0, 1.0 }; // �����������, ��� � ����� �������� 1 ������
    const double beta[numberOfFlows] = { 0.6, 0.65 }; // ����������� ���������� ������ �������

    double getP(const size_t curState, const size_t nextState, const uint8 flow);

    double getQ(const size_t curState, const size_t nextState, const uint8 flow);

    double getG(const size_t curState, const size_t nextState, const uint8 flow);

    double getH(const size_t curState, const size_t nextState, const uint8 flow);

    class Flow {
    public:
        Flow(uint8 flow);
        ~Flow();
        double* Powermatrix(double* const mat, const uint8 power) const;
        void PFormation();
        void QFormation();
        void GFormation();
        void HFormation();
        static void changeState(const size_t md);
        void printP();
        void printQ();
        void printG();
        void printH();
        void getMarginalProbability(const double* marginalVector, double* const res) const;
        double expectedValue(const double* marginalProbabilities) const;
        size_t getDeltaCount(const size_t number) const;
        void VectorH(double* currentH) const;
        double variance(const double* marginalProbabilities) const;
        friend double covariance(const Flow& f1, const Flow& f2, const double* marginalVector);
        friend void transitionMatrix(const Flow& f1, const Flow& f2, double* const res);
    private:
        double* P;
        double* Q;
        double* G;
        double* H;
        uint8 flow;
        static State* st;
    };

    double covariance(const Flow& f1, const Flow& f2, const double* marginalVector);

    void transitionMatrix(const Flow& f1, const Flow& f2, double* const res);

    void ExpectedValueOfRequestTime(const double* vectorH1, const double* vectorH2, double* Z);

    void VectorOfMarginalProbability(double* const trm, double* const res);

    void WeightDetermination(double* P, const double* Z, double* u);

    bool SolutionImprovement(double** P, double** Z, double* u, uint8* d);

    void HowardAlgorithm(uint8* const modes);
}

template <class T>
void printMatrix(const T* mat, const size_t size) {
    std::cout << "\t";
    for (size_t i = 0; i < size; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i < size; ++i) {
        std::cout << i << "\t";
        for (size_t j = 0; j < size; ++j)
            std::cout << mat[i * size + j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

template <>
void printMatrix<uint8>(const uint8* mat, const size_t size) {
    //std::cout << "  ";
    //for (size_t i = 0; i < size; ++i)
    //    std::cout << i << "  ";
    std::cout << std::endl;
    for (size_t i = 0; i < size; ++i) {
        //std::cout << i << "  ";
        for (size_t j = 0; j < size; ++j)
            std::cout << static_cast<int>(mat[i * size + j]) << "  ";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

#endif // !PROBABILITY