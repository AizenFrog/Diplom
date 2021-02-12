#ifndef PROBABILITY
#define PROBABILITY

typedef unsigned char uint8;

const uint8 nomberOfFlows = 2;
const size_t maxCarsInFlows[nomberOfFlows] = { 30, 30 }; // емкость очередей
const uint8 nomberOfStates = 2 * nomberOfFlows; // количество состояний обслуживающего устройства
const size_t statesCount = (static_cast<size_t>(maxCarsInFlows[0]) + 1) * (maxCarsInFlows[1] + 1);

enum mode { // режимы работы обслуживающего устройства
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

void matrixMultiplication(const double* mat1, const double* mat2, double* const res, const size_t size);

void printMatrix(const double* mat, const size_t size);

void printVector(const double* vec, const size_t size);

namespace probability {
    const double alpha[nomberOfFlows] = { 0.15, 0.2 }; // вероятность прибытия клиента в поток i
    const double c[nomberOfFlows] = { 0.6, 0.55 }; // вероятность, что в поток прибудет 1 клиент
    const double beta[nomberOfFlows] = { 0.8, 0.7 }; // вероятность завершения услуги клиента

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
        static void changeState(mode md);
        void printP();
        void printQ();
        void printG();
        void printH();
        void getMarginalProbability(const double* marginalVector, double* const res) const;
        double expectedValue(const double* marginalProbabilities) const;
        void VectorH(double* currentH) const;
        double variance(const double* marginalProbabilities) const;
        friend double сovariance(const Flow& f1, const Flow& f2, const double* marginalVector);
        friend void transitionMatrix(const Flow& f1, const Flow& f2, double* const res);
    private:
        double* P;
        double* Q;
        double* G;
        double* H;
        uint8 flow;
        static State* st;
    };

    double сovariance(const Flow& f1, const Flow& f2, const double* marginalVector);

    void transitionMatrix(const Flow& f1, const Flow& f2, double* const res);

    void ExpectedValueOfRequestTime(const double* vectorH1, const double* vectorH2, double* Z);

    void VectorOfMarginalProbability(double* const trm, double* const res);
}

#endif // !PROBABILITY