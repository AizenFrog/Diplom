#include "probability.h"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <assert.h>

#include <cmath>
#include <memory>
#include <iostream>

// ----------CUDA variables----------

cusolverDnHandle_t cusolverH = nullptr;
cublasHandle_t cublasH       = nullptr;

cublasStatus_t cublas_status     = CUBLAS_STATUS_SUCCESS;
cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

cudaError_t cudaStat1 = cudaSuccess;
cudaError_t cudaStat2 = cudaSuccess;
cudaError_t cudaStat3 = cudaSuccess;
cudaError_t cudaStat4 = cudaSuccess;

double* d_A     = nullptr;
double* d_tau   = nullptr;
double* d_B     = nullptr;
int* devInfo    = nullptr;
double* d_work  = nullptr;
int lwork_geqrf = 0;
int lwork_ormqr = 0;
int lwork       = 0;
int info_gpu    = 0;

// ----------------------------------

double probability::getP(const size_t curState, const size_t nextState, const uint8 flow) {
    if (curState == nextState && curState == 0)
        return (1 - alpha[flow]) + alpha[flow] * beta[flow] * c[flow];
    else if (curState == nextState && curState == maxCarsInFlows[flow])
        return (1 - beta[flow]) * (1 - alpha[flow]) + alpha[flow] * beta[flow] * c[flow] + alpha[flow] * beta[flow] * (1 - c[flow])
            + alpha[flow] * (1 - beta[flow]) * c[flow] + alpha[flow] * (1 - beta[flow]) * (1 - c[flow]);
    else if (curState == nextState)
        return (1 - alpha[flow]) * (1 - beta[flow]) + alpha[flow] * beta[flow] * c[flow];
    else if (curState == nextState - 1 && nextState == maxCarsInFlows[flow])
        return alpha[flow] * (1 - beta[flow]) * c[flow] + alpha[flow] * beta[flow] * (1 - c[flow])
            + alpha[flow] * (1 - beta[flow]) * (1 - c[flow]);
    else if (curState == nextState - 1)
        return alpha[flow] * (1 - beta[flow]) * c[flow] + alpha[flow] * beta[flow] * (1 - c[flow]);
    else if (curState == nextState - 2)
        return alpha[flow] * (1 - beta[flow]) * (1 - c[flow]);
    else if (curState == nextState + 1)
        return (1 - alpha[flow]) * beta[flow];
    return 0.0;
}

double probability::getQ(const size_t curState, const size_t nextState, const uint8 flow) {
    if (curState == nextState && curState == maxCarsInFlows[flow])
        return 1;
    else if (curState == nextState)
        return 1 - alpha[flow];
    else if (curState == nextState - 1 && nextState == maxCarsInFlows[flow])
        return alpha[flow] * c[flow] + alpha[flow] * (1 - c[flow]);
    else if (curState == nextState - 1)
        return alpha[flow] * c[flow];
    else if (curState == nextState - 2)
        return alpha[flow] * (1 - c[flow]);
    return 0.0;
}

double probability::getG(const size_t curState, const size_t nextState, const uint8 flow) {
    //if (curState == nextState && curState == maxCarsInFlows[flow])
    //    return 1;
    if (curState == nextState && curState == 0)
        return 0.0;
    else if (curState == nextState)
        return curState;
    else if (curState == nextState - 1)
        return static_cast<double>(curState) + 0.5;
    else if (curState == nextState - 2)
        return curState + 1;
    return 0.0;
}

double probability::getH(const size_t curState, const size_t nextState, const uint8 flow) {
    if (curState == nextState && curState == 0)
        return (alpha[flow] * c[flow] * beta[flow]) / (4 * (1 - alpha[flow] + alpha[flow] * c[flow] * beta[flow]));
    else if (curState == nextState)
        return curState;
    else if (curState == nextState - 1)
        return static_cast<double>(curState) + 0.5;
    else if (curState == nextState - 2)
        return curState + 1;
    else if (curState == nextState + 1)
        return static_cast<double>(curState) - 0.5;
    return 0.0;
}

State* probability::Flow::st = new ThirdS;

probability::Flow::Flow(uint8 flow) : P(new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)]),
    Q(new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)]), flow(flow),
    G(new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)]),
    H(new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)]) {}

probability::Flow::~Flow() {
    delete[] P;
    delete[] Q;
    delete[] G;
    delete[] H;
}

double* probability::Flow::Powermatrix(double* const mat, const uint8 power) const {
    uint8 locPower = 1;
    double* res = new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)]{ 0.0 };
    double* tmp = new double[(static_cast<size_t>(maxCarsInFlows[flow]) + 1) * (maxCarsInFlows[flow] + 1)];
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j)
            tmp[i * (maxCarsInFlows[flow] + 1) + j] = mat[i * (maxCarsInFlows[flow] + 1) + j];
    while (locPower != power) {
        if (locPower != 1) {
            double* temple = tmp;
            tmp = res;
            res = temple;
        }
        for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
            for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j) {
                res[i * (maxCarsInFlows[flow] + 1) + j] = 0.0;
                for (size_t k = 0; k <= maxCarsInFlows[flow]; ++k) {
                    res[i * (maxCarsInFlows[flow] + 1) + j] += tmp[i * (maxCarsInFlows[flow] + 1) + k]
                        * mat[k * (maxCarsInFlows[flow] + 1) + j];
                }
            }
        ++locPower;
    }
    if (power > 1) {
        delete[] tmp;
        return res;
    } else {
        delete[] res;
        return tmp;
    }
}

void probability::Flow::PFormation() {
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j) {
            P[(maxCarsInFlows[flow] + 1) * i + j] = getP(i, j, flow);
        }
    }
}

void probability::Flow::QFormation() {
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j) {
            Q[(maxCarsInFlows[flow] + 1) * i + j] = getQ(i, j, flow);
        }
    }
}

void probability::Flow::GFormation() {
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j) {
            G[(maxCarsInFlows[flow] + 1) * i + j] = getG(i, j, flow);
        }
    }
}

void probability::Flow::HFormation() {
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j) {
            H[(maxCarsInFlows[flow] + 1) * i + j] = getH(i, j, flow);
        }
    }
}

void probability::Flow::changeState(mode md) {
    delete st;
    switch (md) {
    case first:
        st = new FirstS;
        break;
    case second:
        st = new SecondS;
        break;
    case third:
        st = new ThirdS;
        break;
    default:
        st = new ThirdS;
        break;
    }
}

void probability::Flow::printP() {
    std::cout << "\t";
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        std::cout << i << "\t";
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j)
            std::cout << P[i * (maxCarsInFlows[flow] + 1) + j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void probability::Flow::printQ() {
    std::cout << "\t";
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        std::cout << i << "\t";
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j)
            std::cout << Q[i * (maxCarsInFlows[flow] + 1) + j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void probability::Flow::printG() {
    std::cout << "\t";
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        std::cout << i << "\t";
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j)
            std::cout << G[i * (maxCarsInFlows[flow] + 1) + j] << "\t";
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void probability::Flow::printH() {
    std::cout << "\t";
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
        std::cout << i << "\t";
        for (size_t j = 0; j <= maxCarsInFlows[flow]; ++j)
            std::cout << H[i * (maxCarsInFlows[flow] + 1) + j] << "\t";
        std::cout << std::endl;
    }
}

void probability::Flow::getMarginalProbability(const double* marginalVector, double* const res) const {
    if (flow == 0) {
        for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
            res[i] = 0.0;
            for (size_t j = i; j < statesCount; j += maxCarsInFlows[flow] + 1)
                res[i] += marginalVector[j];
        }
    } else {
        for (size_t i = 0; i <= maxCarsInFlows[flow]; ++i) {
            res[i] = 0.0;
            for (size_t j = i * (maxCarsInFlows[0] + 1); j < (i + 1) * (maxCarsInFlows[0] + 1); ++j)
                res[i] += marginalVector[j];
        }
    }
}

double probability::Flow::expectedValue(const double* marginalProbabilities) const {
    double res = 0.0;
    for (size_t i = 0; i < maxCarsInFlows[flow] + 1; ++i)
        res += marginalProbabilities[i] * i;
    return res;
}

void probability::Flow::VectorH(double* currentH) const {
    double* tempH = new double[maxCarsInFlows[flow] + 1] {};
    size_t s1 = flow == 0 ? Flow::st->getQ1() + Flow::st->getP2() + Flow::st->getQ2() : Flow::st->getQ2(),
           s2 = flow == 0 ? Flow::st->getP1() : Flow::st->getP2(),
           s3 = flow == 0 ? 0 : Flow::st->getP1() + Flow::st->getQ1();
    for (size_t i = 0; i < s1; ++i) {
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l)
            for (size_t k = 0; k <= maxCarsInFlows[flow]; ++k)
                tempH[l] += Q[(maxCarsInFlows[flow] + 1) * l + k] * (G[(maxCarsInFlows[flow] + 1) * l + k] + currentH[k]);
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l) {
            currentH[l] = tempH[l];
            tempH[l] = 0.0;
        }
    }
    //printVector(currentH, maxCarsInFlows[0] + 1);
    for (size_t i = 0; i < s2; ++i) {
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l)
            for (size_t k = 0; k <= (maxCarsInFlows[flow]); ++k)
                tempH[l] += P[(maxCarsInFlows[flow] + 1) * l + k] * (H[(maxCarsInFlows[flow] + 1) * l + k] + currentH[k]);
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l) {
            currentH[l] = tempH[l];
            tempH[l] = 0.0;
        }
    }
    //printVector(currentH, maxCarsInFlows[0] + 1);
    for (size_t i = 0; i < s3; ++i) {
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l)
            for (size_t k = 0; k <= maxCarsInFlows[flow]; ++k)
                tempH[l] += Q[(maxCarsInFlows[flow] + 1) * l + k] * (G[(maxCarsInFlows[flow] + 1) * l + k] + currentH[k]);
        for (size_t l = 0; l <= maxCarsInFlows[flow]; ++l) {
            currentH[l] = tempH[l];
            tempH[l] = 0.0;
        }
    }
    //printVector(currentH, maxCarsInFlows[0] + 1);
    delete[] tempH;
}

double probability::Flow::variance(const double* marginalProbabilities) const {
    double ev = expectedValue(marginalProbabilities);
    double res = 0.0;
    for (size_t i = 0; i < maxCarsInFlows[flow] + 1; ++i)
        res += marginalProbabilities[i] * ((i - ev) * (i - ev));
    return res;
}

double probability::covariance(const Flow& f1, const Flow& f2, const double* marginalVector) {
    double res = 0.0;
    double* mp1 = new double[maxCarsInFlows[f1.flow] + 1];
    double* mp2 = new double[maxCarsInFlows[f2.flow] + 1];
    f1.getMarginalProbability(marginalVector, mp1);
    f2.getMarginalProbability(marginalVector, mp2);
    double ev1 = f1.expectedValue(mp1);
    double ev2 = f2.expectedValue(mp2);
    double* value = new double[statesCount];
    for (size_t i = 0; i < maxCarsInFlows[f1.flow] + 1; ++i)
        for (size_t j = 0; j < maxCarsInFlows[f2.flow] + 1; ++j)
            value[i * (maxCarsInFlows[f2.flow] + 1) + j] = (i - ev1) * (j - ev2);
    for (size_t i = 0; i < statesCount; ++i)
        res += marginalVector[i] * value[i];
    delete[] value;
    delete[] mp2;
    delete[] mp1;
    return res;
}

void probability::transitionMatrix(const Flow& f1, const Flow& f2, double* res) {
    std::unique_ptr<double[]> P1(f1.Powermatrix(f1.P, Flow::st->getP1()));
    std::unique_ptr<double[]> Q1(f1.Powermatrix(f1.Q, Flow::st->getQ1() + Flow::st->getP2() + Flow::st->getQ2()));
    std::unique_ptr<double[]> P2(f2.Powermatrix(f2.P, Flow::st->getP2()));
    std::unique_ptr<double[]> Q2(f2.Powermatrix(f2.Q, Flow::st->getP1() + Flow::st->getQ1()));
    double* res1 = new double[(static_cast<size_t>(maxCarsInFlows[0]) + 1) * (maxCarsInFlows[0] + 1)]{ 0.0 };
    matrixMultiplication(P1.get(), Q1.get(), res1, maxCarsInFlows[0] + 1);
    double* tmp = new double[(static_cast<size_t>(maxCarsInFlows[1]) + 1) * (maxCarsInFlows[1] + 1)]{ 0.0 };
    matrixMultiplication(Q2.get(), P2.get(), tmp, maxCarsInFlows[1] + 1);
    double* res2 = new double[(static_cast<size_t>(maxCarsInFlows[1]) + 1) * (maxCarsInFlows[1] + 1)]{ 0.0 };
    Q2.reset(f2.Powermatrix(f2.Q, Flow::st->getQ2()));
    matrixMultiplication(tmp, Q2.get(), res2, maxCarsInFlows[1] + 1);
    delete[] tmp;

    for (size_t i = 0; i <= maxCarsInFlows[0]; ++i) {
        for (size_t j = 0; j <= maxCarsInFlows[1]; ++j) {
            for (size_t k = 0; k <= maxCarsInFlows[0]; ++k) {
                for (size_t g = 0; g <= maxCarsInFlows[1]; ++g) {
                    res[(i * (maxCarsInFlows[1] + 1) + j) * statesCount + (g * (maxCarsInFlows[0] + 1) + k)] =
                        res1[i * (maxCarsInFlows[0] + 1) + k] * res2[j * (maxCarsInFlows[1] + 1) + g];
                }
            }
        }
    }

    delete[] res1;
    delete[] res2;
}

double isZero(const double number) {
    if (std::abs(number) < 0.0000000001) {
        return 0.0;
    }
    else {
        return number;
    }
}

void probability::VectorOfMarginalProbability(double* const trm, double* const res) {
    double* coef = new double[statesCount - 1] { 0.0 };
    double* transposedMatrix = new double[statesCount * statesCount];
    for (size_t i = 0; i < statesCount; ++i)
        for (size_t j = 0; j < statesCount; ++j)
            transposedMatrix[j * statesCount + i] = trm[i * statesCount + j];
    for (int k = 0; k < statesCount - 1; ++k) {
        transposedMatrix[k * statesCount + k] -= 1.0;
        double leaderElem = transposedMatrix[k * statesCount + k];
        for (int rows = k + 1; rows < statesCount - 1; ++rows) {
            double leaderRow = transposedMatrix[rows * statesCount + k];
            for (int cols = k; cols < statesCount; ++cols) {
                transposedMatrix[rows * statesCount + cols] = isZero(transposedMatrix[rows * statesCount + cols] +
                    transposedMatrix[k * statesCount + cols] * std::abs(leaderRow / leaderElem));
            }
        }
    }
    coef[statesCount - 2] = isZero(-(transposedMatrix[statesCount * (statesCount - 2) + (statesCount - 1)] /
        transposedMatrix[statesCount * (statesCount - 2) + (statesCount - 2)]));
    for (long long i = statesCount - 3; i >= 0; --i) {
        double result = 0.0;
        for (size_t j = statesCount - 2; j >= i + 1; --j)
            result += isZero(coef[j] * transposedMatrix[statesCount * i + j]);
        result += transposedMatrix[statesCount * i + (statesCount - 1)];
        coef[i] = isZero(-(result / transposedMatrix[i * statesCount + i]));
    }
    for (size_t i = 0; i < statesCount - 1; ++i)
        res[statesCount - 1] += coef[i];
    res[statesCount - 1] += 1.0;
    res[statesCount - 1] = 1 / res[statesCount - 1];
    double check = res[statesCount - 1];
    for (size_t i = 0; i < statesCount - 1; ++i) {
        res[i] = coef[i] * res[statesCount - 1];
        check += res[i];
    }
    // std::cout << check << std::endl;
    delete[] coef;
    delete[] transposedMatrix;
}

void probability::WeightDetermination(double* P, const double* Z, double* u) {
    double* Q = new double[statesCount];
    for (size_t i = 0; i < statesCount; ++i) {
        P[i * statesCount + i] -= 1.0;
        P[(i + 1) * statesCount - 1] = -1.0;
        Q[i] = -Z[i];
    }

    LinearSolver(P, Q, u);
    delete[] Q;
}

bool probability::SolutionImprovement(double** P, double** Z, double* u, uint8* d) {
    double min = 0.0;
    bool change = false;
    for (size_t i = 0; i < statesCount; ++i) {
        uint8 index = 0;
        for (uint8 k = 0; k < modeCount; ++k) {
            double sum = 0.0;
            for (size_t j = 0; j < statesCount; ++j)
                sum += P[k][i * statesCount + j] * u[j];
            sum += Z[k][i];
            if (sum < min) {
                min = sum;
                index = k;
            }
        }
        if (d[i] != index) {
            d[i] = index;
            change = true;
        }
    }
    return change;
}

void probability::HowardAlgorithm(uint8* const modes) {
    Flow* allFlows = static_cast<Flow*>(operator new(sizeof(Flow) * (uint32_t)numberOfFlows * modeCount));
    for (size_t i = 0; i < (uint32_t)numberOfFlows * modeCount; ++i)
        new(allFlows + i) Flow(i % 2);

    double** trm = new double* [modeCount];
    double** Z   = new double* [modeCount];

    for (size_t i = 0; i < modeCount; ++i) {
        trm[i] = new double[statesCount * statesCount];
        Z[i] = new double[statesCount];

        // первый поток при определенном режиме
        allFlows[i * 2].PFormation();
        allFlows[i * 2].QFormation();
        allFlows[i * 2].GFormation();
        allFlows[i * 2].HFormation();

        // второй поток при определенном рижиме
        allFlows[i * 2 + 1].PFormation();
        allFlows[i * 2 + 1].QFormation();
        allFlows[i * 2 + 1].GFormation();
        allFlows[i * 2 + 1].HFormation();

        transitionMatrix(allFlows[i * 2], allFlows[i * 2 + 1], trm[i]);
        double* vectorH11 = new double[maxCarsInFlows[0] + 1]{};
        double* vectorH12 = new double[maxCarsInFlows[1] + 1]{};
        allFlows[i * 2].VectorH(vectorH11);
        allFlows[i * 2 + 1].VectorH(vectorH12);
        ExpectedValueOfRequestTime(vectorH11, vectorH12, Z[i]);

        delete[] Z[i];
        delete[] trm[i];
    }

    for (size_t i = 0; i < statesCount; ++i) {
        double min = Z[0][i];
        uint8 min_index = 0;
        for (uint8 j = 1; j < modeCount; ++j)
            if (Z[j][i] < min) {
                min_index = j;
                min = Z[j][i];
            }
        
        new(modes + i) uint8(min_index);
    }

    double* u = new double[statesCount];

    double* P = new double[statesCount * statesCount];
    double* P_cur = new double[statesCount * statesCount];
    double* q = new double[statesCount];

    for (size_t i = 0; i < statesCount; ++i) {
        q = Z[modes[i]];
        for (size_t j = 0; j < statesCount; ++j)
            P[i * statesCount + j] = P_cur[i * statesCount + j] = trm[modes[i]][i * statesCount + j];
    }

    // create handle
    cusolver_status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    cublas_status = cublasCreate(&cublasH);
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);

    // copy memory to gpu
    cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * statesCount * statesCount);
    cudaStat2 = cudaMalloc((void**)&d_tau, sizeof(double) * statesCount);
    cudaStat3 = cudaMalloc((void**)&d_B, sizeof(double) * statesCount);
    cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);
    assert(cudaSuccess == cudaStat3);
    assert(cudaSuccess == cudaStat4);

    uint8* d = new uint8[statesCount];
    bool iter = true;
    do {
        WeightDetermination(P_cur, q, u);
        iter = SolutionImprovement(trm, Z, u, d);
    } while (iter);

    if (d_A) cudaFree(d_A);
    if (d_tau) cudaFree(d_tau);
    if (d_B) cudaFree(d_B);
    if (devInfo) cudaFree(devInfo);
    if (d_work) cudaFree(d_work);
    if (cublasH) cublasDestroy(cublasH);
    if (cusolverH) cusolverDnDestroy(cusolverH);
    cudaDeviceReset();

    delete[] d;
    delete[] P;
    delete[] P_cur;
    delete[] u;
    delete[] Z;
    delete[] trm;

    for (size_t i = 0; i < (uint32_t)numberOfFlows * modeCount; ++i)
        (allFlows + i)->~Flow();
    operator delete(allFlows);
}

void probability::ExpectedValueOfRequestTime(const double* vectorH1, const double* vectorH2, double* Z) {
    for (size_t i = 0; i <= maxCarsInFlows[0]; ++i)
        for (size_t j = 0; j <= maxCarsInFlows[1]; ++j)
            Z[i * (maxCarsInFlows[0] + 1) + j] = vectorH1[i] + vectorH2[j];
}

State::State() : p1(0), q1(0), p2(0), q2(0) {}

uint8 State::getP1() {
    return p1;
}

uint8 State::getQ1() {
    return q1;
}

uint8 State::getP2() {
    return p2;
}

uint8 State::getQ2() {
    return q2;
}

State::State(const uint8 _p1, const uint8 _q1, const uint8 _p2, const uint8 _q2)
    : p1(_p1), q1(_q1), p2(_p2), q2(_q2) {}

State::~State() {}

FirstS::FirstS() : State(3, 1, 2, 1) {}

SecondS::SecondS() : State(2, 1, 3, 1) {}

ThirdS::ThirdS() : State(2, 1, 2, 1) {}

void LinearSolver(const double* Acopy, const double* b, double* x) {
    /*int bufferSize = 0;
    int* info = nullptr;
    double* buffer = nullptr;
    double* A = nullptr;
    int* ipiv = nullptr; // pivoting sequence
    int h_info = 0;
    double time_solve;
    cusolverStatus_t status;
    cudaError_t error;

    status = cusolverDnDgetrf_bufferSize(handler, n, n, const_cast<double*>(A), n, &bufferSize);

    error = cudaMalloc(&info, sizeof(int));
    error = cudaMalloc(&buffer, sizeof(double) * bufferSize);
    error = cudaMalloc(&A, sizeof(double) * n * n);
    error = cudaMalloc(&ipiv, sizeof(int) * n);

    error = cudaMemcpy(A, Acopy, sizeof(double) * n * n, cudaMemcpyDeviceToDevice);
    error = cudaMemset(info, 0, sizeof(int));

    status = cusolverDnDgetrf(handler, n, n, A, n, buffer, ipiv, info);
    error = cudaMemcpy(&h_info, info, sizeof(int), cudaMemcpyDeviceToHost);
    if (0 != h_info) {
        fprintf(stderr, "Error: LU factorization failed\n");
        printf("INFO_VALUE = %d\n", h_info);
    }

    error = cudaMemcpy(x, b, sizeof(double) * n, cudaMemcpyDeviceToDevice);
    status = cusolverDnDgetrs(handler, CUBLAS_OP_N, n, 1, A, n, ipiv, x, n, info);
    cudaDeviceSynchronize();

    if (info) cudaFree(info);
    if (buffer) cudaFree(buffer);
    if (A) cudaFree(A);
    if (ipiv) cudaFree(ipiv);*/

    cudaStat1 = cudaMemcpy(d_A, Acopy, sizeof(double) * statesCount * statesCount, cudaMemcpyHostToDevice);
    cudaStat2 = cudaMemcpy(d_B, b, sizeof(double) * statesCount, cudaMemcpyHostToDevice);
    assert(cudaSuccess == cudaStat1);
    assert(cudaSuccess == cudaStat2);

    cusolver_status = cusolverDnDgeqrf_bufferSize(cusolverH, statesCount, statesCount, d_A, statesCount, &lwork_geqrf);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    cusolver_status = cusolverDnDormqr_bufferSize(cusolverH, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, statesCount, 1,
                                                  statesCount, d_A, statesCount, d_tau, d_B, statesCount, &lwork_ormqr);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
    lwork = (lwork_geqrf > lwork_ormqr) ? lwork_geqrf : lwork_ormqr;

    cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double) * lwork);
    assert(cudaSuccess == cudaStat1);

    cusolver_status = cusolverDnDgeqrf(cusolverH, statesCount, statesCount, d_A, statesCount, d_tau, d_work, lwork, devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);

    /* check if QR is good or not */
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    printf("after geqrf: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);

    cusolver_status = cusolverDnDormqr(cusolverH, CUBLAS_SIDE_LEFT, CUBLAS_OP_T, statesCount, 1, statesCount,
                                       d_A, statesCount, d_tau, d_B, statesCount, d_work, lwork, devInfo);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat1);
    /* check if QR is good or not */
    cudaStat1 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);

    printf("after ormqr: info_gpu = %d\n", info_gpu);
    assert(0 == info_gpu);

    const double one = 1.0;
    cublas_status = cublasDtrsm(cublasH, CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT,
                                statesCount, 1, &one, d_A, statesCount, d_B, statesCount);
    cudaStat1 = cudaDeviceSynchronize();
    assert(CUBLAS_STATUS_SUCCESS == cublas_status);
    assert(cudaSuccess == cudaStat1);

    cudaStat1 = cudaMemcpy(x, d_B, sizeof(double) * statesCount, cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat1);
}

void matrixMultiplication(const double* mat1, const double* mat2, double* const res, const size_t size) {
    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            for (size_t k = 0; k < size; ++k)
                res[i * size + j] += mat1[i * size + k] * mat2[k * size + j];
}

void printMatrix(const double* mat, const size_t size) {
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

void printVector(const double* vec, const size_t size) {
    for (size_t i = 0; i < size; ++i)
        std::cout << i << "\t";
    std::cout << std::endl;
    for (size_t i = 0; i < size; ++i)
        std::cout << vec[i] << "\t";
    std::cout << std::endl;
}
