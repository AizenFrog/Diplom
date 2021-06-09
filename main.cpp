#include "probability.h"
#include <iostream>
#include <fstream>
#include <string>

class HTMLoutput {
public:
    HTMLoutput() = delete;
    HTMLoutput(const HTMLoutput&) = delete;
    HTMLoutput& operator=(const HTMLoutput&) = delete;

    HTMLoutput(const std::string& file_name) : file(file_name + ".html", std::ios_base::app) {
        if (file.is_open() == false) {
            std::cerr << "Output file is not open!" << std::endl;
            return;
        }
    }

    void InitFileStruct() {
        file << "<!DOCTYPE html>\n"
             << "<html>\n"
             << "<head>\n"
             <<    "\t<meta charset = \"UTF-8\">\n"
             <<    "\t<title>Howard output</title>\n"
             <<    "\t<style type=\"text/css\">\n"
             <<        "\t\tth {background: #ccc;}\n"
             <<        "\t\t.green {background: rgb(173, 255, 47);}\n"
             <<        "\t\t.yellow {background: rgb(255, 255, 0);}\n"
             <<        "\t\t.aqua {background: rgb(0, 255, 255);}\n"
             <<    "\t</style>\n"
             << "</head>\n"
             << "<body>\n"
             << "</body>\n"
             << "</html>\n";
    }

    void SetBreak() {
        file << "<br>" << std::endl;
    }

    void SetHeaderH3(const std::string& header) {
        file << "\t<h3><b>" << header << "</b></h3>\n";
    }

    void SetLine() {
        file << "\t<hr>\n";
    }

    void SetMatrix(const uint8* modes, const size_t size) {

        file << "<table border=\"1\">\n"
             << "\t<tr>\n"
             << "\t\t<th>2\\1</th>\n";
        for (size_t i = 0; i < size; ++i)
            file << "\t\t<th>" << i << "</th>\n";
        file << "\t</tr>\n";

        for (size_t i = 0; i < size; ++i) {
            file << "\t<tr>"
                 << "<th>" << i << "</th>";
            for (size_t j = 0; j < size; ++j)
                if (static_cast<int>(modes[i * size + j]) == 0)
                    file << "<td class=\"green\">" << static_cast<int>(modes[i * size + j]) << "</td>";
                else if (static_cast<int>(modes[i * size + j]) == 1)
                    file << "<td class=\"aqua\">" << static_cast<int>(modes[i * size + j]) << "</td>";
                else if (static_cast<int>(modes[i * size + j]) == 2)
                    file << "<td class=\"yellow\">" << static_cast<int>(modes[i * size + j]) << "</td>";
            file << "</tr>\n";
        }
        file << "</table>\n";
    }

    ~HTMLoutput() {
        file << "</body>\n"
             << "</html>";
        file.close();
    }

private:
    std::ofstream file;
};

bool outputInFile(const uint8* modes, const size_t iterationCount, const size_t numberOfExperiment,
                  const int param, const double paramValue) {
    std::ofstream file("programm_output.html", std::ios_base::app);
    if (file.is_open() == false) {
        return false;
    }
    file << "Number of experiment: " << numberOfExperiment << std::endl;
    file << "Number of iterations: " << iterationCount << std::endl;
    switch (param) {
    case 0:
        file << "alpha: " << paramValue << std::endl;
        break;
    case 1:
        file << "beta: " << paramValue << std::endl;
        break;
    case 2:
        file << "c: " << paramValue << std::endl;
        break;
    }
    file << std::endl;
    for (size_t i = 0; i < maxCarsInFlows[0] + 1; ++i) {
        for (size_t j = 0; j < maxCarsInFlows[1] + 1; ++j)
            file << static_cast<int>(modes[i * (maxCarsInFlows[0] + 1) + j]) << "  ";
        file << std::endl;
    }
    file << "--------------------------------------------" << std::endl;
    file.close();
    return true;
}

void OutputInHTML(const std::string& file_name, const uint8* modes) {
    HTMLoutput out(file_name);
    out.InitFileStruct();
    out.SetHeaderH3("lambda1 = " + std::to_string(probability::alpha[0]) + ", lambda2 = " + std::to_string(probability::alpha[1]) + '\n');
    out.SetHeaderH3("beta1 = " + std::to_string(probability::beta[0]) + ", beta2 = " + std::to_string(probability::beta[1]) + '\n');
    out.SetHeaderH3("c1 = " + std::to_string(probability::c[0]) + ", c2 = " + std::to_string(probability::c[1]) + '\n' + '\n');
    out.SetHeaderH3("Politics:\n");
    out.SetHeaderH3("1: " + std::to_string(probability::modesValue[0][0]) + "-" + std::to_string(probability::modesValue[0][1]) +
                    "-" + std::to_string(probability::modesValue[0][2]) + "-" + std::to_string(probability::modesValue[0][3]) + "\n");
    out.SetHeaderH3("2: " + std::to_string(probability::modesValue[1][0]) + "-" + std::to_string(probability::modesValue[1][1]) +
                    "-" + std::to_string(probability::modesValue[1][2]) + "-" + std::to_string(probability::modesValue[1][3]) + "\n");
    out.SetHeaderH3("3: " + std::to_string(probability::modesValue[2][0]) + "-" + std::to_string(probability::modesValue[2][1]) +
                    "-" + std::to_string(probability::modesValue[2][2]) + "-" + std::to_string(probability::modesValue[2][3]) + "\n");
    out.SetMatrix(modes, maxCarsInFlows[0] + 1);
    out.SetLine();
}

std::pair<double, double> computeProbabilitys(const double alpha, const double c) {
    double expectedValue = alpha * c * 1.0 + alpha * (1.0 - c) * 2.0;
    double resAlpha = -1.0, resC = -1.0;
    if (c >= 0.5) {
        double newAlpha = alpha - 0.001;
        int iterCount = 1;
        for (; newAlpha > 0.0; newAlpha -= 0.001) {
            double newC = c - 0.001;
            for (; newC > 0.0; newC -= 0.001) {
                double newExpectedValue = newAlpha * newC * 1.0 + newAlpha * (1.0 - newC) * 2.0;
                if (std::abs(newExpectedValue - expectedValue) <= 0.000001) {
                    resAlpha = newAlpha;
                    resC = newC;
                    //break;
                }
            }
            if (std::abs(resAlpha - resC) > 0.0000001)
                break;
            ++iterCount;
        }
    } else {
        double newAlpha = alpha + 0.001;
        int iterCount = 1;
        for (; newAlpha < 1.0; newAlpha += 0.001) {
            double newC = c + 0.001;
            for (; newC < 1.0; newC += 0.001) {
                double newExpectedValue = newAlpha * newC * 1.0 + newAlpha * (1.0 - newC) * 2.0;
                if (std::abs(newExpectedValue - expectedValue) <= 0.0001) {
                    resAlpha = newAlpha;
                    resC = newC;
                    //break;
                }
            }
            if (std::abs(resAlpha - resC) > 0.0000001)
                break;
            ++iterCount;
        }
    }
    return std::pair<double, double>(resAlpha, resC);
}

double computeArrivalProbability(const double alpha, const double c) {
    return alpha * c * 1.0 + alpha * (1.0 - c) * 2.0;
}

int main(int argc, char** argv) {
    /*std::cout << "----------------First mode--------------" << std::endl;
    // первый режим обслуживания
    probability::Flow f11(0);
    probability::Flow f12(1);
    probability::Flow::changeState(0);

    f11.PFormation();
    f11.QFormation();
    f11.GFormation();
    f11.HFormation();
    f12.PFormation();
    f12.QFormation();
    f12.GFormation();
    f12.HFormation();

    std::cout << "flow 1, matrix P:" << std::endl;
    f11.printP();
    std::cout << "flow 1, matrix Q:" << std::endl;
    f11.printQ();
    std::cout << "flow 2, matrix P:" << std::endl;
    f12.printP();
    std::cout << "flow 2, matrix Q:" << std::endl;
    f12.printQ();

    std::cout << "flow 1, matrix G:" << std::endl;
    f11.printG();
    std::cout << "flow 1, matrix H:" << std::endl;
    f11.printH();
    std::cout << "flow 2, matrix G:" << std::endl;
    f12.printG();
    std::cout << "flow 2, matrix H:" << std::endl;
    f12.printH();

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
    probability::Flow::changeState(1);

    f21.PFormation();
    f21.QFormation();
    f21.GFormation();
    f21.HFormation();
    f22.PFormation();
    f22.QFormation();
    f22.GFormation();
    f22.HFormation();

    std::cout << "flow 1, matrix P:" << std::endl;
    f21.printP();
    std::cout << "flow 1, matrix Q:" << std::endl;
    f21.printQ();
    std::cout << "flow 2, matrix P:" << std::endl;
    f22.printP();
    std::cout << "flow 2, matrix Q:" << std::endl;
    f22.printQ();

    std::cout << "flow 1, matrix G:" << std::endl;
    f21.printG();
    std::cout << "flow 1, matrix H:" << std::endl;
    f21.printH();
    std::cout << "flow 2, matrix G:" << std::endl;
    f22.printG();
    std::cout << "flow 2, matrix H:" << std::endl;
    f22.printH();

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

    std::cout << std::endl << std::endl << "----------------Third mode--------------" << std::endl;
    // третий режим обслуживания
    probability::Flow f31(0);
    probability::Flow f32(1);
    probability::Flow::changeState(2);

    f31.PFormation();
    f31.QFormation();
    f31.GFormation();
    f31.HFormation();
    f32.PFormation();
    f32.QFormation();
    f32.GFormation();
    f32.HFormation();

    std::cout << "flow 1, matrix P:" << std::endl;
    f31.printP();
    std::cout << "flow 1, matrix Q:" << std::endl;
    f31.printQ();
    std::cout << "flow 2, matrix P:" << std::endl;
    f32.printP();
    std::cout << "flow 2, matrix Q:" << std::endl;
    f32.printQ();

    std::cout << "flow 1, matrix G:" << std::endl;
    f31.printG();
    std::cout << "flow 1, matrix H:" << std::endl;
    f31.printH();
    std::cout << "flow 2, matrix G:" << std::endl;
    f32.printG();
    std::cout << "flow 2, matrix H:" << std::endl;
    f32.printH();

    trm = new double[statesCount * statesCount]{ 0.0 };
    probability::transitionMatrix(f31, f32, trm);
    printMatrix(trm, statesCount);
    mv = new double[statesCount] { 0.0 };
    probability::VectorOfMarginalProbability(trm, mv);;
    mv1 = new double[maxCarsInFlows[0] + 1];
    mv2 = new double[maxCarsInFlows[1] + 1];
    f31.getMarginalProbability(mv, mv1);
    f32.getMarginalProbability(mv, mv2);
    double* vectorH31 = new double[maxCarsInFlows[0] + 1]{};
    double* vectorH32 = new double[maxCarsInFlows[1] + 1]{};
    double* Z3 = new double[statesCount];
    f31.VectorH(vectorH31);
    printVector(vectorH31, maxCarsInFlows[0] + 1);
    f32.VectorH(vectorH32);
    printVector(vectorH32, maxCarsInFlows[1] + 1);
    probability::ExpectedValueOfRequestTime(vectorH31, vectorH32, Z3);
    printMatrix(Z3, maxCarsInFlows[0] + 1);
    std::cout << "Expected Value of the number of customers in 1-st flow " << f31.expectedValue(mv1) << std::endl;
    std::cout << "Expected Value of the number of customers in 2-nd flow " << f32.expectedValue(mv2) << std::endl;
    std::cout << "Variance of the number of customers in 1-st flow " << f31.variance(mv1) << std::endl;
    std::cout << "Variance of the number of customers in 2-nd flow " << f32.variance(mv2) << std::endl;
    std::cout << "Covariance of the number of customers " << probability::covariance(f31, f32, mv) << std::endl;
    delete[] mv2;
    delete[] mv1;
    delete[] mv;
    delete[] trm;
    delete[] vectorH31;
    delete[] vectorH32;
    delete[] Z3;

    return 0;*/

    int flowNumber;
    size_t maxCars[2]{};
    std::string experimentType;
    std::cout << "Write system parameters:" << std::endl
              << "Different sizes of the input bundle [y/n/custom]: "; std::cin >> experimentType;
    std::cout << "Number of flow: "; std::cin >> flowNumber;
    std::cout << "Max car in flow 1: "; std::cin >> maxCars[0];
    std::cout << "Max car in flow 2: "; std::cin >> maxCars[1];
    if (experimentType == "n") {
        std::string leadingParam;
        double start_alpha;
        double step_alpha;
        std::string is_new_politics;
        std::cout << "Name of the parameter: "; std::cin >> leadingParam;
        std::cout << "Write start param probability: "; std::cin >> start_alpha;
        std::cout << "Write param probability step: "; std::cin >> step_alpha;
        std::cout << "Do you want to change politics [y/n]: "; std::cin >> is_new_politics;
        if (is_new_politics == "y") {
            std::cout << "First policy: ";  std::cin >> probability::modesValue[0][0],
                std::cin >> probability::modesValue[0][1],
                std::cin >> probability::modesValue[0][2],
                std::cin >> probability::modesValue[0][3];
            std::cout << "Second policy: "; std::cin >> probability::modesValue[1][0],
                std::cin >> probability::modesValue[1][1],
                std::cin >> probability::modesValue[1][2],
                std::cin >> probability::modesValue[1][3];
            std::cout << "Third policy: ";  std::cin >> probability::modesValue[2][0],
                std::cin >> probability::modesValue[2][1],
                std::cin >> probability::modesValue[2][2],
                std::cin >> probability::modesValue[2][3];
        }
        if (start_alpha > 1.0) {
            std::cerr << "alpha > 1" << std::endl;
            return -1;
        }
        size_t iter = 0;
        size_t max_iter_count = (1.0 - start_alpha) / step_alpha;
        uint8* previousModes = new uint8[statesCount];
        double a[2]{};
        if (leadingParam == "alpha")
            a[1 - flowNumber] = probability::alpha[1 - flowNumber];
        else if (leadingParam == "beta")
            a[1 - flowNumber] = probability::beta[1 - flowNumber];
        else if (leadingParam == "c")
            a[1 - flowNumber] = probability::c[1 - flowNumber];
        for (; start_alpha < 1.0; start_alpha += step_alpha) {
            a[flowNumber] = start_alpha;
            if (leadingParam == "alpha")
                probability::setSystemConstants(2, maxCars, a, probability::beta, probability::c);
            else if (leadingParam == "beta")
                probability::setSystemConstants(2, maxCars, probability::alpha, a, probability::c);
            else if (leadingParam == "c")
                probability::setSystemConstants(2, maxCars, probability::alpha, probability::beta, a);

            uint8* modes = reinterpret_cast<uint8*>(operator new(sizeof(uint8) * statesCount));
            size_t iter_count = probability::HowardAlgorithm(modes, 0);
            printMatrix(modes, maxCarsInFlows[0] + 1);
            for (size_t i = 0; i < statesCount; ++i) {
                if (previousModes[i] != modes[i])
                    std::cout << "i - " << i << ", prev - "
                    << static_cast<int>(previousModes[i]) << ", curret - "
                    << static_cast<int>(modes[i]) << std::endl;
            }
            if (iter < max_iter_count - 1)
                std::memcpy(previousModes, modes, statesCount);
            if (leadingParam == "alpha")
                outputInFile(modes, iter_count, iter, 0, start_alpha);
            else if (leadingParam == "beta")
                outputInFile(modes, iter_count, iter, 1, start_alpha);
            else if (leadingParam == "c")
                outputInFile(modes, iter_count, iter, 2, start_alpha);

            for (size_t i = 0; i < statesCount; ++i)
                (modes + i)->~uint8();
            operator delete(modes);
            ++iter;
        }
        delete[] previousModes;
    } else if (experimentType == "y") {
        double alpha[4];
        double beta[2];
        double c[4];
        std::cout << "Write probability of client arrival: "; std::cin >> alpha[0], std::cin >> alpha[1];
        std::cout << "Write probability that client will be served: "; std::cin >> beta[0], std::cin >> beta[1];
        std::cout << "Write probability that client comes alone: "; std::cin >> c[0], std::cin >> c[1];
        //auto params1 = computeProbabilitys(alpha[0], c[0]);
        //auto params2 = computeProbabilitys(alpha[1], c[1]);
        double param1 = computeArrivalProbability(alpha[0], c[0]);
        double param2 = computeArrivalProbability(alpha[1], c[1]);
        alpha[2] = param1;// params1.first;
        alpha[3] = param2;// params2.first;
        c[2] = 1.0;//params1.second;
        c[3] = 1.0;//params2.second;
        std::cout << "Colculated new probability of client arrival: " << alpha[2] << "  " << alpha[3] << std::endl;
        std::cout << "Colculated new probability that client comes alone: " << c[2] << "  " << c[3] << std::endl;
        probability::setSystemConstants(flowNumber, maxCars, alpha, beta, c);

        uint8* modes = reinterpret_cast<uint8*>(operator new(sizeof(uint8) * statesCount));
        size_t iter_count = probability::HowardAlgorithm(modes, 0);
        printMatrix(modes, maxCarsInFlows[0] + 1);
        for (size_t i = 0; i < statesCount; ++i)
            (modes + i)->~uint8();
        OutputInHTML("DifferentTypes", modes);
        operator delete(modes);

        probability::setSystemConstants(flowNumber, maxCars, alpha + 2, beta, c + 2);
        modes = reinterpret_cast<uint8*>(operator new(sizeof(uint8) * statesCount));
        iter_count = probability::HowardAlgorithm(modes, 1);
        printMatrix(modes, maxCarsInFlows[0] + 1);
        for (size_t i = 0; i < statesCount; ++i)
            (modes + i)->~uint8();
        OutputInHTML("DifferentTypes", modes);
        operator delete(modes);
    } else if (experimentType == "custom") {
        double a[2]{}, b[2]{}, cc[2]{};
        std::string is_new_politics;
        std::cout << "alpha in flow 1: "; std::cin >> a[0];
        std::cout << "alpha in flow 2: "; std::cin >> a[1];
        std::cout << "beta in flow 1: "; std::cin >> b[0];
        std::cout << "beta in flow 2: "; std::cin >> b[1];
        std::cout << "c in flow 1: "; std::cin >> cc[0];
        std::cout << "c in flow 2: "; std::cin >> cc[1];
        std::cout << "Do you want to change politics [y/n]: "; std::cin >> is_new_politics;
        if (is_new_politics == "y") {
            std::cout << "First policy: ";  std::cin >> probability::modesValue[0][0],
                                            std::cin >> probability::modesValue[0][1],
                                            std::cin >> probability::modesValue[0][2],
                                            std::cin >> probability::modesValue[0][3];
            std::cout << "Second policy: "; std::cin >> probability::modesValue[1][0],
                                            std::cin >> probability::modesValue[1][1],
                                            std::cin >> probability::modesValue[1][2],
                                            std::cin >> probability::modesValue[1][3];
            std::cout << "Third policy: ";  std::cin >> probability::modesValue[2][0],
                                            std::cin >> probability::modesValue[2][1],
                                            std::cin >> probability::modesValue[2][2],
                                            std::cin >> probability::modesValue[2][3];
        }
        probability::setSystemConstants(2, maxCars, a, b, cc);

        uint8* modes = reinterpret_cast<uint8*>(operator new(sizeof(uint8) * statesCount));
        probability::HowardAlgorithm(modes, 0);
        printMatrix(modes, maxCarsInFlows[0] + 1);
        for (size_t i = 0; i < statesCount; ++i)
            (modes + i)->~uint8();
        OutputInHTML("Output", modes);
        operator delete(modes);
    }
    
    return 0;
}