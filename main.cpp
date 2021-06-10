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
        double param1 = computeArrivalProbability(alpha[0], c[0]);
        double param2 = computeArrivalProbability(alpha[1], c[1]);
        alpha[2] = param1;
        alpha[3] = param2;
        c[2] = 1.0;
        c[3] = 1.0;
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