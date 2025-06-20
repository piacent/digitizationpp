#include "Utils.h"
#include <sstream>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <filesystem>

namespace Utils {

std::vector<std::string> splitString(const std::string& input, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream stream(input);
    std::string token;
    while (getline(stream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

double angleBetween(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != 3 || v2.size() != 3) throw std::invalid_argument("Both vectors must be 3D.");
    double dot = std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
    double len1 = std::sqrt(std::inner_product(v1.begin(), v1.end(), v1.begin(), 0.0));
    double len2 = std::sqrt(std::inner_product(v2.begin(), v2.end(), v2.begin(), 0.0));
    return std::acos(dot / (len1 * len2));
}

std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != 3 || b.size() != 3) throw std::invalid_argument("Both vectors must be 3D.");
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

std::vector<double> rotateByAngleAndAxis(const std::vector<double>& vec, double angle, const std::vector<double>& axis) {
    if (vec.size() != 3 || axis.size() != 3) throw std::invalid_argument("Both vectors must be 3D.");
    
    std::vector<double> cross = crossProduct(axis, vec);
    double dot = std::inner_product(axis.begin(), axis.end(), vec.begin(), 0.0);

    std::vector<double> result(3);
    for (int i = 0; i < 3; ++i) {
        result[i] = vec[i] * cos(angle) +
                    cross[i] * sin(angle) +
                    axis[i] * dot * (1.0 - cos(angle));
    }
    return result;
}

double roundUpToEven(double value) {
    int intVal = static_cast<int>(std::ceil(value));
    return (intVal % 2 == 0) ? intVal : intVal + 1;
}

std::string resolvePath(const std::string& relativePath) {
    return std::filesystem::absolute(std::filesystem::path(relativePath)).string();
}

}
