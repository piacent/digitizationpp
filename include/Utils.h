#pragma once

#include <vector>
#include <string>
#include <filesystem>
#include <algorithm>

namespace Utils {

std::string resolvePath(const std::string& relativePath);

std::vector<std::string> splitString(const std::string& input, char delimiter);

double angleBetween(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> rotateByAngleAndAxis(const std::vector<double>& vec, double angle, const std::vector<double>& axis);

double roundUpToEven(double value);
std::vector<double> arange(double start, double stop, double step);

}
