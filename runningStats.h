#pragma once
#include <vector>
#include <cmath>
#include <iostream>

void updateStatistics(const std::vector<double>& values, std::vector<double>& sum, std::vector<double>& sum_of_squares, size_t& count) {
    count++;
    for (size_t i = 0; i < values.size(); ++i) {
        sum[i] += values[i];
        sum_of_squares[i] += values[i] * values[i];
    }
}

std::vector<double> calculateMean(const std::vector<double>& sum, size_t count) {
    std::vector<double> mean_values(sum.size(), 0.0);
    for (size_t i = 0; i < sum.size(); ++i) {
        mean_values[i] = sum[i] / count;
    } 
    return mean_values;
}       

std::vector<double> calculateStandardError(const std::vector<double>& sum, const std::vector<double>& sum_of_squares, size_t count) {
    std::vector<double> std_err(sum.size(), 0.0);
    for (size_t i = 0; i < sum.size(); ++i) {
        double mean = sum[i] / count;
        double variance = (sum_of_squares[i] / count) - (mean * mean);
        std_err[i] = std::sqrt(variance / count);
    }
    return std_err;
}
