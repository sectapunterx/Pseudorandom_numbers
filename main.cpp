#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <algorithm>

// Функция для генерации последовательности
std::vector<int> generate_sequence(int lambda, int mu, int m, int x0, int length = 600) {
    std::vector<int> sequence(length);
    sequence[0] = x0;
    for (int i = 1; i < length; ++i) {
        sequence[i] = (lambda * sequence[i - 1] + mu) % m;
    }
    return sequence;
}

// Функция для нормализации последовательности к интервалу (0, 1)
std::vector<double> normalize_sequence(const std::vector<int>& sequence, int m) {
    std::vector<double> normalized(sequence.size());
    for (size_t i = 0; i < sequence.size(); ++i) {
        normalized[i] = static_cast<double>(sequence[i]) / m;
    }
    return normalized;
}

// Функции для проверки тестов
bool is_uniform(const std::vector<double>& normalized_sequence) {
    std::vector<int> counts(12, 0); // 12 интервалов
    for (double num : normalized_sequence) {
        counts[static_cast<int>(num * 10)]++;
    }
    return std::all_of(counts.begin(), counts.end(), [&](int count) { return count == counts[0]; });
}

bool is_random(const std::vector<double>& normalized_sequence) {
    int increasing = 0, decreasing = 0;
    for (size_t i = 0; i < normalized_sequence.size() - 1; ++i) {
        if (normalized_sequence[i] < normalized_sequence[i + 1]) increasing++;
        if (normalized_sequence[i] > normalized_sequence[i + 1]) decreasing++;
    }
    return std::abs(increasing - decreasing) <= 5;
}

double calculate_correlation(const std::vector<double>& normalized_sequence) {
    double sum_xi_xi1 = 0.0;
    double sum_xi = std::accumulate(normalized_sequence.begin(), normalized_sequence.end(), 0.0);
    double sum_xi_squared = 0.0;
    for (size_t i = 0; i < normalized_sequence.size() - 1; ++i) {
        sum_xi_xi1 += normalized_sequence[i] * normalized_sequence[i + 1];
        sum_xi_squared += normalized_sequence[i] * normalized_sequence[i];
    }
    sum_xi_squared += normalized_sequence.back() * normalized_sequence.back(); // Add the last term

    double numerator = sum_xi_xi1 - sum_xi * sum_xi / normalized_sequence.size();
    double denominator = sum_xi_squared - sum_xi * sum_xi / normalized_sequence.size();

    return numerator / denominator;
}

bool is_independent(const std::vector<double>& normalized_sequence) {
    double correlation = calculate_correlation(normalized_sequence);
    return std::abs(correlation) < 0.05;
}

int main() {
    //int x0 = 1;

    for(int x0 = 2; x0 <= 100; ++x0) {
        for (int lambda = 1; lambda <= 5; ++lambda) { // Примерные диапазоны
            for (int mu = 100; mu <= 200; ++mu) {
                for (int m = 100; m <= 300; ++m) {
                    auto sequence = generate_sequence(lambda, mu, m, x0);
                    auto normalized_sequence = normalize_sequence(sequence, m);

                    if (is_uniform(normalized_sequence) && is_random(normalized_sequence) &&
                        is_independent(normalized_sequence)) {
                        std::cout << "Found valid parameters: " << std::endl;
                        std::cout << "X0: " << x0 << ", Lambda: " << lambda << ", Mu: " << mu << ", M: " << m << std::endl;
                        std::cout << "Uniform: " << is_uniform(normalized_sequence)
                                  << std::endl; // Uniform - равномерный
                        std::cout << "Random: " << is_random(normalized_sequence) << std::endl; // Random - случайный
                        std::cout << "Independent: " << is_independent(normalized_sequence)
                                  << std::endl; // Independent - независимый
                        return 0;
                    }
                }
            }
        }
    }

    std::cout << "No valid parameters found." << std::endl;
    return 0;
}