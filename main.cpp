#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <atomic>
#include <numeric>
#include <cmath>
#include <algorithm>
#include <map>

std::mutex mtx; // Мьютекс для синхронизации вывода

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
        counts[static_cast<int>(num * 12)]++;
    }
    return std::all_of(counts.begin(), counts.end(), [&](int count) { return count == counts[0]; });
}

bool is_random(const std::vector<double>& normalized_sequence) {
    int increasing = 0, decreasing = 0;
    for (size_t i = 0; i < normalized_sequence.size() - 1; ++i) {
        if (normalized_sequence[i] < normalized_sequence[i + 1]) increasing++;
        if (normalized_sequence[i] > normalized_sequence[i + 1]) decreasing++;
    }
    return std::abs(increasing - decreasing) <= 5; // 5 - это примерное количество переходов, чтобы увеличить точность, нужно увеличить количество элементов в последовательности
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

double calculate_chi_square_for_randomness(const std::vector<int>& sequence) {
    // Подсчет серий разной длины
    std::map<int, int> series_count; // Словарь для хранения количества серий каждой длины
    int current_length = 1;
    for (size_t i = 1; i < sequence.size(); ++i) {
        if (sequence[i] == sequence[i - 1]) {
            current_length++;
        } else {
            series_count[current_length]++;
            current_length = 1;
        }
    }
    series_count[current_length]++; // Добавить последнюю серию

    // Вычисление хи-квадрат статистики
    double chi_square = 0.0;
    int total_series = sequence.size() - 1; // Общее количество серий
    for (const auto& [length, count] : series_count) {
        double expected_count = 2.0 / (length * (length + 1)) * total_series; // Ожидаемое количество серий данной длины
        chi_square += std::pow(count - expected_count, 2) / expected_count;
    }

    return chi_square;
}

// Функция для проверки равномерности - вычисление хи-квадрат статистики
double calculate_chi_square_for_uniformity(const std::vector<double>& normalized_sequence, int num_intervals) {
    std::vector<int> interval_counts(num_intervals, 0); // Счетчики для каждого интервала
    for (double value : normalized_sequence) {
        int index = static_cast<int>(value * num_intervals);
        if (index == num_intervals) index--; // Обработка краевого случая, когда значение равно 1
        interval_counts[index]++;
    }

    double expected_count = static_cast<double>(normalized_sequence.size()) / num_intervals; // Ожидаемое количество в каждом интервале
    double chi_square = 0.0;
    for (int count : interval_counts) {
        chi_square += std::pow(count - expected_count, 2) / expected_count;
    }

    return chi_square;
}

void find_parameters(int lambda_start, int lambda_end, int mu_start, int mu_end, int m_start, int m_end, int x0_start, int x0_end, std::atomic<bool>& found) {
    for (int x0 = x0_start; x0 <= x0_end && !found.load(); ++x0) {
        for (int lambda = lambda_start; lambda <= lambda_end && !found.load(); ++lambda) {
            for (int mu = mu_start; mu <= mu_end && !found.load(); ++mu) {
                for (int m = m_start; m <= m_end && !found.load(); ++m) {
                    auto sequence = generate_sequence(lambda, mu, m, x0, 900);
                    auto normalized_sequence = normalize_sequence(sequence, m);
                    if (is_uniform(normalized_sequence) && is_random(normalized_sequence) && is_independent(normalized_sequence)) {
                        mtx.lock();
                        std::cout << "Found valid parameters: Lambda: " << lambda << ", Mu: " << mu << ", M: " << m << ", X0: " << x0 << std::endl;
                        mtx.unlock();
                        found.store(true);
                        return;
                    }
                }
            }
        }
    }
}

int main() {
    const int num_threads = 5; // Количество потоков
    std::vector<std::thread> threads;
    std::atomic<bool> found(false);

    // Разделение диапазона параметров между потоками
    int lambda_range = 5; // Примерный диапазон для lambda
    int lambda_per_thread = lambda_range / num_threads;

    for (int i = 0; i < num_threads; ++i) {
        int lambda_start = 1 + i * lambda_per_thread;
        int lambda_end = (i + 1) * lambda_per_thread;
        threads.emplace_back(find_parameters, lambda_start, lambda_end, 5, 50, 5, 50, 1, 5, std::ref(found));
    }

    for (auto& t : threads) {
        t.join();
    }

    if (!found.load()) {
        std::cout << "No valid parameters found." << std::endl;
    }
    
    //------------------------------------------------------------------------------------------------------------------
    // Проверка тестов
    int lambda = 5, mu = 29, m = 37, x0 = 1;
    auto sequence = generate_sequence(lambda, mu, m, x0, 900);
    auto normalized_sequence = normalize_sequence(sequence, m);

    double chi_square_randomness = calculate_chi_square_for_randomness(sequence);
    double chi_square_uniformity = calculate_chi_square_for_uniformity(normalized_sequence, 12);

    std::cout << "Chi-Square for Randomness: " << chi_square_randomness << std::endl;
    std::cout << "Chi-Square for Uniformity: " << chi_square_uniformity << std::endl;


    return 0;
}
