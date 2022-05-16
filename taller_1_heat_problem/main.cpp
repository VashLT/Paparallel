#include <iostream>
#include <cmath>

using namespace std;

const int MIN_TEMP = 0;
const int MAX_TEMP = 100;
int L, STEPS, TIME;

/**
 * @brief Prints bar temperature points
 * 
 * @param bar 
 */
void print_bar(double* bar) {
    for (int i = 0; i < L; i++) {
        printf("%f ", bar[i]);
    }
    printf("\n");
}

/**
 * @brief get value of flags
 * 
 * @param argc number of flags
 * @param argsv flags values
 * @param flag target flag
 * @return int 
 */
int parse_flag(int argc, char* argsv[], string flag) {
    bool arg = false;
    int n = -1;

    for (int i = 0; i < argc; i++) {
        if (flag.compare(argsv[i]) == 0 && !arg) {
            n = stoi(argsv[i + 1]);
            arg = true;
        } else {
            continue;
        }
    }

    if (n == -1) {
        throw invalid_argument("Missing flag (" + flag + ").");
    }

    return n;
}

/**
 * @brief generates initial points for bar
 * 
 * @param fn temperature function that defines initial conditions
 * @return double* 
 */
double* gen_bar(double (*fn)(int, double)) {
    double* bar = (double*)malloc(L * sizeof(double));

    for (int i = 0; i < L; i++) {
        bar[i] = fn(i, 0);
    }

    return bar;
}

/**
 * @brief Heat function
 * 
 * @param x 
 * @param t 
 * @return double 
 */
double T(int x, double _t) {
    if (x == 0) {
        return MAX_TEMP;
    } else if (x == L - 1) {
        return sqrt(MAX_TEMP);
    } else {
        return MIN_TEMP;
    }
}

/**
 * @brief apply heat formula to bar points
 * 
 * @param src_bar current bar state
 * @return double* 
 */
double* compute_heat(double* src_bar) {
    double* bar = (double*)malloc(L * sizeof(double));

    for (int i = 0; i < L; i++) {
        bar[i] = (src_bar[i - 1] + src_bar[i + 1]) / 2;
    }

    free(src_bar);

    return bar;
}

/**
 * @brief iteratively solve heat equation
 * 
 * @param src_bar bar with set initial conditions
 * @param fn heat funciton
 * @return double* 
 */
double* solve_heat_problem(double* src_bar, double(*fn)(int, double)) {
    double* bar = src_bar;

    const double dt = TIME / STEPS;

    for (double t = 0; t < TIME; t = t + dt) {
        bar = compute_heat(bar);

        for (int j = 0; j < L; j++) {
            double temperature = fn(j, t);
            if (temperature > 0) {
                bar[j] = temperature;
            }
        }
    }

    return bar;
}

int main(int argc, char* argsv[]) {
    const string size_flag = "-n";
    const string time_flag = "-t";
    const string steps_flag = "-s";

    L = parse_flag(argc, argsv, size_flag);
    TIME = parse_flag(argc, argsv, time_flag);
    STEPS = parse_flag(argc, argsv, steps_flag);

    double *bar = gen_bar(&T);

    printf("Initial bar's heat map:\n");
    print_bar(bar);

    printf("Final bar's heat map:\n");
    print_bar(solve_heat_problem(bar, &T));

    return 0;
}