#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <chrono>
#include <ratio>
#include <string>
#include <vector>

using namespace std;
using vec2d = std::vector<std::vector<double>>;

std::vector<double> gen_vec(int size) {
  std::vector<double> mx;
  mx.reserve(size * size);

  for (int i = 0; i < size * size; ++i)
  {
    mx.push_back((double)rand() / (double)RAND_MAX);
  }

  return mx;
}

vec2d gen_vec2d(int size) {
  vec2d mx;
  mx.reserve(size);

  for (int i = 0; i < size; ++i){
    std::vector<double> v;
    v.reserve(size);
    for (int j = 0; j < size; ++j)
    {
        v.push_back((double)rand() / (double)RAND_MAX);
    }

    // push back above one-dimensional vector
    mx.push_back(v);
  }

  return mx;
}

double * gen_vec_double(int size, string mem_type) {
  double *mx;
  if (mem_type.compare("malloc") == 0) {
    mx = (double *)malloc((size * size) * sizeof(double));
  } else if (mem_type.compare("new") == 0) {
    mx = new double[size * size];
  } else {
      throw std::invalid_argument("Invalid value for 'implement'");
  }

  for (int i = 0; i < size * size; ++i)
  {
    mx[i] = (double)rand() / (double)RAND_MAX;
  }

  return mx;
}

double * mult_matrix_double(double *m1, double *m2, int size, string mem_type) {
  double *mult;

  if (mem_type.compare("malloc") == 0) {
    mult = (double *)malloc((size * size) * sizeof(double));
  } else if (mem_type.compare("new") == 0) {
    mult = new double[size * size];
  } else {
      throw std::invalid_argument("Invalid value for 'implement'");
  }

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
      mult[i * size + j] = 0;
      for (int k = 0; k < size; ++k)
      {
        mult[i * size + j] += m1[i * size + k] * m2[k * size + j];
      }
    }
  }

  return mult;
}

std::vector<double> mul_vec_matrix(std::vector<double> m1, std::vector<double> m2, int size) {
  std::vector<double> mult;
  mult.reserve(size * size);

  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j) {
    //   mult[i * size + j] = 0;
      mult.push_back(0);
      for (int k = 0; k < size; ++k)
      {
        mult.at(i * size + j) += m1.at(i * size + k) * m2.at(k * size + j);
      }
    }
  }

  return mult;
}

vec2d mul_vec2d_matrix(const vec2d& m1, const vec2d& m2) {
    int n = m1.size();
    vec2d mult(n, std::vector<double>(n));

    // allocate the result matrix
    mult.resize(n);
    for(vec2d::iterator it = mult.begin(); it != mult.end(); ++it){
        it->resize(m1[0].size());
    }

    for (int i = 0; i < n; i++){
        for (int j = 0 ; j < n ; j++){
            for (int k = 0; k < n; k++)
            {
                mult[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }

    return mult;
}

int evaluate_multiplication(int N, std::string mem_type) {
    if (mem_type.compare("malloc") == 0 || mem_type.compare("new") == 0) {
        double *mx1, *mx2, *mx_mult;
        /**
         * start the cronometer
         */
        auto t1 = std::chrono::high_resolution_clock::now();
        mx1 = gen_vec_double(N, mem_type);
        mx2 = gen_vec_double(N, mem_type);
        mx_mult = mult_matrix_double(mx1, mx2, N, mem_type);

        /**
         * end the cronometer
         */
        auto t2 = std::chrono::high_resolution_clock::now();
        /* Getting number of milliseconds as an integer. */
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        if (mem_type.compare("malloc") == 0) {
          cout << "[MALLOC] Runtime: " << duration.count() << " milliseconds"<< endl;
        } else {
          cout << "[NEW] Runtime: " << duration.count() << " milliseconds"<< endl;
        }
    }
    else if (mem_type.compare("vector") == 0)
    {
        std::vector<double> mx1, mx2, mx_mult;
        /**
         * start the cronometer
         */
        auto t1 = std::chrono::high_resolution_clock::now();
        mx1 = gen_vec(N);
        mx2 = gen_vec(N);
        mx_mult = mul_vec_matrix(mx1, mx2, N);
        /**
         * end the cronometer
         */
        auto t2 = std::chrono::high_resolution_clock::now();
        /* Getting number of milliseconds as an integer. */
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        cout << "[VECTOR] Runtime: " << duration.count() << " milliseconds"<< endl;
    }
    else if (mem_type.compare("vector2d") == 0)
    {
        /**
         * start the cronometer
         */
        auto t1 = std::chrono::high_resolution_clock::now();
        vec2d mx1, mx2, mx_mult;
        mx1 = gen_vec2d(N);
        mx2 = gen_vec2d(N);
        mx_mult = mul_vec2d_matrix(mx1, mx2);
        /**
         * end the cronometer
         */
        auto t2 = std::chrono::high_resolution_clock::now();
        /* Getting number of milliseconds as an integer. */
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

        cout << "[VECTOR2D] Runtime: " << duration.count() << " milliseconds"<< endl;
    }
    else
    {
        throw std::invalid_argument("Invalid value for 'mem_type'");
    }
    // mult_malloc(N);
    return 0;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
      cout << "Only N (matrix size) is expected.\n";
      return 0;
  }

  int N = atoi(argv[1]);

  if (N == 0) {
      cout << "N has invalid value.\n";
      return 0;
  }

  evaluate_multiplication(N, "malloc");
  evaluate_multiplication(N, "new");
  evaluate_multiplication(N, "vector");
  evaluate_multiplication(N, "vector2d");

  return 0;
}