//////////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2020 OpenAtom developers.
//
// File developed by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//
//
// File created by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numeric>
#include <iostream>
#include <vector>
#include "./util.h"
#include "./vector_util.h"

std::vector<double> linspace(const double a, const double b, const int N) {
    assert(N>1); 
    double h = (b - a) / (1.0*(N-1));
    // std::vector<double> xs(N);
    // std::vector<double>::iterator x;
    // double val;
    // for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    //     *x = val;
    // return xs;
    std::vector<double> result;
    result.clear();
    for (int i = 0; i < N; i++) {
        result.push_back(a+(i*1.0)*h);
    }
    return result;
}

void print_vector(const std::vector<double> v) {
    for (int i = 0; i < v.size(); i++) {
        printf("%f ", v[i]);
    }
    printf("\n");
}

void print_vector(const std::vector<int> v) {
    for (int i = 0; i < v.size(); i++) {
        printf("%d ", v[i]);
    }
    printf("\n");
}

std::vector<int> range(const int start, const int stop, const int step) {
    if (step == 0) {
        printf("step for range must be non-zero");
        exit(1);
    }
    std::vector<int> result;
    int i = start;
    while ((step > 0) ? (i < stop) : (i > stop)) {
        result.push_back(i);
        i += step;
    }
    return result;
}

std::vector<int> range(const int start, const int stop) {
    return range(start, stop, 1);
}

std::vector<int> range(const int stop) {
    return range(0, stop, 1);
}

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    assert(a.size() == b.size());
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::plus<double>());
    return result;
}
std::vector<double> operator+(const std::vector<double>& a, const double & b) {
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), std::back_inserter(result), std::bind(std::plus<double>(), std::placeholders::_1, b));
    return result;
}

std::vector<double> operator+(const double & b, const std::vector<double>& a) {
    return a + b;
}

std::vector<double> operator* (const std::vector<std::vector<double>> &a, const std::vector<double> &b) {
    // Matrix-array multiplication, makes m by 1 array
    int m = a.size();
    int n = a[0].size();
    int k = b.size();
    assert(n == k);
    std::vector<double> result(m, 0.);
    for (int i: range(0, m)) {
        for (int j: range(0, n)) {
            result[i] += a[i][j]*b[j];
        }
    }
    return result;
}

std::vector<double> operator* (const std::vector<double> &a, const std::vector<double> &b) {
    assert(a.size() == b.size());
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(result), std::multiplies<double>());
    return result;
}
std::vector<double> operator* (const double &b, const std::vector<double> &a) {
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), std::back_inserter(result), std::bind(std::multiplies<double>(), std::placeholders::_1, b));
    return result;
}
std::vector<double> operator* (const std::vector<double> &a, const double &b) {
    std::vector<double> result;
    result = b * a;
    return result;
}

std::vector<double> operator/(const std::vector<double> &a, const std::vector<double> &b) {
    assert(a.size() == b.size());
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::divides<double>());
    return result;
}
std::vector<double> operator/(const std::vector<double> &a, const double &b) {
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), std::back_inserter(result), std::bind(std::divides<double>(), std::placeholders::_1, b));
    return result;
}
std::vector<double> operator/(const double a, const std::vector<double> &b) {
    std::vector<double> result;
    result.reserve(b.size());
    std::transform(b.begin(), b.end(), std::back_inserter(result), std::bind(std::divides<double>(), std::placeholders::_1, a));
    return result;
    //return (1./a) * b;
}


std::vector<double> operator-(const std::vector<double> &a, const std::vector<double>&b) {
    assert(a.size() == b.size());
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(), std::back_inserter(result), std::minus<double>());
    return result;
}
std::vector<double> operator-(const std::vector<double>& a, const double & b) {
    std::vector<double> result;
    result.reserve(a.size());
    std::transform(a.begin(), a.end(), std::back_inserter(result), std::bind(std::minus<double>(), std::placeholders::_1, b));
    return result;
}

std::vector<double> operator-(const double & b, const std::vector<double>& a) {
    return -1.0 * (a - b);
}

std::vector<double> zeros(const int len) {
    if (len <= 0) {
        printf("length of a vector can't be <=0");
        exit(1);
    }
    std::vector<double> result;
    result.reserve(len);
    return (0.0 * result);
}

double norm(const std::vector<double>& a) {
    double sq_sum = 0.0;
    int len = a.size();
    for (int i = 0; i < len; i++) {
        sq_sum += pow(a[i], 2);
    }
    double result = sqrt(sq_sum);
    return result;
}

double sum(const std::vector<double>& a) {
    return accumulate(a.begin(), a.end(), 0., std::plus<double>());
}

std::vector<double> abs(const std::vector<double> &a) {
    int len = a.size();
    std::vector<double> result;
    // result.reserve(len);
    for (int i = 0; i < len; i++) {
        result.push_back(abs(a[i]));
    }
    return result;
}


std::vector<double> sum(const std::vector<std::vector<double>>& a) {
    int size = a.size();
    std::vector<double> result;
    for (std::vector<double> row:a) {
        result = result + row/size;
    }
    return result;
}

std::vector<double> exp(const std::vector<double> &a) {
    int len = a.size();
    std::vector<double> result;
    result.reserve(len);
    for (int i = 0; i < len; i++) {
        result.push_back(exp(a[i]));
    }
    return result;
}

// std::vector<double> log10(const std::vector<double> &a) {
//     int len = a.size();
//     std::vector<double> result;
//     // result.reserve(len);
//     for (int i = 0; i < len; i++) {
//         result.push_back(log10(a[i]));
//     }
//     return result;
// }

// std::vector<double> log(const std::vector<double> &a) {
//     int len = a.size();
//     std::vector<double> result;
//     // result.reserve(len);
//     for (int i = 0; i < len; i++) {
//         result.push_back(log(a[i]));
//     }
//     return result;
// }

std::vector<double> sqrt(const std::vector<double> &a) {
    int len = a.size();
    std::vector<double> result;
    result.reserve(len);
    for (int i = 0; i < len; i++) {
        result.push_back((a[i]));
    }
    return result;
}

std::vector<double> pow(const std::vector<double> &a, const double &x) {
    int len = a.size();
    std::vector<double> result;
    result.reserve(len);
    for (int i = 0; i < len; i++) {
        result.push_back(pow(a[i], x));
    }
    return result;
}

std::vector<double> pow(const double &x, const std::vector<double> &a) {
    int len = a.size();
    std::vector<double> result;
    result.reserve(len);
    for (int i = 0; i < len; i++) {
        result.push_back(pow(x, a[i]));
    }
    return result;
}
    
std::vector<double> da2v(double* x, int n) {
    // Double array to vector
    std::vector<double> result;
    for (int i = 0; i < n; i++) {
        result.push_back(x[i]);
    }
    return result;
}

std::vector<double> cross(const std::vector<double> &a, const std::vector<double> &b) {
    assert(a.size() == 3);
    assert(b.size() == 3);
    std::vector<double> result;
    result.push_back(a[1] * b[2] - a[2] * b[1]);
    result.push_back(-(a[0] * b[2] - a[2] * b[0]));
    result.push_back(a[0] * b[1] - a[1] * b[0]);
    return result;
}

double dot(const std::vector<double> &a, const std::vector<double> &b) {
    double result;
    result = sum(a * b);
    return result;
}

// std::vector<double> mult(double b, const std::vector<double> &a) {
//     std::vector<double> result;
//     result.reserve(a.size());
//     // std::transform(a.begin(), a.end(), a.begin(),
//     //        std::back_inserter(result), std::bind(std::multiplies<double>(), std::placeholders::_1, b));
//     for (int i = 0; i < a.size(); i++) {
//         result.push_back(b*a[i]);
//     }
//     return result;
// }

// std::vector<double> mult(const std::vector<double> &a, const std::vector<double> &b) {
//     assert(a.size() == b.size());
//     std::vector<double> result;
//     result.reserve(a.size());
//     // std::transform(a.begin(), a.end(), a.begin(),
//     //        std::back_inserter(result), std::bind(std::multiplies<double>(), std::placeholders::_1, b));
//     for (int i =0; i < a.size(); i++) {
//         result.push_back(b[i]*a[i]);
//     }
//     return result;
// }
