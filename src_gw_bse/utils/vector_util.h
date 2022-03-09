//
// Copyright (c) 2020 OpenAtom developers.
//
// File developed by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//
//
// File created by: Kayahan Saritas, saritaskayahan@gmail.com, Yale
//////////////////////////////////////////////////////////////////////////////////////

#ifndef VECTOR_UTIL_H_
#define VECTOR_UTIL_H_
#include <vector>

std::vector<double> linspace(const double, const double, const int);
std::vector<int> range(const int, const int, const int);
std::vector<int> range(const int, const int);
std::vector<int> range(const int);

// Elementwise operations
std::vector<double> operator+(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator*(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator*(const std::vector<std::vector<double>>&, const std::vector<double>&);
std::vector<double> operator/(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator-(const std::vector<double>&, const std::vector<double>&);
std::vector<double> exp(const std::vector<double>&);
std::vector<double> pow(const std::vector<double>&, const double&);
std::vector<double> pow(const double&, const std::vector<double>&);
// std::vector<double> log10(const std::vector<double>&);
// std::vector<double> log(const std::vector<double>&);
std::vector<double> sqrt(const std::vector<double>&);
std::vector<double> zeros(const int);

// Operations with a constant
std::vector<double> operator+(const std::vector<double>&, const double&);
std::vector<double> operator+(const double&, const std::vector<double>&);
std::vector<double> operator*(const std::vector<double>&, const double&);
std::vector<double> operator*(const double&, const std::vector<double>&);
std::vector<double> operator/(const std::vector<double>&, const double&);
std::vector<double> operator/(const double&, const std::vector<double>&);
std::vector<double> operator-(const std::vector<double>&, const double&);
std::vector<double> operator-(const double&, const std::vector<double>&);

// Basic vector math
double dot(const std::vector<double> &a, const std::vector<double> &b);
std::vector<double> cross(const std::vector<double> &a, const std::vector<double> &b);
double norm(const std::vector<double> &a);
double sum(const std::vector<double> &a);
std::vector<double> abs(const std::vector<double>&);
std::vector<double> sum(const std::vector<std::vector<double> > &a, int axis);

// Utilities
void print_vector(std::vector<double>);
void print_vector(std::vector<int>);
std::vector<double> da2v(double* x, int size);  // double array to vector

// std::vector<double> mult(double, const std::vector<double>&);
// std::vector<double> mult(const std::vector<double>&, const std::vector<double>&);

#endif  // VECTOR_UTIL_H_

