
#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>

class point
{
public:
    double x,y;
    point();
    point(double x, double y);
    double dist(const point &other);
    std::string print();
    point operator +(const point &other) const;
    point operator -(const point &other) const;
    point operator *(double d)       const;
//    point& operator =(point&& other);
//    const point& operator =(const point& other);
//    point(const point& other);
//    point(point&& other);
};

class vector_2f
{
    point begin;
    point end;
public:
    point vec;
    double len;
    double x_len;
    double y_len;
    vector_2f();
    vector_2f(point begin_inp, point end_inp);
};

struct derivative_2d_data{
    double x_coord_value;
    double y_coord_value;
    double time;
    point pos;
    int step;
//    derivative_2d_data();
//    derivative_2d_data(const derivative_2d_data& other);
//    derivative_2d_data(derivative_2d_data&& other);
};


#endif  /* UTILS_HPP */
