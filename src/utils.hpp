
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
};


#endif  /* UTILS_HPP */
