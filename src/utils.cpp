#include "utils.hpp"
#include <cmath>
#include <sstream>

point::point(): x(0), y(0) {
}

point::point(double x, double y): x(x), y(y) {}

double point::dist(const point& other){
    double delta;
    delta = sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
    return fabs(delta);
}

point point::operator +(const point &other) const { return point(x+other.x, y+other.y); }

point point::operator -(const point &other) const { return point(x-other.x, y-other.y); }

point point::operator *(double d) const { return point(x*d, y*d); }

std::string point::print(){
    std::stringstream str;
    str << "(" << x << ", " << y << ")";
    return str.str();
}

vector_2f::vector_2f():begin(0,0), end(0,0), vec(0,0), len(0), x_len(0), y_len(0){}

vector_2f::vector_2f(point begin_inp, point end_inp) : begin(begin_inp), end(end_inp) {
    len = begin.dist(end);
    vec = end - begin;
    x_len = vec.x;
    y_len = vec.y;
}
