#include "utils.hpp"
#include <cmath>
#include <sstream>

point::point(): x(0), y(0) {
}

point::point(double x, double y): x(x), y(y) {}

double point::dist(point other){
    double delta;
    delta = sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
    return fabs(delta);
}

point point::operator +(const point &other) const { return point(x+other.x, y+other.y); }

point point::operator -(const point &other) const { return point(x-other.x, y-other.y); }

point point::operator *(double d) const { return point(x*d, y*d); }

const point &point::operator =(const point &other) {
    x = other.x;
    y = other.y;
    return *this;
}

point &point::operator =(point &&other) {
    x = other.x;
    y = other.y;
    other.x = 0;
    other.y = 0;
    return *this;
}

point::point(point &&other) {
    x = other.x;
    y = other.y;
}

point::point(const point &other) {
    x = other.x;
    y = other.y;
}

std::string point::print(){
    std::stringstream str;
    str << "(" << x << ", " << y << ")";
    return str.str();
}

vector_2f::vector_2f():begin(0,0), end(0,0), vec(0,0), len(0), x_len(0), y_len(0){}

vector_2f::vector_2f(point begin, point end) : begin(begin), end(end) {
    len = begin.dist(end);
    vec = end - begin;
    x_len = vec.x;
    y_len = vec.y;
}

derivative_2d_data::derivative_2d_data():pos(0,0), x_coord_value(0), y_coord_value(0), step(0){
    step = 0;
}

derivative_2d_data::derivative_2d_data(const derivative_2d_data &other) {
    pos = other.pos;
    x_coord_value = other.x_coord_value;
    step = other.step;
    y_coord_value = other.y_coord_value;
}

derivative_2d_data::derivative_2d_data(derivative_2d_data &&other) {
    pos = other.pos;
    x_coord_value = other.x_coord_value;
    step = other.step;
    y_coord_value = other.y_coord_value;
}
