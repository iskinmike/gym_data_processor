#ifndef TRACK_HANDLER_HPP
#define TRACK_HANDLER_HPP


#include <vector>
#include <string>
#include <map>
#include <math.h>
#include <algorithm>
#include <set>

#include <sstream>
#include <iostream>
#include <fstream>

//#include "marker_tracker.hpp"
//#include "geometry/geometry.h"


class point
{
public:
    double x,y;
    point(): x(0), y(0) {}
    point(double x, double y): x(x), y(y) {}
    double dist(point other);
    std::string print();

    point operator +(const point &other) const;
    point operator -(const point &other) const;
    point operator *(double d)       const;
    point& operator =(point&& other);
    const point& operator =(const point& other);
    point(const point &other);
    point(point&& other);
    ~point() {}
};

struct point_data{
    point pos;
    int64_t step;
    int64_t time;
    std::vector<float> matrix;
    point_data() : pos(0,0){}
    point_data(const point_data& other) {
        pos = other.pos;
        step = other.step;
        time = other.time;
        matrix = other.matrix;
    }
};

struct delta_data{
    double delta;
    int line;
    int point;
};


class vector_2f
{
    point begin;
    point end;
    point vec;
    double len;
    double x_len;
    double y_len;
public:
    vector_2f():begin(0,0), end(0,0), vec(0,0){}
    explicit vector_2f(point begin, point end) : begin(begin), end(end) {
        len = begin.dist(end);
        vec = end - begin;
        x_len = vec.x;
        y_len = vec.y;
    }
//    vector_2f(){

//    }

//    point get_begin() const;
//    point get_end() const;
//    double get_len() const;
//    double get_x_len() const;
//    double get_y_len() const;
};


class oth {
public:
    double x,y;
    oth(): x(0), y(0){}
    oth(double x, double y): x(x), y(y){}
    ~oth(){}
};

class testc{
//    point begin;
//    point end;
//    double len;
//    double x_len;
//    double y_len;
    oth begin;
    oth end;
    oth vec;
    double tt;
public:
//    testc(): begin(0,0), end(0,0), vec(0,0), len(0) {}
    testc(): begin(0,0), end(0,0), vec(0,0), tt(0.0f) {}
//    testc(): begin(0,0), end(0,0), tt(0.0f) {}
//    testc(): begin(0,0){}
};


struct test{
    int i;
};

struct shift_data{
//    vector_2f shift_vec;
    testc tt;
    double shift;
    int start_pos;
    int end_pos;
    double cs_pos_x_1;
    double cs_pos_x_2;
    double cs_pos_y_1;
    double cs_pos_y_2;
//    shift_data() : shift_vec() {}
};

struct derivative_2d_data{
//    point pos;
    double x_coord_value;
    double y_coord_value;
    double step;
//    explicit derivative_2d_data(): pos(0,0){
//        x_coord_value = 0;
//        y_coord_value = 0;
//        step = 0;
//    }
//    derivative_2d_data& operator =(const derivative_2d_data& other) {
//        pos = other.pos;
//        step = other.step;
//        x_coord_value = other.x_coord_value;
//        y_coord_value = other.y_coord_value;
//        return *this;
//    }
//    derivative_2d_data& operator =(derivative_2d_data&& other) {
//        pos = other.pos;
//        step = other.step;
//        x_coord_value = other.x_coord_value;
//        y_coord_value = other.y_coord_value;
//        return *this;
//    }
//    derivative_2d_data(const derivative_2d_data& other) {
//        pos = other.pos;
//        step = other.step;
//        x_coord_value = other.x_coord_value;
//        y_coord_value = other.y_coord_value;
//    }
//    derivative_2d_data(derivative_2d_data&& other) {
//        pos = other.pos;
//        step = other.step;
//        x_coord_value = other.x_coord_value;
//        y_coord_value = other.y_coord_value;
//    }
};




struct marker_description {
    int id;
    point pos;
//    double pos_x;
//    double pos_y;
    int64_t frame_time;
    int count;
    bool empty;
    std::vector<float> matrix;
};

struct coord_system_origin{
    marker_description first_marker;
    marker_description second_marker;
};

struct time_slice {
//    double time;
    int64_t frame_time;
    // может лучше строка как идентификатор?
    std::map<int, marker_description> markers;
    coord_system_origin origin;
};

struct marker_amp_mm_data{
    int pos;
    marker_amp_mm_data(int pos) : pos(pos){}
};

enum class extremum_type {
    min, max, plateu, valley
};

enum class coord{
    x,y
};

struct extremum_point {
    extremum_type type;
    int pos;
    extremum_point(int pos, extremum_type type) : pos(pos), type(type) {}
};

struct marker_amplitude_data {
//    std::vector<marker_amp_mm_data> minimums;
//    std::vector<marker_amp_mm_data> maximums;
    coord coord_type;
    int minima_count;
    int maxima_count;
    std::vector<extremum_point> extremums;
    std::vector<shift_data> shifts;
    marker_amplitude_data(): minima_count(0), maxima_count(0){}
};

struct marker_amplitudes_xy{
    marker_amplitude_data amp_x;
    marker_amplitude_data amp_y;
};

class track_handler
{
    std::vector<std::vector<point_data>> plot;
    std::vector<std::vector<point_data>> timeline;

    std::vector<coord_system_origin> origin_line;

    double tracker_evaluation_error;
public:
    track_handler() {}

    // Добавить задание данных о системе отсчета. должно быть достаточно 2-х меток
    void setup_origing_track();

    void append_origin(std::vector<point_data> origin_points);


//    double calc_shift(point_data p1, point_data p2);

    void append(std::vector<point_data> points_on_step);

    size_t get_lines_count();

    std::vector<point_data> get_line(int n);

    void dump_all_to_file(const std::string& path);

    int get_right_pos(int i, const std::vector<point_data>& data);

    void calc_amp();

    marker_amplitude_data calc_amplitude_for_marker(const std::vector<point_data> &marker_track, coord coordinate);
    marker_amplitude_data calc_amplitude_for_marker_from_derivative(const std::vector<derivative_2d_data>& derivative_track, const std::vector<point_data> &marker_track, coord coordinate);
    marker_amplitudes_xy calc_marker_amplitudes_xy(const std::vector<point_data> &marker_track);

    void calc_derivatives();
};





#endif  /* TRACK_HANDLER_HPP */
