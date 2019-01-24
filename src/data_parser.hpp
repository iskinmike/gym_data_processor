#ifndef DATA_PARSER_HPP
#define DATA_PARSER_HPP


#include <string>
#include <vector>

#include "utils.hpp"

struct shift_data{
    vector_2f shift_vec;
    double shift;
    int start_pos;
    int end_pos;
    shift_data() :/*shift_vec(), */shift(0), start_pos(0), end_pos(0){}
};

struct markers_data {
    point pos;
    int64_t time;
    int64_t step;
    markers_data(): pos(0,0), time(0), step(0) {}
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
    extremum_point(int pos, extremum_type type) : type(type), pos(pos) {}
    extremum_point():type(extremum_type::valley), pos(0) {}
};

struct marker_amplitude_data {
    coord coord_type;
    int minima_count;
    int maxima_count;
    std::vector<extremum_point> extremums;
    std::vector<shift_data> shifts;
    marker_amplitude_data(): coord_type(coord::x), minima_count(0), maxima_count(0) , extremums{}, shifts{}{}
};

struct marker_amplitudes_xy{
    marker_amplitude_data amp_x;
    marker_amplitude_data amp_y;
    marker_amplitudes_xy(): amp_x{}, amp_y{}{}
};

struct marker_track_data {
    marker_amplitudes_xy amplitudes;
    std::vector<derivative_2d_data> first_derivative;
    std::vector<derivative_2d_data> second_derivative;
    std::vector<markers_data> track;
    marker_track_data(): amplitudes{}, first_derivative{},second_derivative{}, track{} {}
};




class data_parser
{
//    std::string filepath;
    std::vector<std::vector<markers_data>> marker_line;
    std::vector<markers_data> zero_marker_line;
    std::vector<markers_data> direction_marker_line;
public:
    data_parser();
    ~data_parser();

    void parse_data(std::string filepath);
    std::vector<markers_data> parse_line(const std::string& line);
    markers_data parse_line_to_point(const std::string& line);

    void add_marker_line(std::string filepath);
    void add_zero_line(std::string filepath);
    void add_direction_line(std::string filepath);

    void append(const std::vector<markers_data> &points_on_step);
    void append_marker_points(const std::vector<markers_data> &points);

    std::vector<marker_track_data> do_calc();
    marker_track_data do_calc_for_track(const std::vector<markers_data>& line_track);
};







#endif  /* DATA_PARSER_HPP */
