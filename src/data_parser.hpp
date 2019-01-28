#ifndef DATA_PARSER_HPP
#define DATA_PARSER_HPP


#include <string>
#include <vector>

#include "utils.hpp"


struct point_data{
    point pos;
    int64_t step;
    double time;
    std::vector<float> matrix;
    point_data() : pos(0,0){}
    point_data(const point_data& other) {
        pos = other.pos;
        step = other.step;
        time = other.time;
        matrix = other.matrix;
    }
};

enum class extremum_type {
    min, max, plateau, valley
};

enum class coord {
    x,y
};

enum class shift_type {
    down, up, plateau
};

struct shift_data{
    shift_type type;
    vector_2f shift_vec;
    double shift;
    double start_pos;
    double end_pos;
    int start_step;
    int end_step;
    double start_time;
    double end_time;
    double median_velocity_x;
    double median_velocity_y;
    static std::string get_text_type(shift_type type) {
        switch (type) {
        case shift_type::down:
            return std::string{"down"};
        case shift_type::up:
            return std::string{"up"};
        case shift_type::plateau:
            return std::string{"plateau"};
        }
    }
};

struct markers_data {
    point pos;
    double time;
    int64_t step;
};


struct extremum_point {
    extremum_type type;
    int pos;
    extremum_point(int pos, extremum_type type) : type(type), pos(pos) {}
};

struct marker_amplitude_data {
    coord coord_type;
    int minima_count;
    int maxima_count;
    int plateu_count;
    int valley_count;
    std::vector<extremum_point> extremums;
    std::vector<shift_data> shifts;
    marker_amplitude_data(): coord_type(coord::x),
        minima_count(0), maxima_count(0),
        plateu_count(0), valley_count(0) {}
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
};


class data_parser
{
    std::vector<std::vector<markers_data>> marker_line;
    std::vector<markers_data> zero_marker_line;
    std::vector<markers_data> direction_marker_line;
    double pixel_scale;
    double scale_mm_len;
public:
    data_parser();
    ~data_parser();

    markers_data parse_line_to_point(const std::string& line);

    std::vector<markers_data> add_marker_line(const std::string& filepath, int limit);
    void append_zero_points(const std::string& filepath);
    void append_direction_points(const std::string& filepath);
    void append_marker_points(const std::string& filepath);

    void eval_pixel_scale();

    std::vector<marker_track_data> do_calc();
    marker_track_data do_calc_for_track(const std::vector<markers_data>& line_track);
    void set_scale_mm_len(double value);


    derivative_2d_data calc_first_derivative(int pos, const std::vector<markers_data>& data);
    derivative_2d_data calc_second_derivative(int pos, const std::vector<derivative_2d_data>& data);
    marker_amplitude_data calc_amplitude_for_marker_from_derivative(const std::vector<derivative_2d_data>& derivative_track, coord coordinate);
    marker_amplitudes_xy calc_marker_amplitudes_xy(const std::vector<derivative_2d_data>& derivative_track);
    double first_derivative(double p1, double p2, double delta);

    void calc_median_velocity(shift_data& shift, const std::vector<derivative_2d_data> &derivative);

    void dump_to_files(const std::string& dirpath);
    void dump_amplitudes_data(const std::string& filepath, const marker_amplitude_data& data);
    void dump_derivative_to_file(const std::string& filepath, const std::vector<derivative_2d_data>& derivative);
    void dump_marker_track_to_file(const std::string& filepath, const std::vector<markers_data>& track);
};







#endif  /* DATA_PARSER_HPP */
