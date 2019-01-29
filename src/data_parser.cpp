

#include "data_parser.hpp"

#include <algorithm>
#include <set>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <regex>
#include <sstream>

struct delta_data{
    double delta;
    int line;
    int point;
    delta_data():delta(0), line(0), point(0){}
};

void data_parser::set_scale_meter_lenght(double value) {
    scale_meter_len = value;
}



double data_parser::first_derivative(double p1, double p2, double delta){
//    std::cout << "first der: " << p1 << " | " << p2 << " | " << delta << std::endl;
    double delta_prc = (delta == 0) ? 1.0f : delta;
    return (double)(p2-p1)/delta_prc;
}

void data_parser::calc_median_velocity(shift_data &shift, const std::vector<derivative_2d_data>& derivative){
    double median_velocity_x = 0;
    double median_velocity_y = 0;
    int count = shift.end_pos - shift.start_pos;
    for (int i = shift.start_pos; i < shift.end_pos; ++i) {
        median_velocity_x += derivative[i].x_coord_value;
        median_velocity_y += derivative[i].y_coord_value;
    }
    shift.median_velocity_x = median_velocity_x/count;
    shift.median_velocity_y = median_velocity_y/count;
}

const double g = 9.8f;
double calc_work_y(double mass, double shift_len_projection, double acceleration, double dt){
    if (0 != shift_len_projection) {
        return mass*shift_len_projection*(acceleration + g);
    }
    return mass*g*g*dt*dt/2.0f;
}
double calc_work_x(double mass, double acceleration, double dt){
    return mass*acceleration*acceleration*dt*dt/2.0f;
}

void data_parser::calc_work_on_shift(shift_data &shift, const std::vector<derivative_2d_data> &second_derivative) {
    shift.work_y = 0;
    shift.work_x = 0;
    for (int i = shift.start_pos; i < shift.end_pos; ++i) {
        int left = i;
        int right = i+1;
//        double shift_projection_x = fabs(second_derivative[right].pos.x - second_derivative[left].pos.x);
        double shift_projection_y = fabs(second_derivative[right].pos.y - second_derivative[left].pos.y);
        double dt = second_derivative[right].time - second_derivative[left].time;
        shift.work_x += calc_work_x(1.0f, second_derivative[right].x_coord_value, dt);
        shift.work_y += calc_work_y(1.0f, shift_projection_y, second_derivative[right].y_coord_value, dt);
    }
}

void data_parser::dump_to_files(const std::string &dirpath){
    std::vector<marker_track_data> data = do_calc();
    int counter = 0;
    std::string sufix = ".dat";
    for (auto& el : data) {
        std::string num = std::to_string(counter++);
        dump_amplitudes_data(dirpath + "/shifts_x_" + num + sufix, el.amplitudes.amp_x);
        dump_amplitudes_data(dirpath + "/shifts_y_" + num + sufix, el.amplitudes.amp_y);
        dump_derivative_to_file(dirpath + "/first_derivative_" + num + sufix, el.first_derivative);
        dump_derivative_to_file(dirpath + "/second_derivative_" + num + sufix, el.second_derivative);
        dump_marker_track_to_file(dirpath + "/track_" + num + sufix, el.track);
    }
}

void data_parser::dump_amplitudes_data(const std::string &filepath, const marker_amplitude_data &data) {
    std::fstream file(filepath, std::ios_base::out);
    if (file.is_open()) {
        // Сдвиг сохраняем так: type x_shift, y_shift, start_time, end_time, median_velocity_x, median_velocity_y
        file << "type, x_shift, y_shift, start_time, end_time, duration, start_frame, end_frame, median_velocity_x, median_velocity_y, work_x, work_y\n";
        for (auto& el: data.shifts) {
            file << shift_data::get_text_type(el.type) << "," // type
                 << el.shift_vec.x_len << ","               // x shift
                 << el.shift_vec.y_len << ","               // y shift
                 << el.start_time << ","                    // start time
                 << el.end_time << ","                      // end time
                 << (el.end_time - el.start_time) << ","    // duration
                 << el.start_step << ","                    // start frame number
                 << el.end_step << ","                      // end frame number
                 << el.median_velocity_x << ","             // median velocity projection to x
                 << el.median_velocity_y << ","             // median velocity projection to y
                 << el.work_x << ","                        // work to move on x
                 << el.work_y << "\n";                      // work to move on y
        }
        std::cout << "file opened [" << filepath << "]" << std::endl;
    }
}

void data_parser::dump_derivative_to_file(const std::string &filepath, const std::vector<derivative_2d_data> &derivative){
    std::fstream file(filepath, std::ios_base::out);
    if (file.is_open()) {
        std::stringstream stream;
        for (auto& el: derivative) {
            stream << el.step << " " << el.time << " " << el.x_coord_value << " " << el.y_coord_value << "\n";
        }
        std::cout << "file opened [" << filepath << "]" << std::endl;
        file << stream.str();
    }
}

void data_parser::dump_marker_track_to_file(const std::string &filepath, const std::vector<markers_data> &track){
    std::fstream file(filepath, std::ios_base::out);
    if (file.is_open()) {
        std::stringstream stream;
        for (auto& el: track) {
            stream << el.step << " " << el.time << " " << el.pos.x << " " << el.pos.y << "\n";
        }
        file << stream.str();
    }
}


derivative_2d_data data_parser::calc_first_derivative(int pos, const std::vector<markers_data> &data){
    int left = (pos != 0) ? pos-1 : 0;
    int right = pos;
    double delta = (pos != 0) ? (data[right].time - data[left].time) : 1;
    derivative_2d_data der;
    der.x_coord_value = first_derivative(data[left].pos.x, data[right].pos.x, delta);
    der.y_coord_value = first_derivative(data[left].pos.y, data[right].pos.y, delta);
    der.time = data[right].time;
    der.pos = data[right].pos;
    der.step = data[right].step;
    return der;
}

derivative_2d_data data_parser::calc_second_derivative(int pos, const std::vector<derivative_2d_data> &data){
    int left = (pos != 0) ? pos-1 : 0;
    int right = pos;
    double delta = (pos != 0) ? (data[right].time - data[left].time) : 1;
    derivative_2d_data der;
    der.x_coord_value = first_derivative(data[left].x_coord_value, data[right].x_coord_value,delta);
    der.y_coord_value = first_derivative(data[left].y_coord_value, data[right].y_coord_value, delta);
    der.time = data[right].time;
    der.pos = data[right].pos;
    der.step = data[right].step;
    return der;
}

marker_amplitude_data data_parser::calc_amplitude_for_marker_from_derivative(const std::vector<derivative_2d_data> &derivative_track, const std::vector<derivative_2d_data> &second_derivative_track, coord coordinate){
    // считаем что у нас не теряются данные по точкам с производными.
    marker_amplitude_data amp_data;
    amp_data.coord_type = coordinate;

    double left_val = 0;
    double right_val = 0;
    int size_max = static_cast<int>(derivative_track.size());
    for (int i = 0; i < size_max; ++i) {
        // общий случай
        int left = i;
        int right = (size_max-1 > i) ? i+1 : i;

        if (coord::y == coordinate) {
            left_val = derivative_track[left].y_coord_value;
            right_val = derivative_track[right].y_coord_value;
        } else {
            left_val = derivative_track[left].x_coord_value;
            right_val = derivative_track[right].x_coord_value;
        }

        if (left_val <= 0 && right_val > 0) { // смена с падения на рост роста или функция с нуля растет.
            amp_data.extremums.emplace_back(i, extremum_type::min);
            amp_data.minima_count++;
        }
        if (left_val >= 0 && right_val < 0) { // смена с роста на падение или с нуля падает
            amp_data.extremums.emplace_back(i, extremum_type::max);
            amp_data.maxima_count++;
        }
        if (left_val < 0 && right_val == 0) {
            amp_data.extremums.emplace_back(i, extremum_type::valley);
            amp_data.valley_count++;
        }
        if (left_val > 0 && right_val == 0) {
            amp_data.extremums.emplace_back(i, extremum_type::plateau);
            amp_data.plateu_count++;
        }

    }
    std::cout << "calc_amplitude_for_marker_from_derivative size : " << derivative_track.size() << std::endl;
    // Теперь мы готовы по порядку считать разницу в пикселях.
    for (int i=1; i< amp_data.extremums.size(); ++i) {
        shift_data shift;
        int left_extr = i-1;
        int right_extr = i;
        int left = amp_data.extremums[left_extr].pos;
        int right = amp_data.extremums[right_extr].pos;
        extremum_type left_type = amp_data.extremums[left_extr].type;
        extremum_type right_type = amp_data.extremums[right_extr].type;
        point left_marker = derivative_track[left].pos;
        point right_marker = derivative_track[right].pos;
        // Оставим пока то что мы берем вектор но надо по типу амплитуды разделять какую проекцию вектора мы хотим использовать
        // и выбирать соответствующую координатной оси
        shift_type type = shift_type::plateau;
        if (left_type == extremum_type::max && right_type == extremum_type::min ||
            left_type == extremum_type::max && right_type == extremum_type::valley) {
            type = shift_type::down;
        } else if (left_type == extremum_type::min && right_type == extremum_type::max ||
                   left_type == extremum_type::min && right_type == extremum_type::plateau) {
            type = shift_type::up;
        }
        shift.start_pos = left;
        shift.end_pos = right;
        shift.start_step = derivative_track[left].step;
        shift.end_step = derivative_track[right].step;
        shift.start_time = derivative_track[left].time;
        shift.end_time = derivative_track[right].time;
        shift.shift = right_marker.x - left_marker.x;
        shift.shift_vec = vector_2f(left_marker, right_marker);
        shift.type = type;
        calc_median_velocity(shift, derivative_track);
        calc_work_on_shift(shift, second_derivative_track);
        amp_data.shifts.push_back(shift);
    }

    return amp_data;
}

marker_amplitudes_xy data_parser::calc_marker_amplitudes_xy(const std::vector<derivative_2d_data> &derivative_track, const std::vector<derivative_2d_data> &second_derivative_track){
    marker_amplitudes_xy amps;
    amps.amp_x = calc_amplitude_for_marker_from_derivative(derivative_track, second_derivative_track, coord::x);
    amps.amp_y = calc_amplitude_for_marker_from_derivative(derivative_track, second_derivative_track, coord::y);
    return amps;
}

data_parser::data_parser(): pixel_scale(1), scale_meter_len(100){}

data_parser::~data_parser(){}

// Данные вида point_count step time x y x y x y ...
markers_data data_parser::parse_line_to_point(const std::string &line){
    //    const std::regex count_step_time("(\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+)");
    const std::regex count_step_time("([\\d]+[\\.]?\\d*) ([\\d]+[\\.]?\\d*) ([\\d]+[\\.]?\\d*) ([\\d]+[\\.]?\\d*) ([\\d]+[\\.]?\\d*)");
    std::smatch match;
    markers_data result_point;

    if (std::regex_search(line, match, count_step_time)) {
        int match_pos = 2;
        int64_t step = std::atol(match[match_pos++].str().c_str());
        double time = std::atof(match[match_pos++].str().c_str());
        int64_t x = std::atol(match[match_pos++].str().c_str());
        int64_t y = std::atol(match[match_pos++].str().c_str());
        markers_data data;
        data.pos = point(x,y);
        data.step = step;
        data.time = time;
        return data;
    } else {
        std::cout << "not found" << std::endl;
    }
    return result_point;
}

std::vector<markers_data> data_parser::add_marker_line(const std::string& filepath, int limit){
    std::fstream file(filepath, std::ios_base::in);
    std::vector<markers_data> points;
    if (file.is_open()) {
        std::string line;
        int counter = 0;
        while (std::getline(file, line)) {
            points.push_back(parse_line_to_point(line));
            if (-1 != limit) {
                counter++;
                if (counter > limit) break;
            }
        }
        file.close();
    }
    return points;
}

void data_parser::append_marker_points(const std::string &filepath){
    marker_line.emplace_back();
    std::vector<markers_data> line = add_marker_line(filepath, -1);
    for (auto& el : line){
        el.pos.x *= pixel_scale;
        el.pos.y *= pixel_scale;
    }
    marker_line.back() = line;
}

void data_parser::eval_pixel_scale(){
    markers_data zero_point;
    markers_data direction_point;
    bool exit = false;
    for (auto& i : zero_marker_line) {
        for (auto& j : direction_marker_line) {
            if (i.step == j.step) {
                zero_point = i;
                direction_point = j;
                exit = true;
                break;
            }
        }
        if (exit) break;
    }

    if (exit) {
        double dist = zero_point.pos.dist(direction_point.pos);
        pixel_scale = scale_meter_len/dist;

    }
    std::cout << "exit: " << exit << std::endl;
    std::cout << "zero_point: " << zero_point.pos.print() << std::endl;
    std::cout << "direction_point: " << direction_point.pos.print() << std::endl;
    std::cout << "pixel scale: " << pixel_scale << std::endl;

}
void data_parser::append_zero_points(const std::string &filepath){
    zero_marker_line = add_marker_line(filepath, -1);
}
void data_parser::append_direction_points(const std::string &filepath){
    direction_marker_line = add_marker_line(filepath, -1);
}

std::vector<marker_track_data> data_parser::do_calc(){
    std::vector<marker_track_data> result;
    for (const auto& el: marker_line) {
        result.push_back(do_calc_for_track(el));
    }
    return result;
}

marker_track_data data_parser::do_calc_for_track(const std::vector<markers_data> &line_track) {
    marker_track_data result;
    result.first_derivative.clear();
    for (int i = 0; i < line_track.size(); ++i) {
        derivative_2d_data first_der = calc_first_derivative(i, line_track);
        result.first_derivative.push_back(first_der);
    }

    for (int i = 0; i < result.first_derivative.size(); ++i) {
        derivative_2d_data second_der = calc_second_derivative(i, result.first_derivative);
        result.second_derivative.push_back(second_der);
    }

    result.amplitudes = calc_marker_amplitudes_xy(result.first_derivative, result.second_derivative);

    marker_amplitude_data amp = result.amplitudes.amp_y;
    std::cout << std::dec << "mins " << amp.minima_count << std::endl;
    std::cout << "maxes " << amp.maxima_count << std::endl;
    std::cout << "valley " << amp.valley_count << std::endl;
    std::cout << "plateu " << amp.plateu_count << std::endl;

    for (auto& el : amp.shifts) {
        std::cout << "shifts len by Y [" << el.shift_vec.y_len << "] "
                  << " start step [" << el.start_pos << " - " << el.end_pos << "]"
                  << " velocity [" << el.median_velocity_y << "]"
                  << std::endl;
    }

    result.track = line_track;
    return result;
}






















