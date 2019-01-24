

#include "data_parser.hpp"

#include <algorithm>
#include <set>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <regex>



struct delta_data{
    double delta;
    int line;
    int point;
    delta_data():delta(0), line(0), point(0){}
};



namespace {
double first_derivative(double p1, double p2, double delta){
    return (p2-p1)/delta;
}

derivative_2d_data calc_first_derivative(int pos, const std::vector<markers_data>& data){
    int left = (pos != 0) ? pos-1 : 0;
    int right = pos;
    double delta = (pos != 0) ? (data[right].step - data[left].step) : 1;
    derivative_2d_data der;
    der.x_coord_value = first_derivative(data[left].pos.x, data[right].pos.x, delta);
    der.y_coord_value = first_derivative(data[left].pos.y, data[right].pos.y, delta);
    der.step = (pos != 0) ? data[left].step : 1;
    der.pos = data[left].pos;
    return der;
}

derivative_2d_data calc_second_derivative(int pos, const std::vector<derivative_2d_data>& data){
    int left = (pos != 0) ? pos-1 : 0;
    int right = pos;
    double delta = (pos != 0) ? (data[right].step - data[left].step) : 1;
    derivative_2d_data der;
    der.x_coord_value = first_derivative(data[left].x_coord_value, data[right].x_coord_value,delta);
    der.y_coord_value = first_derivative(data[left].y_coord_value, data[right].y_coord_value, delta);
    der.step = data[left].step;
    der.pos = data[left].pos;
    return der;
}

marker_amplitude_data calc_amplitude_for_marker_from_derivative(const std::vector<derivative_2d_data>& derivative_track, coord coordinate){
    // считаем что у нас не теряются данные по точкам с производными.
    marker_amplitude_data amp_data;
    amp_data.coord_type = coordinate;

    double left_val = 0;
    double right_val = 0;
    int size_max = static_cast<int>(derivative_track.size());
    for (int i = 0; i < size_max; ++i) {
        // общий случай
        int left = i;
        int right = (size_max > i) ? i+1 : i;

        if (coord::y == coordinate) {
            left_val = derivative_track[left].y_coord_value;
            right_val = derivative_track[right].y_coord_value;
        } else {
            left_val = derivative_track[left].x_coord_value;
            right_val = derivative_track[right].x_coord_value;
        }

        if (left_val <= 0 && right_val > 0) { // смена с падения на рост роста или функция с нуля растет.
            amp_data.extremums.push_back(extremum_point(i, extremum_type::min));
            amp_data.minima_count++;
        }
        if (left_val >= 0 && right_val < 0) { // смена с роста на падение или с нуля падает
            amp_data.extremums.push_back(extremum_point(i, extremum_type::max));
            amp_data.maxima_count++;
        }
        if (left_val < 0 && right_val == 0) {
            amp_data.extremums.push_back(extremum_point(i, extremum_type::plateu));
            amp_data.maxima_count++;
        }
        if (left_val > 0.0f && right_val == 0.0f) {
//            amp_data.extremums.push_back(extremum_point(i, extremum_type::valley));
            amp_data.extremums.emplace_back(i, extremum_type::valley);
            amp_data.minima_count++;
        }

    }
    // Теперь мы готовы по порядку считать разницу в пикселях.
    for (int i=1; i< amp_data.extremums.size(); ++i) {
        shift_data shift;
        int left_extr = i-1;
        int right_extr = i;
        int left = amp_data.extremums[left_extr].pos;
        int right = amp_data.extremums[right_extr].pos;
        point left_marker = derivative_track[left].pos;
        point right_marker = derivative_track[right].pos;
        // Оставим пока то что мы берем вектор но надо по типу амплитуды разделять какую проекцию вектора мы хотим использовать
        // и выбирать соответствующую координатной оси
        shift.start_pos = derivative_track[left].step;
        shift.end_pos = derivative_track[right].step;
        shift.shift = right_marker.x - left_marker.x; //.dist(right_marker);
        shift.shift_vec = vector_2f(left_marker, right_marker);
//        std::cout << std::dec << "extr left pos " << derivative_track[left].step << "; extr right pos " << derivative_track[right].step << std::endl;
//        std::cout << std::dec << "shift left pos " << shift.start_pos << "; extr right pos " << shift.end_pos << std::endl;
//        std::cout << std::dec << "markers " << left_marker.print() << "; extr right pos " << right_marker.print() << std::endl;
//        std::cout << std::dec << "shift " << shift.shift/* << "; calced " << shift.shift_vec.get_len()*/ << std::endl;
        amp_data.shifts.push_back(shift);
    }

    return amp_data;
}

marker_amplitudes_xy calc_marker_amplitudes_xy(const std::vector<derivative_2d_data>& derivative_track){
    marker_amplitudes_xy amps;
    amps.amp_x = calc_amplitude_for_marker_from_derivative(derivative_track, coord::x);
    amps.amp_y = calc_amplitude_for_marker_from_derivative(derivative_track, coord::y);
    return amps;
}




}




data_parser::data_parser(){}

data_parser::~data_parser(){}

// Данные вида point_count step time x y x y x y ...
void data_parser::parse_data(std::string filepath) {
    std::fstream file(filepath, std::ios_base::in);
    if (file.is_open()) {
        std::string line;
        int counter = 0;
        while (std::getline(file, line)) {
//            std::cout << "line [" << line << "]" << std::endl;
            std::vector<markers_data> step_points = parse_line(line);
            append(step_points);
            counter++;
            if (counter > 50) break;
        }

        file.close();
    }
}

std::vector<markers_data> data_parser::parse_line(const std::string &line){
    const std::regex count_step_time("(\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+)");
    std::smatch match;
    std::vector<markers_data> result_points;

    if (std::regex_search(line, match, count_step_time)) {
//        std::cout << "full match " << match[0] << std::endl;
        int64_t count = std::atol(match[1].str().c_str());
        int64_t step = std::atol(match[2].str().c_str());
        int64_t time = std::atol(match[3].str().c_str());
//        std::cout << "cnt " << count
//                  << " step " << match[2]
//                  << " time " << time
//                  << std::endl;
        if (0 < count) {
            for (int i = 0; i < count; i++) {
                int pos_x = 2*(i+2);
                int pos_y = pos_x + 1;
                int64_t x = std::atol(match[pos_x].str().c_str());
                int64_t y = std::atol(match[pos_y].str().c_str());
                markers_data data;
                data.pos = point(x,y);
                data.step = step;
                data.time = time;
//                std::cout << "point " << x << "," <<y <<std::endl;
                result_points.push_back(data);
            }
        }
    } else {
        std::cout << "not found" << std::endl;
    }
    return result_points;
}

markers_data data_parser::parse_line_to_point(const std::string &line){
    const std::regex count_step_time("(\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+) (\\d+\\.?\\d+)");
    std::smatch match;
    markers_data result_point;

    if (std::regex_search(line, match, count_step_time)) {
//        int64_t count = std::atol(match[1].str().c_str());
        int64_t step = std::atol(match[2].str().c_str());
        int64_t time = std::atol(match[3].str().c_str());
        int64_t x = std::atol(match[4].str().c_str());
        int64_t y = std::atol(match[5].str().c_str());
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

void data_parser::add_marker_line(std::string filepath){
    std::fstream file(filepath, std::ios_base::in);
    if (file.is_open()) {
        std::string line;
        int counter = 0;
        std::vector<markers_data> points;
        while (std::getline(file, line)) {
//            std::vector<markers_data> step_points = parse_line_to_point(line);
            points.push_back(parse_line_to_point(line));
            counter++;
            if (counter > 50) break;
        }
        append_marker_points(points);

        file.close();
    }
}

void data_parser::append(const std::vector<markers_data>& points_on_step) {
    std::vector<markers_data> copy_points = points_on_step;
    // Если нет еще ни одной точки просто заполняем
    if (marker_line.empty()) {
        for (auto& el: points_on_step) {
            marker_line.emplace_back();
            marker_line.back().push_back(el);
        }
        return;
    }
    // Начинаем сравнивать с уже имеющимися на предыдущем шаге
    // сначала соберем сравниваемые точки.
    std::vector<markers_data> tmp;
    for (auto& el : marker_line) {
        tmp.push_back(el.back());
    }

    std::vector<delta_data> dldt;
    for (size_t i = 0; i < tmp.size(); ++i) {
        for (size_t j = 0; j < points_on_step.size(); ++j) {
            delta_data dat;
            dat.delta = tmp[i].pos.dist(points_on_step[j].pos);// calc_shift(tmp[i], points_on_step[j]);
            dat.line = i;
            dat.point = j;
            dldt.push_back(dat);
        }
    }

    std::set<int> points_to_setup;
    while (!dldt.empty()) {
        auto min_element = std::min_element(dldt.begin(), dldt.end(), [] (const delta_data &a, const delta_data &b){
            return a.delta < b.delta;
        });

        int line_to_setup = min_element->line;
        int point_to_setup = min_element->point;

        marker_line[line_to_setup].push_back(points_on_step[point_to_setup]);
//        std::cout << "[*****] marker_line append: [" << line_to_setup << "]" << std::endl;
        points_to_setup.insert(point_to_setup);
        dldt.erase(std::remove_if(dldt.begin(), dldt.end(),
                [&line_to_setup, &point_to_setup](const delta_data &a){
                    if (a.line == line_to_setup || a.point == point_to_setup) {
                        return true;
                    }
                    return false;
                }),
                dldt.end());
    }

    if (points_to_setup.size() != points_on_step.size()) {
        for (size_t i = 0; i < points_on_step.size(); ++i) {
            if (!points_to_setup.count(i)) {
                marker_line.emplace_back();
                marker_line.back().push_back(points_on_step[i]);
            }
        }
    }
}

void data_parser::append_marker_points(const std::vector<markers_data>& points){
    marker_line.emplace_back();
    marker_line.back() = points;
}

std::vector<marker_track_data> data_parser::do_calc(){
    std::vector<marker_track_data> result;
    for (const auto& el: marker_line) {
        result.push_back(do_calc_for_track(el));
    }
    return result;
}

marker_track_data data_parser::do_calc_for_track(const std::vector<markers_data>& line_track) {
    std::cout << "line_track size " << line_track.size() << std::endl;

    marker_track_data result;
//    result.first_derivative.reserve(line_track.size());
//    result.second_derivative.reserve(line_track.size());
//    result.track = line_track;
    for (int i = 0; i < line_track.size(); ++i) {
        derivative_2d_data first_der = calc_first_derivative(i, line_track);
        result.first_derivative.push_back(first_der);
    }

//    for (auto& el : result.first_derivative) {
//        std::cout << "dx/dt " << el.x_coord_value << " dy/dt " << el.y_coord_value << std::endl;
//    }

    for (int i = 0; i < result.first_derivative.size(); ++i) {
        derivative_2d_data second_der = calc_second_derivative(i, result.first_derivative);
        result.second_derivative.push_back(second_der);
    }

//    for (auto& el : result.second_derivative) {
//        std::cout << "ddx/ddt " << el.x_coord_value << " ddy/ddt " << el.y_coord_value << std::endl;
//    }


    std::cout << "Test minmax from derivatives." << std::endl;
    marker_amplitudes_xy amps = calc_marker_amplitudes_xy(result.first_derivative);


    marker_amplitude_data amp = amps.amp_y;
    std::cout << std::dec << "mins " << amp.minima_count << std::endl;
    std::cout << "maxes " << amp.maxima_count << std::endl;
    double mm_len = 63;
    double dist = 269;

    // pixel scale is mm_len/dist
    const double pix_scale = mm_len/dist;
    for (auto& el : amp.shifts) {
        std::cout << "shifts len[" << el.shift_vec.x_len << " | " << el.shift_vec.y_len << " | " << el.shift
                  << " pix;  start step [" << el.start_pos << "] end step [" << el.end_pos << "]"
                  << std::endl;
    }

    return result;
}






















