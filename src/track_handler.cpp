
#include "track_handler.hpp"
#include <sys/time.h>
#include <sys/resource.h>
#include <errno.h>

namespace {
    int gen_origin_id(){
        static int orig_id = 1;
        return orig_id++;
    }

}

void track_handler::append_origin(std::vector<point_data> origin_points){
    // Работаем только как минимум с 2-мя точками.
    if (2 > origin_points.size()) return;
    if (origin_line.empty()) {
        if (2 == origin_points.size()){
            point_data p1 = origin_points[0];
            point_data p2 = origin_points[1];
            coord_system_origin origin;
            origin.first_marker.frame_time = p1.time;
            origin.first_marker.pos = p1.pos;
            origin.first_marker.id = gen_origin_id();
            origin.first_marker.matrix = p1.matrix;
            origin.first_marker.empty = false;
            origin.second_marker.frame_time = p2.time;
            origin.second_marker.pos = p2.pos;
            origin.second_marker.id = gen_origin_id();
            origin.second_marker.matrix = p2.matrix;
            origin.second_marker.empty = false;
            origin_line.push_back(std::move(origin));
        }
        return;
    }
    // теперь надо разлисить этои точки.. а хотя так ли это нужно нам все равно важно знать расстояние между ними..
    // хотя впоследствии нужно будет узнать  направление осей и тогда потребуется так что дучше сделать сейчас.
    int pos_p1 = 0;
    int pos_p2 = 1;

    auto& firts_marker = origin_line.back().first_marker;
    auto& second_marker = origin_line.back().second_marker;

    double dist_id1_to_point_1 = firts_marker.pos.dist(origin_points[pos_p1].pos);
    double dist_id1_to_point_2 = firts_marker.pos.dist(origin_points[pos_p2].pos);

    if (dist_id1_to_point_2 < dist_id1_to_point_1) {
        pos_p1 = 1;
        pos_p2 = 0;
    }

    coord_system_origin origin;
    origin.first_marker.frame_time = origin_points[pos_p1].time;
    origin.first_marker.pos = origin_points[pos_p1].pos;
    origin.first_marker.id = firts_marker.id;
    origin.first_marker.matrix = origin_points[pos_p1].matrix;
    origin.first_marker.empty = false;
    origin.second_marker.frame_time = origin_points[pos_p2].time;
    origin.second_marker.pos = origin_points[pos_p2].pos;
    origin.second_marker.id = second_marker.id;
    origin.second_marker.matrix = origin_points[pos_p2].matrix;
    origin.second_marker.empty = false;
    origin_line.push_back(std::move(origin));
}

void track_handler::append(std::vector<point_data> points_on_step) {
    std::vector<point_data> copy_points = points_on_step;
    // Если нет еще ни одной точки просто заполняем
    if (plot.empty()) {
        for (auto& el: points_on_step) {
            plot.emplace_back();
            plot.back().push_back(el);
        }
        return;
    }
    // Начинаем сравнивать с уже имеющимися на предыдущем шаге
    // сначала соберем сравниваемые точки.
    std::vector<point_data> tmp;
    for (auto& el : plot) {
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

        plot[line_to_setup].push_back(points_on_step[point_to_setup]);
        std::cout << "[*****] plot append: [" << line_to_setup << "]" << std::endl;
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
                plot.emplace_back();
                plot.back().push_back(points_on_step[i]);
            }
        }
    }
}

size_t track_handler::get_lines_count(){
    return plot.size();
}

std::vector<point_data> track_handler::get_line(int n){
    return plot[n];
}

void track_handler::dump_all_to_file(const std::string &path){
    int count = 0;
    for (auto& line: plot) {
        std::stringstream stream;
//        stream << "step, time, x, y\n";
        for (auto& el: line) {
            stream << el.step << " " /*<< el.time*/ << " " << el.pos.x << " " << el.pos.y << "\n";
        }
        std::string fpath = path + std::to_string(count++) + ".plot";
        std::fstream file(fpath, std::ios_base::out);
        if (file.is_open()) {
            file << stream.str();
            file.close();
        }
    }
}

int track_handler::get_right_pos(int i, const std::vector<point_data> &data){
    int max_i = data.size()-1;
    if (i >= max_i) return max_i;
    int right_pos = i+1;
    while (data[i].pos.y == data[right_pos].pos.y) {
        ++right_pos;
        if (right_pos == max_i) {
//            return -1;
            break;
        }
    }
    return right_pos;
}


void track_handler::calc_amp(){
    // Пока работаем с 3-мя нашими точками.
    // определим какие 2 из 3-х являются системой координат.
    int coord_id_1 = 0;
    bool is_1_setted = false;
    int coord_id_2 = 0;
    int proc_id = 0;
    int pos = 0;
    for (auto& el : plot) {
        if (el.begin()->pos.x > 300) {
            proc_id = pos;
        } else if (!is_1_setted){
            coord_id_1 = pos;
            is_1_setted = true;
        } else {
            coord_id_2 = pos;
        }
        pos++;
    }

    // Берем процесс.
    // ищем все минимумы и максимумы. по порядку как они друг за другом идут.
    // Определяем между ними время и т.д.
    std::vector<int> min_poses;
    std::vector<int> max_poses;

    auto& process = plot[proc_id];
    for (int i = 1; i < process.size()-1; ++i) {
        int left = i-1;
        int right = get_right_pos(i, process);
        std::cout << std::dec << "right pos " << right << std::endl;
        if (-1 == right) break;
        // find minimums
        if (process[i].pos.y < process[right].pos.y && process[i].pos.y < process[left].pos.y) {
            min_poses.push_back(i);
        } else {
            // find maximums
            if (process[i].pos.y > process[right].pos.y && process[i].pos.y > process[left].pos.y) {
                max_poses.push_back(i);
            }
        }
        i = right-1;
    }

    std::cout << "mins " << min_poses.size() << std::endl;
    std::cout << "maxes " << max_poses.size() << std::endl;
    std::cout << "process " <<  process.size() << std::endl;
    std::cout << "proc_id " << proc_id << std::endl;

    for(auto& pos : max_poses) {
        std::cout << "max " << process[pos].step << " " << process[pos].pos.y << std::endl;
    }

    // тепер после того как нашли min и max надо сгруппировать их.
    // затем вычислить размер перемещения в пикселях
    std::vector<point_data> mm_points;
    std::vector<int> min_maxes;
    min_maxes.insert(min_maxes.end(), min_poses.begin(), min_poses.end());
    min_maxes.insert(min_maxes.end(), max_poses.begin(), max_poses.end());
    std::sort(min_maxes.begin(),min_maxes.end());

    for (auto& el : min_maxes) {
        mm_points.push_back(process[el]);
    }

    // Теперь мы готовы по порядку считать разницу в пикселях.
    std::vector<shift_data> shifts;
    for (int i=1; i< mm_points.size(); ++i) {
        shift_data shift;
        if (mm_points[0].pos.y < mm_points[1].pos.y) { // начинаем с минимума
            shift.shift = mm_points[i-1].pos.y - mm_points[i].pos.y;
        } else { // начинаем с максимума
            shift.shift = mm_points[i].pos.y - mm_points[i-1].pos.y;
        }
        shifts.push_back(shift);
    }
    for(auto& el : shifts) {
        std::cout << "shifts: [" << el.shift << "]" << std::endl;
    }

    // Известна разница в пикселях.
    // Теперь нужно это в миллиметры перевести.

    // сначала найдем расстояние между метками Системы координат (СК)
    // берем данные из первого кадра.
    point cs1(plot[coord_id_1].begin()->pos.x, plot[coord_id_1].begin()->pos.y);
    point cs2(plot[coord_id_2].begin()->pos.x, plot[coord_id_2].begin()->pos.y);

    // посчитаем погрешность определения координат таким образом проверим алгоритм.
    std::cout << "*** calc artool detection error" << std::endl;
    double sum = 0;
    double mean = 0;
    // считаем среднее
    std::vector<double> dists;
    std::vector<point_data>& coord_1 = plot[coord_id_1];
    std::vector<point_data>& coord_2 = plot[coord_id_2];

    double coord_1_sum = 0;
    double coord_2_sum = 0;
    for (auto& i : coord_1) {
        coord_1_sum += i.pos.x;
    }
    for (auto& i : coord_2) {
        coord_2_sum += i.pos.x;
    }
    double coord_1_mean = coord_1_sum/coord_1.size();
    double coord_2_mean = coord_2_sum/coord_2.size();

    std::vector<double> errors;
    for (auto& i : coord_1) {
        errors.push_back(fabs(i.pos.x-coord_1_mean));
    }
    double coord_1_error = *std::max_element(errors.begin(), errors.end());
    errors.clear();
    for (auto& i : coord_2) {
        errors.push_back(fabs(i.pos.x-coord_2_mean));
    }
    double coord_2_error = *std::max_element(errors.begin(), errors.end());


    for (auto& i : coord_1) {
        for (auto& j : coord_2) {
            if (i.step == j.step) {
                point p1(i.pos.x, i.pos.y);
                point p2(j.pos.x, j.pos.y);
                sum += p1.dist(p2);
                dists.push_back(p1.dist(p2));
                break;
            }
        }
    }

    mean = sum/dists.size();
    std::cout << "*** dists.size() " << dists.size() << std::endl;
    std::cout << "*** sum " << sum << std::endl;
    std::cout << "*** mean " << mean << std::endl;
    std::vector<double> deltas;
//    deltas.reserve(dists.size());
    for (int i = 0; i < dists.size(); ++i) {
        deltas.push_back(fabs(dists[i] - mean));
//        std::cout << "dists " << dists[i] << std::endl;
    }

//    for (auto& el : deltas) {
//        std::cout << "delta " << el << std::endl;
//    }

    auto max = std::max_element(deltas.begin(), deltas.end());
    double artool_error = *max;
    std::cout << "ARTool error: [" << artool_error << "]" << std::endl;


    std::cout << "coord_1_error: [" << coord_1_error << "]" << std::endl;
    std::cout << "coord_2_error: [" << coord_2_error << "]" << std::endl;


    double dist = cs1.dist(cs2);
    const double mm_len = 63;

    double rel_mm_len_error = 1/mm_len;
    double rel_pix_dist_error = 6/dist;

    // pixel scale is mm_len/dist
    const double pix_scale = mm_len/dist;
    std::cout << "cs dist: [" << dist << "]" << std::endl;
    std::cout << "shift scale: [" << pix_scale << "]" << std::endl;
    std::cout << "rel_mm_len_error: [" << rel_mm_len_error << "]" << std::endl;
    std::cout << "rel_pix_dist_error: [" << rel_pix_dist_error << "]" << std::endl;
    std::cout << "rel error: [" << rel_mm_len_error + rel_pix_dist_error << "]" << std::endl;
    double clced_mm_error = mm_len*0.04;//(rel_pix_dist_error + rel_mm_len_error);
    std::cout << "clced_mm_error: [" << clced_mm_error << "]" << std::endl;
    for(auto& el : shifts) {
        std::cout << "shifts(mm): [" << el.shift*pix_scale << "] +- [" << clced_mm_error << "]" << std::endl;
    }
}

marker_amplitude_data track_handler::calc_amplitude_for_marker(const std::vector<point_data>& marker_track, coord coordinate){
    // Берем процесс.
    // ищем все минимумы и максимумы. по порядку как они друг за другом идут.
    // Определяем между ними время и т.д.

    // TODO: рефакторим таким образом чтобы надо было только 1 раз эту функцию вызвать и уже сразу все получить
    marker_amplitude_data amp_data;
    amp_data.coord_type = coordinate;
    double left_val = 0;
    double middle_val = 0;
    double right_val = 0;

//    bool mid_less_left = false;
//    bool mid_less_right = false;
//    int last_mid = 0;
    for (int i = 0; i < marker_track.size(); ++i) {
//        if () {

//        }


        int left = i-1;
        int right = get_right_pos(i, marker_track);


        if (0 == i) {
            i = right-1;
            left = i;
        }

        if (marker_track.size() == right) {
//            last_mid = i;
            continue;
        }


        if (coord::y == coordinate) {
            left_val = marker_track[left].pos.y;
            middle_val = marker_track[i].pos.y;
            right_val = marker_track[right].pos.y;
//            mid_less_left = marker_track[i].pos.y < marker_track[left].pos.y;
//            mid_less_right = marker_track[i].pos.y < marker_track[right].pos.y;
        } else {
            left_val = marker_track[left].pos.x;
            middle_val = marker_track[i].pos.x;
            right_val = marker_track[right].pos.x;
//            mid_less_left = marker_track[i].pos.x < marker_track[left].pos.x;
//            mid_less_right = marker_track[i].pos.x < marker_track[right].pos.x;
        }

        if (0 == i) {

        } else if (marker_track.size() == i || right){

        }


        // find minimums
        if (middle_val < right_val && middle_val < left_val) {
            amp_data.extremums.push_back(extremum_point(i, extremum_type::min));
            amp_data.minima_count++;
        } else {
            // find maximums
            if (middle_val > right_val && middle_val > left_val) {
                amp_data.extremums.push_back(extremum_point(i, extremum_type::max));
                amp_data.maxima_count++;
            }
        }
        i = right-1;
    }

    std::cout << "mins " << amp_data.minima_count << std::endl;
    std::cout << "maxes " << amp_data.maxima_count << std::endl;
    std::cout << "marker_track " <<  marker_track.size() << std::endl;

    // Теперь мы готовы по порядку считать разницу в пикселях.
    for (int i=1; i< amp_data.extremums.size(); ++i) {
        shift_data shift;
        int left_extr = i-1;
        int right_extr = i;
        int left = amp_data.extremums[left_extr].pos;
        int right = amp_data.extremums[right_extr].pos;
        point_data left_marker = marker_track[left];
        point_data right_marker = marker_track[right];
        // Оставим пока то что мы берем вектор но надо по типу амплитуды разделять какую проекцию вектора мы хотим использовать
        // и выбирать соответствующую координатной оси
//        shift.shift_vec = vector_2f(left_marker.pos, right_marker.pos);
        amp_data.shifts.push_back(shift);
    }

    return amp_data;
}

marker_amplitude_data track_handler::calc_amplitude_for_marker_from_derivative(const std::vector<derivative_2d_data> &derivative_track, const std::vector<point_data> &marker_track, coord coordinate){
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
        if (left_val > 0 && right_val == 0) {
            amp_data.extremums.push_back(extremum_point(i, extremum_type::valley));
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
        point left_marker = marker_track[left].pos;
        point right_marker = marker_track[right].pos;
        // Оставим пока то что мы берем вектор но надо по типу амплитуды разделять какую проекцию вектора мы хотим использовать
        // и выбирать соответствующую координатной оси
        shift.start_pos = marker_track[left].step;
        shift.end_pos = marker_track[right].step;
        shift.shift = right_marker.y - left_marker.y; //.dist(right_marker);
//        shift.shift_vec = vector_2f(left_marker, right_marker);
        std::cout << std::dec << "extr left pos " << marker_track[left].step << "; extr right pos " << marker_track[right].step << std::endl;
        std::cout << std::dec << "shift left pos " << shift.start_pos << "; extr right pos " << shift.end_pos << std::endl;
        std::cout << std::dec << "markers " << left_marker.print() << "; extr right pos " << right_marker.print() << std::endl;
        std::cout << std::dec << "shift " << shift.shift/* << "; calced " << shift.shift_vec.get_len()*/ << std::endl;
        amp_data.shifts.push_back(shift);
    }

    return amp_data;
}

marker_amplitudes_xy track_handler::calc_marker_amplitudes_xy(const std::vector<point_data> &marker_track){
    marker_amplitudes_xy amps;
    amps.amp_x = calc_amplitude_for_marker(marker_track, coord::x);
    amps.amp_y = calc_amplitude_for_marker(marker_track, coord::y);
    return amps;
}

namespace {

double first_derivative(double p1, double p2, double delta){
    return (p2-p1)/delta;
}

derivative_2d_data calc_first_derivative(int pos, const std::vector<point_data>& data){
    int left = (pos != 0) ? pos-1 : 0;
    int right = pos;
    double delta = data[right].step - data[left].step;
    derivative_2d_data der;
    der.x_coord_value = first_derivative(data[left].pos.x, data[right].pos.x, delta);
    der.y_coord_value = first_derivative(data[left].pos.y, data[right].pos.y, delta);
    der.step = data[left].step;
//    der.pos = data[left].pos;
    return der;
}

double second_derivative(double p1, double p2, double delta){
    return (p2-p1)/delta;
}

derivative_2d_data calc_second_derivative(int pos, const std::vector<derivative_2d_data>& data){
    int left = pos;
    int right = (data.size() > pos) ? pos+1 : pos;
    double delta = data[right].step - data[left].step;
    derivative_2d_data der;
    der.x_coord_value = second_derivative(data[left].x_coord_value, data[right].x_coord_value,delta);
    der.y_coord_value = second_derivative(data[left].y_coord_value, data[right].y_coord_value, delta);
    der.step = data[left].step;
    return der;
}

void dump_derivative_to_file(const std::vector<derivative_2d_data>& derivative, std::string filepath){
    std::fstream file(filepath, std::ios_base::out);
    if (file.is_open()) {
        std::stringstream stream;
//    stream << "step, x, y\n";
        for (auto& el: derivative) {
            stream << el.step << " " << el.x_coord_value << " " << el.y_coord_value << "\n";
//            std::cout  << el.step << " " << el.x_coord_value << " " << el.y_coord_value << "\n";
        }
        std::cout << "file opened [" << filepath << "]" << std::endl;
        file << stream.str();
//        file.close();
    }
}

double func(double t){
    return pow(t,2);
    return 25*t;
    return sin(t);
}

} // namespace

void track_handler::calc_derivatives(){
    // Пока работаем с 3-мя нашими точками.
    // определим какие 2 из 3-х являются системой координат.
//    int coord_id_1 = 0;
//    bool is_1_setted = false;
//    int coord_id_2 = 0;
    int proc_id = 0;
//    int pos = 0;
//    for (auto& el : plot) {
//        if (el.begin()->pos.x > 300) {
//            proc_id = pos;
//        } else if (!is_1_setted){
//            coord_id_1 = pos;
//            is_1_setted = true;
//        } else {
//            coord_id_2 = pos;
//        }
//        pos++;
//    }
    // Берем процесс, ищем производные
    std::vector<derivative_2d_data> first_derivative;
    std::vector<derivative_2d_data> second_derivative;

    std::cout << std::dec << "--=== plot size. [" << plot.size() << "]" << std::endl;
    std::vector<point_data> process = plot.back();
//    std::vector<point_data> process;
//    process.reserve(plot.size());
//    for (auto el : plot.back()) {
//        process.push_back(el);
//    }
    std::cout << "--=== process size. [" << process.size() << "]" << std::endl;
//    std::vector<point_data> process = plot[proc_id];
    for (int i = 0; i < process.size(); ++i) {
        derivative_2d_data first_der = calc_first_derivative(i, process);
//        std::cout  << first_der.step << " " << first_der.x_coord_value << " " << first_der.y_coord_value << "\n";
        first_derivative.push_back(first_der);
//        first_derivative.emplace_back(first_der);
    }
    std::cout << "--=== Firt DERIVATIVE VECTOR. [" << first_derivative.size() << "]" << std::endl;
    for (auto& el: first_derivative){
//        std::cout  << el.step << " " << el.x_coord_value << " " << el.y_coord_value << "\n";
    }

    for (int i = 0; i < first_derivative.size(); ++i) {
        derivative_2d_data second_der = calc_second_derivative(i, first_derivative);
        second_derivative.push_back(second_der);
    }


    // Добавить дамп в файл и проверить на известном процессе.
    // dump to file
    dump_derivative_to_file(first_derivative, "first_derivative.dat");
    dump_derivative_to_file(second_derivative, "second_derivative.dat");

    std::cout << "Test minmax from derivatives." << std::endl;
    marker_amplitude_data amp = calc_amplitude_for_marker_from_derivative(first_derivative, process, coord::y);
    std::cout << std::dec << "mins " << amp.minima_count << std::endl;
    std::cout << "maxes " << amp.maxima_count << std::endl;
    double mm_len = 63;
    double dist = 269;
//    double rel_mm_len_error = 1/mm_len;
//    double rel_pix_dist_error = 6/dist;

    // pixel scale is mm_len/dist
    const double pix_scale = mm_len/dist;
    for (auto& el : amp.shifts) {
        std::cout << "shifts len[" << pix_scale << " | " << el.shift
//                  << "] begin[" << el.shift_vec.get_begin().print()
//                  << "] end[" << el.shift_vec.get_end().print() << "]"
                  << " start step [" << el.start_pos << "] end step [" << el.end_pos << "]"
                  << std::endl;
    }


    ///// Take test calc
    // create test points
//    std::vector<point_data> sin_points;
//    for (int i = 0; i< 200; ++i) {
//        point_data dat;
//        dat.step = i+1;
//        dat.pos.x = func(i*0.1);
//        dat.pos.y = func(i*0.1);
//        sin_points.push_back(dat);
//    }

//    // calc derivatives
//    first_derivative.clear();
//    second_derivative.clear();
//    for (int i = 0; i < sin_points.size()-1; ++i) {
//        derivative_2d_data first_der = calc_first_derivative(i, sin_points);
//        first_derivative.push_back(first_der);
//    }

//    for (int i = 0; i < first_derivative.size()-1; ++i) {
//        derivative_2d_data second_der = calc_second_derivative(i, first_derivative);
//        second_derivative.push_back(second_der);
//    }
//    dump_derivative_to_file(first_derivative, "sin_first_derivative.dat");
//    dump_derivative_to_file(second_derivative, "sin_second_derivative.dat");


}




//point vector_2f::get_end() const
//{
//    return end;
//}

//double vector_2f::get_len() const
//{
//    return len;
//}

//double vector_2f::get_x_len() const
//{
//    return x_len;
//}

//double vector_2f::get_y_len() const
//{
//    return y_len;
//}

//point vector_2f::get_begin() const
//{
//    return begin;
//}

point &point::operator =(point &&other) {
    x = other.x;
    y = other.y;
    other.x = 0;
    other.y = 0;
    return *this;
}

double point::dist(point other){
    double delta;
    delta = sqrt(pow(x - other.x, 2) + pow(y - other.y, 2));
    return fabs(delta);
}

std::string point::print(){
    std::stringstream str;
    str << "(" << x << ", " << y << ")";
    return str.str();
}

point point::operator +(const point &other) const { return point(x+other.x, y+other.y); }

point point::operator -(const point &other) const { return point(x-other.x, y-other.y); }

point point::operator *(double d) const { return point(x*d, y*d); }

const point &point::operator =(const point &other) {
    x = other.x;
    y = other.y;
    return *this;
}

point::point(point&& other) {
    x = other.x;
    y = other.y;
    other.x = 0;
    other.y = 0;
}

point::point(const point &other) {
    x = other.x;
    y = other.y;
}
