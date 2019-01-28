#include <iostream>

#include <map>
//#include "data_parser.hpp"
#include "data_processor_api.hpp"

int main(int argc, char const *argv[])
{
    data_processor_api processor;
    std::string zero_path = "/home/mike/workspace/tmp/gym/error_test_processor/resources/zero_marker.dat";
    std::string direction_path = "/home/mike/workspace/tmp/gym/error_test_processor/resources/direction_marker.dat";
    std::string point_path = "/home/mike/workspace/tmp/gym/error_test_processor/resources/data.dat";
    std::string dump_dir = "/home/mike/workspace/tmp/gym/error_test_processor/build/result/";
    processor.setup_coordinates(63, zero_path, direction_path);
    processor.append_point(point_path);
    processor.evaluate(dump_dir);

//    parser.set_scale_mm_len(63);
//    parser.append_zero_points("/home/mike/workspace/tmp/gym/data_processor/resources/zero_marker.dat");
//    parser.append_direction_points("/home/mike/workspace/tmp/gym/data_processor/resources/direction_marker.dat");
//    parser.eval_pixel_scale();
//    processor.append_marker_points("/home/mike/workspace/tmp/gym/data_processor/resources/data.dat");
//    std::vector<marker_track_data> data = parser.do_calc();
//    processor.dump_to_files("/home/mike/workspace/tmp/gym/data_processor/build/result/");
	return 0;
}
