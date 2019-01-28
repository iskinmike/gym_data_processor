#ifndef DATA_PROCESSOR_API_HPP
#define DATA_PROCESSOR_API_HPP

#include "data_parser.hpp"

class data_processor_api
{
	data_parser parser;
public:
        data_processor_api(){}
        ~data_processor_api(){}

	void setup_coordinates(double scale_length, const std::string& zero_point_path, const std::string& direction_point_path) {
	    parser.set_scale_mm_len(scale_length);
	    parser.append_zero_points(zero_point_path);
	    parser.append_direction_points(direction_point_path);
	    parser.eval_pixel_scale();
	}
	void append_point(const std::string& point_path) {
		parser.append_marker_points(point_path);
	}
	void evaluate(const std::string& directory_path){
		parser.dump_to_files(directory_path);
	}
};


#endif  /* DATA_PROCESSOR_API_HPP */
