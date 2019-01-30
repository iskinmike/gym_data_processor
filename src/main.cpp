#include <iostream>

#include <map>
//#include "data_parser.hpp"
#include "data_processor_api.hpp"

int main(int argc, char const *argv[])
{
    data_processor_api processor;
    std::string zero_path = "../../build/dist//zero_marker.dat";
    std::string direction_path = "../../build/dist//direction_marker.dat";
    std::string point_path = "../../build/dist//data.dat";
    std::string dump_dir = "result";
    processor.setup_coordinates(0.063f, zero_path, direction_path);
    processor.append_point(point_path);
    processor.evaluate(dump_dir);

	return 0;
}
