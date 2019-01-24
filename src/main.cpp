#include <iostream>

#include <map>
#include "data_parser.hpp"

int main(int argc, char const *argv[])
{
    data_parser parser;
    parser.add_marker_line("/home/mike/workspace/tmp/gym/data_processor/resources/data.dat");
    parser.do_calc();
	return 0;
}
