#!/bin/sh
# track
gnuplot -e "set terminal png size 640,480; set output 'x_track.png'; plot 'track.dat' u 1:2 with lines"
gnuplot -e "set terminal png size 640,480; set output 'y_track.png'; plot 'track.dat' u 1:3 with lines"
# velocity
gnuplot -e "set terminal png size 640,480; set output 'x_velocity.png'; plot 'first_derivative.dat' u 1:2 with lines"
gnuplot -e "set terminal png size 640,480; set output 'y_velocity.png'; plot 'first_derivative.dat' u 1:3 with lines"
# acceleration
gnuplot -e "set terminal png size 640,480; set output 'x_accel.png'; plot 'second_derivative.dat' u 1:2 with lines"
gnuplot -e "set terminal png size 640,480; set output 'y_accel.png'; plot 'second_derivative.dat' u 1:3 with lines"
