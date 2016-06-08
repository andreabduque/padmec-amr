all:
	g++ -o main.o -c -I/home/abd/interpolation main.cpp 
	g++ -o Edge.o -c -I/home/abd/interpolation Edge.cpp
	g++ -o Vector.o -c -I/home/abd/interpolation Vector.cpp	
	g++ -o Interpolation.o -c -I/home/abd/interpolation Interpolation.cpp
	g++ -o main main.o Edge.o Vector.o Interpolation.o
