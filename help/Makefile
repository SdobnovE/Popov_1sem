a.out: main.o functions.o matrix_3diagonal.o
	g++ main.o functions.o matrix_3diagonal.o
main.o: main.cpp functions.h matrix_3diagonal.h
	g++ -c main.cpp
list.o: matrix_3diagonal.cpp matrix_3diagonal.h
	g++ -c matrix_3diagonal.cpp
student.o: functions.cpp functions.h
	g++ -c functions.cpp
