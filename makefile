main:
	g++ -std=c++17 -g -c main.cpp
	g++ -std=c++17 -o main.exe main.o
	./main.exe
cst:
	g++ -std=c++17 -g -c testCST.cpp
	g++ -std=c++17 -o main.exe testCST.o
	./main.exe
matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o

clean:
	rm -v *.o
	rm -v *.exe
