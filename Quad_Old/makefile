eigenPath= $(CURDIR)/
solverIntef= $(CURDIR)/Solver_Interfaces/
cst:
	g++ -std=c++17 -g -c testCST.cpp -I$(eigenPath) -I$(solverIntef)
	g++ -std=c++17 -o main.exe testCST.o
	./main.exe
matrix_test_1:
	g++ -std=c++17 -g -c */main_MatrixDemo_1.cpp
	g++ -std=c++17 -o main.exe main_MatrixDemo_1.o
quad:
	g++ -std=c++17 -g -c testQUAD.cpp -I$(eigenPath) -I$(solverIntef)
	g++ -std=c++17 -o main.exe testQUAD.o
	./main.exe
clean:
	rm -v *.o
	rm -v *.exe
