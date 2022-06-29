eigenPath= $(CURDIR)/
solverIntef= $(CURDIR)/Solver_Interfaces/
incPath= $(CURDIR)/inc/
quad:
	g++ -std=c++17 -g -c mainLaplaceQuad.cpp -I$(eigenPath) -I$(solverIntef) -I$(incPath)
	g++ -std=c++17 -o main.exe mainLaplaceQuad.o
	./main.exe
quad1:
	g++ -std=c++17 -g -c mainLaplaceQuad1.cpp -I$(eigenPath) -I$(solverIntef) -I$(incPath)
	g++ -std=c++17 -o main.exe mainLaplaceQuad1.o
	./main.exe
GraFEA:
	g++ -std=c++17 -g -c mainGraFEA2D.cpp
	g++ -std=c++17 -o main.exe mainGraFEA2D.o
	./main.exe
clean:
	rm -v *.o
	rm -v *.exe
