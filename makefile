eigenPath= $(CURDIR)/
solverIntef= $(CURDIR)/Solver_Interfaces/
incPath= $(CURDIR)/inc/
quad:
	g++ -std=c++17 -g -c mainLaplaceQuad.cpp -I$(eigenPath) -I$(solverIntef) -I$(incPath)
	g++ -std=c++17 -o main.exe mainLaplaceQuad.o
clean:
	rm -v *.o
	rm -v *.exe
