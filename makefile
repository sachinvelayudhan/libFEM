eigenPath= $(CURDIR)/
incPath= $(CURDIR)/inc/
cst:
	g++ -std=c++17 -g -c mainLaplaceCST.cpp -I$(eigenPath) -I$(incPath)
	g++ -std=c++17 -o main.exe mainLaplaceCST.o
	./main.exe
bar:
	g++ -std=c++17 -g -c mainAxialBar.cpp -I$(eigenPath) -I$(incPath)
	g++ -std=c++17 -o main.exe mainAxialBar.o
	./main.exe
testMat:
	g++ -std=c++17 -g -c mainTestMatrixOperations.cpp -I$(eigenPath) -I$(incPath)
	g++ -std=c++17 -o main.exe mainTestMatrixOperations.o
	./main.exe

clean:
	rm -v *.o
	rm -v *.exe
