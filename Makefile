GXX = g++
FLAG = -std=c++11 -O3 -pthread -g -w
INC = inc
SRC = src
OBJ = main.o GeneticAlgorithm.o Individual.o io.o
EXE = GACR


.PHONY: all
all: $(OBJ)
	$(GXX) $(FLAG) -I$(INC) $(OBJ) -o $(EXE)
	rm *.o

main.o: main.cpp
	$(GXX) $(FLAG) -I$(INC) -c main.cpp

GeneticAlgorithm.o: $(SRC)/GeneticAlgorithm.cpp $(INC)/GeneticAlgorithm.hpp
	$(GXX) $(FLAG) -I$(INC) -c $(SRC)/GeneticAlgorithm.cpp 

Individual.o: $(SRC)/Individual.cpp $(INC)/Individual.hpp
	$(GXX) $(FLAG) -I$(INC) -c $(SRC)/Individual.cpp

io.o: $(SRC)/io.cpp $(INC)/io.hpp
	$(GXX) $(FLAG) -I$(INC) -c $(SRC)/io.cpp

.PHONY: clean
clean:
	rm $(OBJ) $(EXE)