CC = g++
CFLAGS = -Wall -O3 --std=c++11

circle.o : circle.cpp
	mpicxx -c circle.cpp

square.o : square.cpp
	mpicxx -c square.cpp

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

all : square.o circle.o random.o genetic.o
	mpicxx random.o square.o genetic.o -o square.exe
	mpicxx random.o circle.o genetic.o -o circle.exe

run : circle.exe square.exe
	mpiexec -np 4 ./circle.exe
	mpiexec -np 4 ./square.exe
