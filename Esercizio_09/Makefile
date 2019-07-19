CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o genetic.o
	$(CC) random.o main.o genetic.o -o main.exe

main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
