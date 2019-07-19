CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o func.o
	$(CC) random.o main.o func.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
func.o : func.cpp func.h
	$(CC) -c func.cpp -o func.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
