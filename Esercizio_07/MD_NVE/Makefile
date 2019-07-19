CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o random.o func.o clean.sh
	./clean.sh
	$(CC) random.o MolDyn_NVE.o func.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
func.o : func.cpp func.h
	$(CC) -c func.cpp -o func.o $(CFLAGS)
clean :
	rm *.o MolDyn_NVE.exe seed.out
