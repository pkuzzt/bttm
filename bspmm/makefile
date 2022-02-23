CC=g++
FLAG= -Wall -fopenmp -g3

.PHONY: spmv_test

spmv_test: ./src/main.cpp ./build/matrix.o
	$(CC) $(FLAG) ./src/main.cpp ./build/matrix.o -o ./exe/spmv_test

./build/matrix.o: ./src/matrix.cpp ./src/matrix.hpp
	$(CC) $(FLAG) ./src/matrix.cpp -c -o ./build/matrix.o

clean: 
	rm -rf ./exe/* ./build/*


