CMPLR=g++
CFLAGS=-o1.out -Wall -std=c++11
all:
	$(CMPLR) parsing.cpp Matrix.cpp Circuit.cpp main.cpp $(CFLAGS) 

clean:
	rm *.out
