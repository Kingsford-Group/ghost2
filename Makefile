CC = g++
CFLAGS = -O3 -std=c++11
BOOSTDIR ?= /usr/local/lib
LFLAGS = -lboost_iostreams -lboost_system -lboost_thread

GHOST: GHOST.cpp
	$(CC) GHOST.cpp $(CFLAGS) -L $(BOOSTDIR) $(LFLAGS) -o GHOST

clean:
	\rm *.sig.gz *.sdf *.af

tar:
	tar cfv *.sig.gz *.sdf *.af
