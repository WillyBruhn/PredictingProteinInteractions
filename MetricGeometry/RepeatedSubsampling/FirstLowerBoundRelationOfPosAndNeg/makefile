COMP = g++ -std=c++11 -c -g -O3
LINK = g++ -std=c++11 -pthread -lboost_system

COMP_PROF = g++ -std=c++11 -c -g -O3 -pg
COMP_LINK = g++ -std=c++11 -pg -pthread -lboost_system

# -I /usr/include/boost/filesystem/

all: main
clean:
	rm -f *.o rm main main_profiled

main: main.o
	$(LINK) main.o -o main

main.o: main.cpp first_lower_bound.h ../mathFunctions.h ../random_points.h ../reading.h
	$(COMP) main.cpp -o main.o

main_profiled: main_profiled.o
	$(COMP_LINK) main_profiled.o -o main_profiled

main_profiled.o: main.cpp first_lower_bound.h ../mathFunctions.h ../random_points.h ../reading.h
	$(COMP_PROF) main.cpp -o main_profiled.o
