CXX = mpic++
CXXFLAGS = -fopenmp
LFLAGS = -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17 -O3 -mavx2 
FILE = Application

all: Application.cpp
	$(CXX) $(CXXFLAGS) Application.cpp $(LFLAGS) -o app && make run
	
gen: Application.cpp
	$(CXX) $(CXXFLAGS) Application.cpp $(LFLAGS) -o app

run: app hosts
	mpirun --hostfile hosts ./app

clean:
	rm -rf *.o
