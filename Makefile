CXX = g++
CXXFLAGS = -fopenmp
LFLAGS = -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17 -O3 -mavx2 
FILE = Application

all: $(FILE).cpp
	$(CXX) $(CXXFLAGS) $(FILE).cpp $(LFLAGS) -o app && make run
	
gen: $(FILE).cpp
	$(CXX) $(CXXFLAGS) $(FILE).cpp $(LFLAGS) -o app

run: app hosts
	mpirun --hostfile hosts ./app

install:
	sudo apt-get install build-essential libglu1-mesa-dev libpng-dev

clean:
	rm -rf *.o
