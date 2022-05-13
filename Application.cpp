#include <iostream>
#include "mpi.h"
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"


class Example : public olc::PixelGameEngine
{
public:
	Example()
	{
		sAppName = "Example";
	}

public:
	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		// called once per frame
		for (int x = 0; x < ScreenWidth(); x++)
			for (int y = 0; y < ScreenHeight(); y++)
				Draw(x, y, olc::Pixel(rand() % 255, rand() % 255, rand()% 255));	
		return true;
	}
};


int main(int argc, char** argv)
{
    int rank, size;

    MPI::Init(argc, argv);

    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    std::cout << "size: " << size << "\n";
    std::cout << "rank: " << rank << "\n";

    MPI::Finalize();

    Example demo;

    if (demo.Construct(256, 240, 4, 4))
		demo.Start();


    return 0;
}