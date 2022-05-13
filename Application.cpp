#include <iostream>
#include <omp.h>
#include <pthread.h>
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
    char name[MPI_MAX_PROCESSOR_NAME];
    int namet;

    MPI::Init(argc, argv);

    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    MPI::COMM_WORLD.Get_name(name, namet);
    std::cout << "Has " << size << " nodes in " << name << " computer, i'm node " << rank << "\n";

    MPI::Finalize();

    #pragma omp parallel
    {
        std::cout << "Num of threads " << omp_get_num_threads() << "\n";
        std::cout << "i'm thread " << omp_get_thread_num() << "\n";
    }

    if (rank == 0)
    {
        Example demo;
        if (demo.Construct(256, 240, 4, 4))
            demo.Start();
    }


    return 0;
}