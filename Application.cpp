#include <iostream>
#include <omp.h>
#include <pthread.h>
#include "mpi.h"
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

#include <fstream>

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
    MPI::Init(argc, argv);


    int rank, size;
    char name[MPI_MAX_PROCESSOR_NAME];
    int namet;

    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    MPI::Get_processor_name(name, namet);
    std::cout << "Has " << size << " nodes in " << name << " computer, i'm node " << rank << "\n";

    std::ofstream *OutputFile;

    if(rank == 0)
    {
        OutputFile = new std::ofstream("out.ppm");
        *OutputFile << "P3\n" << "256 256\n" << "256\n";

        for (int x = 0; x < 256; x++)
			for (int y = 0; y < 256; y++)
				*OutputFile << (rand() % 255) << " " << (rand() % 255) << " " << (rand()% 255) << std::endl;	
		
    }

    if(rank == 0)
    {
        int b;
        MPI::COMM_WORLD.Recv((void*)&b, 1, MPI::INT, 1, MPI::ANY_TAG);
        std::cout << "I " << rank << " receive " << b << " from node 1!\n";
    }

    else
    {
        int a = 10;
        MPI::COMM_WORLD.Send((void*)&a, 1, MPI::INT, 0, 0);
        std::cout << "I " << rank << " have send " << a << " to node 0!\n";
    }

    MPI::Finalize();

    // TODO: ver se openmp esta mesmo funcionando
    //       tira o olcengine e verifica no gerenciador
    #pragma omp parallel
    {
        std::cout << "Num of threads " << omp_get_num_threads() << "\n";
        std::cout << "i'm thread " << omp_get_thread_num() << "\n";
    }

    if (rank == 0)
    {
        Example demo;
        if (demo.Construct(256, 256, 4, 4))
            demo.Start();
    }


    return 0;
}