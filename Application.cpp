// Standart Stuff
#include <iostream>
#include <chrono>
#include <complex>

// Parallelization Stuff
#include <omp.h>
#include <pthread.h>
// #include "mpi.h"

// Drawing Stuff
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"


class MandelbrotFractal : public olc::PixelGameEngine
{
public:
	MandelbrotFractal()
	{
		sAppName = "Example";
	}

protected:
    int iWidth;
    int iHeight;
    int iMaxIteration;

public:
    float map(float x, float oldLow, float oldHigh, float newLow, float newHigh)
    {
        float oldRange = (x - oldLow)/(oldHigh - oldLow);
        return oldRange * (newHigh - newLow) + newLow;
    }

public:
	bool OnUserCreate() override
	{
		iWidth = ScreenWidth();
        iHeight = ScreenHeight();
        iMaxIteration = 100;
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		// called once per frame
		

        // #pragma omp parallel for
        for (int x = 0; x < iWidth; x++)
        {
            for (int y = 0; y < iHeight; y++)
            {
                float a = map(x, 0, iWidth, -2, 1);
                float b = map(y, 0, iHeight, -2, 1);

                int n = 0;

                float ca = a;
                float cb = b;

                while (n < iMaxIteration)
                {
                    // z1 = z0^2 + c
                    // z2 = c^2 + c
                    //      c^2 = a^2 - b^2 + 2abi

                    // C^2
                    float aa = a*a - b*b;
                    float bb = 2 * a * b;

                    // C^2 + C
                    a = aa + ca;
                    b = bb + cb;

                    // It diverges, or not...
                    if (a + b > 16)
                        break;

                    n++;
                }

                // How to Draw the Fractal?
                float bright = map(n, 0, 100, 0, 255);
                if (n == 100) bright = 0;

                Draw(x, y, olc::Pixel(bright, bright, bright, 255));	
            }
        }
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
        if (demo.Construct(600, 600, 1, 1))
            demo.Start();
    }


    return 0;
}