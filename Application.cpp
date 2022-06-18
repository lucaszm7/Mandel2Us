// Standart Stuff
#include <iostream>
#include <chrono>
#include <complex>

// Parallelization Stuff
#include <omp.h>
#include <pthread.h>
#include "mpi.h"

// Drawing Stuff
#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

// Paralelized sections
double map(double x, double oldLow, double oldHigh, double newLow, double newHigh)
{
    double oldRange = (x - oldLow)/(oldHigh - oldLow);
    return oldRange * (newHigh - newLow) + newLow;
}

void DivideFractal(double** pParam, 
                   const olc::vi2d&pixel_tl, const olc::vi2d& pixel_br, 
                   const olc::vd2d& frac_real, const olc::vd2d& frac_imag,
                   unsigned int nMaxIteration, unsigned int nNodesSize,
                   double exit_code = 0)
{
    olc::vi2d new_pixel_tl = pixel_tl;
    olc::vi2d new_pixel_br = { pixel_br.x / nNodesSize, pixel_br.y};

    double frac_real_part = std::abs(frac_real.y - frac_real.x) / nNodesSize;
    olc::vd2d new_frac_real = {frac_real.x, frac_real.x + frac_real_part};

    for(int i = 0; i < nNodesSize; ++i)
    {
        pParam[i][0] = new_pixel_tl.x;
        pParam[i][1] = new_pixel_tl.y;
        pParam[i][2] = new_pixel_br.x;
        pParam[i][3] = new_pixel_br.y;
        pParam[i][4] = new_frac_real.x;
        pParam[i][5] = new_frac_real.y;
        pParam[i][6] = frac_imag.x;
        pParam[i][7] = frac_imag.y;
        pParam[i][8] = nMaxIteration;
        pParam[i][9] = exit_code;

        new_pixel_tl.x = new_pixel_br.x;
        new_pixel_br.x = new_pixel_br.x + (pixel_br.x / nNodesSize);
        
        new_frac_real = {new_frac_real.y, new_frac_real.y + frac_real_part};
    }
}

void CreateFractal(const olc::vi2d& pixel_tl, const olc::vi2d& pixel_br, 
                       const olc::vd2d& frac_real, const olc::vd2d& frac_imag,
                       int* pFractalIterations, unsigned int nMaxIteration,
                       int nScreenHeightSize = 0)
{

    if (nScreenHeightSize == 0)
        nScreenHeightSize = pixel_br.x;
    
    auto CHUNK = (pixel_br.x - pixel_tl.x) / 16;

    #pragma omp parallel for schedule(dynamic, CHUNK) num_threads(omp_get_num_procs()) 
    for (int x = pixel_tl.x; x < pixel_br.x; x++)
    {
        for (int y = pixel_tl.y; y < pixel_br.y; y++)
        {
            double a = map(x, pixel_tl.x, pixel_br.x, frac_real.x, frac_real.y);
            double b = map(y, pixel_tl.y, pixel_br.y, frac_imag.x, frac_imag.y);

            int n = 0;

            double ca = a;
            double cb = b;

            while (n < nMaxIteration)
            {
                // z1 = z0^2 + c
                // z2 = c^2 + c
                //      c^2 = a^2 - b^2 + 2abi

                // C^2
                double aa = a*a - b*b;
                double bb = 2 * a * b;

                // C^2 + C
                a = aa + ca;
                b = bb + cb;

                // It diverges, or not...
                if (a + b > 16)
                    break;

                n++;
            }
            pFractalIterations[(x * nScreenHeightSize) + y] = n;
        }
    }
}

class MandelbrotFractal : public olc::PixelGameEngine
{
public:
	MandelbrotFractal()
	{
		sAppName = "Mandelbrot Fractal!!";
	}

protected:

    int nWidth;
    int nHeight;

    int nMaxIteration = 32;
    int nMode = 0;

    int* pFractalIterations;

    // MPI Coord stuff
    int nMyRank, nNodesSize;

    double** pNodesParam;

public:
	bool OnUserCreate() override
	{
		nWidth = ScreenWidth();
        nHeight = ScreenHeight();
        pFractalIterations = new int[ScreenWidth() * ScreenHeight()]{ 0 };

        // MPI
        nNodesSize = MPI::COMM_WORLD.Get_size();
        std::cout << "Nodes size: " << nNodesSize << "\n";
        nMyRank = MPI::COMM_WORLD.Get_rank();

        pNodesParam = new double*[nNodesSize];
        for(int i = 0; i < nNodesSize; ++i)
            pNodesParam[i] = new double[10]{-1};

		return true;
	}

    
	bool OnUserUpdate(float fElapsedTime) override
	{
		// Panning and Zoomig, credits to @OneLoneCoder who i'am inpired for
        olc::vd2d vMouse = {(double)GetMouseX(), (double)GetMouseY()};

        // Get the position of the mouse and move the world Final Pos - Inital Pos
        // This make us drag Around the Screen Space, with the OffSet variable
        if(GetMouse(0).bPressed)
        {
            vStartPan = vMouse;
        }

        if(GetMouse(0).bHeld)
        {
            vOffset -= (vMouse - vStartPan) / vScale;
            vStartPan = vMouse;
        }

        olc::vd2d vMouseBeforeZoom;
        ScreenToWorld(vMouse, vMouseBeforeZoom);

        if (GetKey(olc::Key::E).bHeld) vScale *= 1.1;
		if (GetKey(olc::Key::Q).bHeld) vScale *= 0.9;
		
		olc::vd2d vMouseAfterZoom;
		ScreenToWorld(vMouse, vMouseAfterZoom);
		vOffset += (vMouseBeforeZoom - vMouseAfterZoom);

        // Now we have a smaller screen, and want to map to the world coord
        
        olc::vi2d pixel_tl = {  0,  0 };
        olc::vi2d pixel_br = {  ScreenWidth(),  ScreenHeight() };

        olc::vd2d frac_tl = { -2.0, -1.0 };
        olc::vd2d frac_br = {  1.0,  1.0 };

        // Then in the limites we now have the cartesian coords we want to draw
        // cartesian plane starting at top-left to bottom-right

        olc::vd2d frac_real;
        olc::vd2d frac_imag;

        ScreenToFrac(pixel_tl, pixel_br, frac_tl, frac_br, frac_real, frac_imag);

        // Color Option
        if (GetKey(olc::Key::F1).bPressed) nColorMode = 0;
        if (GetKey(olc::Key::F2).bPressed) nColorMode = 1;
        if (GetKey(olc::Key::F3).bPressed) nColorMode = 2;
        if (GetKey(olc::Key::F4).bPressed) nColorMode = 3;
        if (GetKey(olc::Key::F5).bPressed) nColorMode = 4;
        // Modify the max iteration on the fly
        if (GetKey(olc::UP).bPressed) nMaxIteration += 32;
		if (GetKey(olc::DOWN).bPressed) nMaxIteration -= 32;
        if (nMaxIteration < 32) nMaxIteration = 32;

        // Divide Fractal
        DivideFractal(pNodesParam, pixel_tl, pixel_br, frac_real, frac_imag, nMaxIteration, nNodesSize);

        for(int i = 0; i < nNodesSize - 1; i++)
            MPI::COMM_WORLD.Send((void*)pNodesParam[i], 10, MPI::DOUBLE, i+1, 0);
        
        // Cont the time with chrono clock
        auto tStart = std::chrono::high_resolution_clock::now();

        CreateFractal({pNodesParam[nNodesSize-1][0], pNodesParam[nNodesSize-1][1]}, 
                      {pNodesParam[nNodesSize-1][2], pNodesParam[nNodesSize-1][3]}, 
                      {pNodesParam[nNodesSize-1][4], pNodesParam[nNodesSize-1][5]}, 
                      {pNodesParam[nNodesSize-1][6], pNodesParam[nNodesSize-1][7]}, 
                      pFractalIterations, pNodesParam[nNodesSize-1][8]);
        
        auto tEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> fTime = tEnd - tStart;

        for(int i = 0; i < nNodesSize - 1; i++)
        {
            MPI::COMM_WORLD.Recv((void*)(pFractalIterations + ((int)pNodesParam[i][0] * (int)pNodesParam[i][3])), ((ScreenWidth()*ScreenHeight()) / nNodesSize), MPI::INT, i + 1, MPI::ANY_TAG);
        }

        switch (nColorMode)
        {
            case 0: 
                // Render result to screen
                for (int x = 0; x < ScreenWidth(); x++)
                {
                    for (int y = 0; y < ScreenHeight(); y++)
                    {
                        int n = pFractalIterations[x * ScreenWidth() + y];
                        // Coloring Algorithm - Picked from https://solarianprogrammer.com/2013/02/28/mandelbrot-set-cpp-11/
                        double t = (double)n/(double)nMaxIteration;
                        // Use smooth polynomials for r, g, b
                        int rr = (int)(9*(1-t)*t*t*t*255);
                        int rg = (int)(15*(1-t)*(1-t)*t*t*255);
                        int rb =  (int)(8.5*(1-t)*(1-t)*(1-t)*t*255);
                        Draw(x, y, olc::Pixel(rr, rg, rb, 255));	
                    }
                } break;

            case 1:
                for (int x = 0; x < ScreenWidth(); x++)
                {
                    for (int y = 0; y < ScreenHeight(); y++)
                    {
                        int i = pFractalIterations[x * ScreenWidth() + y];
                        float n = (float)i;
                        float a = 0.1f;
                        // Thank you @Eriksonn - Wonderful Magic Fractal Oddball Man
                        Draw(x, y, olc::PixelF(0.5f * sin(a * n) + 0.5f, 0.5f * sin(a * n + 2.094f) + 0.5f,  0.5f * sin(a * n + 4.188f) + 0.5f));
                    }
                } break;
            case 2:
                for (int x = 0; x < ScreenWidth(); x++)
                {
                    for (int y = 0; y < ScreenHeight(); y++)
                    {
                        int n = pFractalIterations[x * ScreenWidth() + y];
                        Draw(x, y, olc::Pixel(n*n, n, n*3, 255));
                    } 
                } break;
        }

        DrawString(0, 30, "Time Taken: " + std::to_string(fTime.count()) + "s", olc::WHITE, 3);
		DrawString(0, 60, "Iterations: " + std::to_string(nMaxIteration), olc::WHITE, 3);
		DrawString(0, 120, "Draw Mode: F" + std::to_string(nColorMode + 1) + "/ F3", olc::WHITE, 3);

		return true;
	}

    ~MandelbrotFractal()
    {
        delete[] pFractalIterations;
        for(int i = 0; i < nNodesSize; ++i)
            delete[] pNodesParam[i];
        delete[] pNodesParam;
    }

// Pan and Zoom Created following tutorials on channel @OneLoneCoder
protected:
    // Pan & Zoom variables
	olc::vd2d vOffset = { 0.0, 0.0 };
	olc::vd2d vStartPan = { 0.0, 0.0 };
	olc::vd2d vScale = { 1.0, 1.0 };

    void ScreenToWorld(const olc::vi2d& n, olc::vd2d& v)
	{
		v.x = (double)(n.x) / vScale.x + vOffset.x;
		v.y = (double)(n.y) / vScale.y + vOffset.y;
	}

    // Converte coords from Screen Space to World Space
    void ScreenToFrac(const olc::vi2d& screen_tl_before, const olc::vi2d& screen_br_before, 
                           const olc::vd2d& world_tl_before, const olc::vd2d& world_br_before, 
                           olc::vd2d& world_real_after, olc::vd2d& world_imag_after)
    {
        olc::vd2d screen_tl_after;
        ScreenToWorld(screen_tl_before, screen_tl_after);
        olc::vd2d screen_br_after;
        ScreenToWorld(screen_br_before, screen_br_after);
        
        world_real_after.x = map(screen_tl_after.x, (double)screen_tl_before.x, (double)screen_br_before.x, world_tl_before.x, world_br_before.x);
        world_real_after.y = map(screen_br_after.x, (double)screen_tl_before.x, (double)screen_br_before.x, world_tl_before.x, world_br_before.x);

        world_imag_after.x = map(screen_tl_after.y, (double)screen_tl_before.y, (double)screen_br_before.y, world_tl_before.y, world_br_before.y);
        world_imag_after.y = map(screen_br_after.y, (double)screen_tl_before.y, (double)screen_br_before.y, world_tl_before.y, world_br_before.y);

    }
};

int main(int argc, char** argv)
{

    int nScreenWidth  = 600;
    int nScreenHeight = 600;

    MPI::Init(argc, argv);

    int nMyRank, nNodesSize;
    char sComputerName[MPI::MAX_PROCESSOR_NAME];
    int nComputerName;

    nNodesSize = MPI::COMM_WORLD.Get_size();
    nMyRank = MPI::COMM_WORLD.Get_rank();

    MPI::Get_processor_name(sComputerName, nComputerName);

    // Just master node create window
    if(nMyRank == 0)
    {
        MandelbrotFractal demo;
        if (demo.Construct(nScreenWidth, nScreenHeight, 1, 1, false, false))
            demo.Start();
        
        double** pFinishCode;
        pFinishCode = new double*[nNodesSize];
        for(int i = 0; i < (nNodesSize - 1); ++i)
        {
            pFinishCode[i] = new double[10];
            pFinishCode[i][9] = -1.0;
            MPI::COMM_WORLD.Send((void*)pFinishCode[i], 10, MPI::DOUBLE, i+1, 0);
            delete[] pFinishCode[i];
        }
        delete[] pFinishCode;

        std::cout << "I node " << nMyRank << " have finish!\n";
    }

    // Other nodes just do the computation
    else
    {
        double pParam[10]{0};
        int* pFractalIterations = new int[nScreenWidth * nScreenHeight]{0};
        
        while(pParam[9] >= 0)
        {
            // Receive
            MPI::COMM_WORLD.Recv((void*)pParam, 10, MPI::DOUBLE, 0, MPI::ANY_TAG);
            if(pParam[9] == -1)
                break;
            
            // Compute
            CreateFractal({pParam[0], pParam[1]}, {pParam[2], pParam[3]}, 
                          {pParam[4], pParam[5]}, {pParam[6], pParam[7]}, 
                          pFractalIterations, pParam[8], pParam[3]);
            
            // Send Back
            MPI::COMM_WORLD.Send((void*)(pFractalIterations + ((int)pParam[0] * (int)pParam[3])), ((nScreenHeight*nScreenWidth)/nNodesSize), MPI::INT, 0, 0);
        }

        delete[] pFractalIterations;
        std::cout << "I node " << nMyRank << " have finish!\n";
    }

    std::cout << "I node " << nMyRank << " am waiting for other nodes to finish...\n";
    MPI::COMM_WORLD.Barrier();
    MPI::Finalize();
    std::cout << "All nodes has finish!\n";

    return 0;
}
