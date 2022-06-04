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

    int nWidth;
    int nHeight;

    int nMaxIteration = 32;
    int nMode = 0;

    int* pFractalIterations;

public:
	bool OnUserCreate() override
	{
		nWidth = ScreenWidth();
        nHeight = ScreenHeight();
        pFractalIterations = new int[nWidth * nHeight]{ 0 };
		return true;
	}

    void CreateFractal(const olc::vi2d& pixel_tl, const olc::vi2d& pixel_br, 
                       const olc::vd2d& frac_real, const olc::vd2d& frac_imag)
    {
        // std::cout << "====================================\n";
        // std::cout << "Pixel Top-Left: " << pixel_tl << "\nPixel Bottom-Right: " << pixel_br << "\n";
        // std::cout << "Frac Real: " << frac_real << "\nFrac Imaginary: " <<  frac_imag << "\n";

        #pragma omp parallel for
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
                pFractalIterations[(x * (int)pixel_br.x) + y] = n;
            }
        }
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

        // Select Mode
        if (GetKey(olc::Key::K0).bPressed) nMode = 0;
        if (GetKey(olc::Key::K1).bPressed) nMode = 1;
        // Modify the max iteration on the fly
        if (GetKey(olc::UP).bPressed) nMaxIteration += 32;
		if (GetKey(olc::DOWN).bPressed) nMaxIteration -= 32;
        if (nMaxIteration < 32) nMaxIteration = 32;

        // Cont the time with chrono clock
        auto tStart = std::chrono::high_resolution_clock::now();

        switch (nMode)
        {
            case 0: CreateFractal(pixel_tl, pixel_br,
                                  frac_real, frac_imag); break;
        }

        auto tEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> fTime = tEnd - tStart;

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
		}

        DrawString(0, 30, "Time Taken: " + std::to_string(fTime.count()) + "s", olc::WHITE, 3);
		DrawString(0, 60, "Iterations: " + std::to_string(nMaxIteration), olc::WHITE, 3);

		return true;
	}

protected:

    double map(double x, double oldLow, double oldHigh, double newLow, double newHigh)
    {
        double oldRange = (x - oldLow)/(oldHigh - oldLow);
        return oldRange * (newHigh - newLow) + newLow;
    }

// Pan and Zoom Created by the channel @OneLoneCoder
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
        // v.x = ((double)(n.x) / vScale.x) + vOffset.x;
        // v.y = ((double)(n.y) / vScale.y) + vOffset.y;

        // std::cout << "OFFSET: " << vOffset << "\nSCALE: " << vScale << "\n";

        olc::vd2d screen_tl_after;
        ScreenToWorld(screen_tl_before, screen_tl_after);
        olc::vd2d screen_br_after;
        ScreenToWorld(screen_br_before, screen_br_after);
        
        // std::cout << "+++++++++++++++++++++++++++++++++++++++++" << "\n";
        // std::cout << "Screen TL Before: " << screen_tl_before << "\n";
        // std::cout << "Screen BR Before: " << screen_br_before << "\n";
        // std::cout << "World TL Before: " << world_tl_before << "\n";
        // std::cout << "World BR Before: " << world_br_before << "\n";
        // std::cout << "Screen TL After: " << screen_tl_after << "\n";
        // std::cout << "Screen BR After: " << screen_br_after << "\n";
        
        world_real_after.x = map(screen_tl_after.x, (double)screen_tl_before.x, (double)screen_br_before.x, world_tl_before.x, world_br_before.x);
        world_real_after.y = map(screen_br_after.x, (double)screen_tl_before.x, (double)screen_br_before.x, world_tl_before.x, world_br_before.x);

        world_imag_after.x = map(screen_tl_after.y, (double)screen_tl_before.y, (double)screen_br_before.y, world_tl_before.y, world_br_before.y);
        world_imag_after.y = map(screen_br_after.y, (double)screen_tl_before.y, (double)screen_br_before.y, world_tl_before.y, world_br_before.y);

        // std::cout << "World Real After: " << world_real_after << "\n";
        // std::cout << "World Imag After: " << world_imag_after << "\n";

    }
};

int main(int argc, char** argv)
{
    MandelbrotFractal demo;
    if (demo.Construct(900, 900, 1, 1, false, false))
        demo.Start();

    return 0;
}
