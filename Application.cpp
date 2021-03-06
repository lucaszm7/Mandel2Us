// Standart Stuff
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <fstream>
#include <cstdlib>

// Parallelization Stuff
#include <omp.h>

#if defined _WIN32
        #include <intrin.h>
        namespace MPI
        {
            void Init(int &, char **& ){}
            void Finalize( void ){}
            struct nan 
            {
                int Get_size( void ) const { return 1; }
                int Get_rank( void ) const { return 0; }
                void Send( const void * v1, int v2, double v3, int v4, int v5 ) const {}
                void Send( const void * v1, int v2,  int v3, int v4, int v5 ) const {}
                void Send( const void * v1, int v2,  float v3, int v4, int v5 ) const {}
                void Recv( void * v1, int v2, int v3, int v4, int v5 ) const {}
                void Barrier( void ) const {}
            };
            struct Datatype {};
            nan COMM_WORLD;
            double DOUBLE;
            int INT;
            float FLOAT;
            int ANY_TAG;
        }
        int setenv (const char *__name, const char *__value, int __replace) { return 0; }
#endif

#ifndef _WIN32
    #include <pthread.h>
    #include "mpi.h"
#endif


bool constexpr UseMPI = true;

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
    if(nNodesSize == 1)
    {
        pParam[0][0] = pixel_tl.x;
        pParam[0][1] = pixel_tl.y;
        pParam[0][2] = pixel_br.x;
        pParam[0][3] = pixel_br.y;
        pParam[0][4] = frac_real.x;
        pParam[0][5] = frac_real.y;
        pParam[0][6] = frac_imag.x;
        pParam[0][7] = frac_imag.y;
        pParam[0][8] = nMaxIteration;
        pParam[0][9] = exit_code;
    }

    else
    {
        olc::vi2d new_pixel_tl = pixel_tl;
        olc::vi2d new_pixel_br = { pixel_br.x / (int)nNodesSize, pixel_br.y};

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
}

void CreateFractalSequential(const olc::vi2d& pixel_tl, const olc::vi2d& pixel_br, 
                       const olc::vd2d& frac_real, const olc::vd2d& frac_imag,
                       int* pFractalIterations, unsigned int nMaxIteration,
                       int nScreenHeightSize = 0)
{
    
    nScreenHeightSize = pixel_br.y;
    for (int x = pixel_tl.x; x < (pixel_br.x - 1); x++)
    {
        for (int y = pixel_tl.y; y < pixel_br.y; y++)
        {
            double a = map(x, pixel_tl.x, pixel_br.x, frac_real.x, frac_real.y);
            double b = map(y, pixel_tl.y, pixel_br.y, frac_imag.x, frac_imag.y);

            int n = 0;

            double ca = a;
            double cb = b;

            while (n < nMaxIteration && (a*a + b*b) < 4 )
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

                n++;
            }
            pFractalIterations[(x * nScreenHeightSize) + y] = n;
        }
    }
}

void CreateFractalParallel(const olc::vi2d& pixel_tl, const olc::vi2d& pixel_br, 
                       const olc::vd2d& frac_real, const olc::vd2d& frac_imag,
                       int* pFractalIterations, unsigned int nMaxIteration,
                       int nScreenHeightSize = 0)
{
    
    nScreenHeightSize = pixel_br.y;

    auto CHUNK = (pixel_br.x - pixel_tl.x) / 128;
    #pragma omp parallel for schedule(dynamic, CHUNK) num_threads(omp_get_num_procs()) 
    for (int x = pixel_tl.x; x < (pixel_br.x - 1); x++)
    {
        for (int y = pixel_tl.y; y < pixel_br.y; y++)
        {
            double a = map(x, pixel_tl.x, pixel_br.x, frac_real.x, frac_real.y);
            double b = map(y, pixel_tl.y, pixel_br.y, frac_imag.x, frac_imag.y);

            int n = 0;

            double ca = a;
            double cb = b;

            while (n < nMaxIteration && (a*a + b*b) < 4 )
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

                n++;
            }
            pFractalIterations[(x * nScreenHeightSize) + y] = n;
        }
    }
}

void CreateFractalParallelAVX(const olc::vi2d& pixel_tl, const olc::vi2d& pixel_br, 
                       const olc::vd2d& frac_real, const olc::vd2d& frac_imag,
                       int* pFractalIterations, int nMaxIteration,
                       int nScreenHeightSize = 0)
{
    
    nScreenHeightSize = pixel_br.y;
    
    double x_scale = (frac_real.y - frac_real.x) / (double(pixel_br.x) - double(pixel_tl.x));
	double y_scale = (frac_imag.y - frac_imag.x) / (double(pixel_br.y) - double(pixel_tl.y));

    double x_pos = frac_real.x;

    // 64-bit "double" registers
    __m256d _aa, _bb, _ca, _cb, _a, _b, _zr2, _zi2, _two, _four, _mask1;

    // 64-bit "int" registers
    __m256i _n, _maxIt, _mask2, _c, _one;

    // start of Y
    __m256d _y_pos_offsets, _y_pos, _y_scale, _y_jump;

    _y_scale = _mm256_set1_pd(y_scale);
    _y_jump = _mm256_set1_pd(y_scale * 4);
    _y_pos_offsets = _mm256_set_pd(0, 1, 2, 3);
    _y_pos_offsets = _mm256_mul_pd(_y_pos_offsets, _y_scale);

    // | 32.0 | 32.0 | 32.0 | 32.0 | 
    _maxIt = _mm256_set1_epi64x(nMaxIteration);

    // | 1.0 | 1.0 | 1.0 | 1.0 | 
    _one = _mm256_set1_epi64x(1);

    // | 2.0 | 2.0 | 2.0 | 2.0 | 
    _two = _mm256_set1_pd(2.0);

    // | 4.0 | 4.0 | 4.0 | 4.0 | 
    _four = _mm256_set1_pd(4.0);

    auto numProcs = omp_get_num_procs();

    auto CHUNK = (pixel_br.x - pixel_tl.x) / 128;
    #pragma omp parallel for schedule(dynamic, CHUNK) num_threads(numProcs) private(_n, _y_pos, x_pos, _c, _ca, _cb, _a, _b, _aa, _bb, _zr2, _zi2, _mask1, _mask2)
    for (int x = pixel_tl.x; x < (pixel_br.x - 1); x++)
    {
        // Calc start x
        x_pos = (frac_real.x + ((x) * x_scale));

        // Reset y position
        _bb =  _mm256_set1_pd(frac_imag.x);
        _y_pos = _mm256_add_pd(_bb, _y_pos_offsets);

        _ca = _mm256_set1_pd(x_pos);

        for (int y = pixel_tl.y; y < pixel_br.y; y += 4)
        {

            _a = _mm256_setzero_pd();
            _b = _mm256_setzero_pd();

            _n = _mm256_setzero_si256();

            _cb = _y_pos;

            repeat:

            // double ca = a;
            // double cb = b;

            // double aa = a*a - b*b;
            // double bb = 2 * a * b;

            // a = aa + ca;
            // b = bb + cb;

            // Multiply 256-bit registers in parallel, as they are doubles

            // a * a
            _zr2 = _mm256_mul_pd(_a, _a); // a * a

            // b * b
            _zi2 = _mm256_mul_pd(_b, _b); // b * b

            // a*a - b*b
            _aa = _mm256_sub_pd(_zr2, _zi2); // (a * a) - (b * a)

            // a * b
            _bb = _mm256_mul_pd(_a, _b); // a * b

            // (bb * 2)
            _b = _mm256_mul_pd(_bb, _two); // ((a * b) * 2)

            // (bb) + cb
            _b = _mm256_add_pd(_b, _cb); // ((a * b) * 2) + cb

            // aa + ca
            _a = _mm256_add_pd(_aa, _ca); // ((a * a) - (b * b)) + ca


            // while ((zr2 + zi2) < 4.0 && n < nMaxIteration)

            // aa = (a * a + b * b)
            _aa = _mm256_add_pd(_zr2, _zi2);

            // m1 = if(aa < 4.0)
            _mask1 = _mm256_cmp_pd(_aa, _four, _CMP_LT_OQ);

            // m2 = (nMaxIteration > n)
            _mask2 = _mm256_cmpgt_epi64(_maxIt, _n);

            // m2 = m1 AND m2 = if(aa < 4.0 && nMaxIterations > n)
            _mask2 = _mm256_and_si256(_mask2, _mm256_castpd_si256(_mask1));

            // mask2 AND 00...01
            // mask2 = |00...00|11...11|00...00|11...11|
            // one   = |00...01|00...01|00...01|00...01|
            // c     = |00...00|00...01|00...00|00...00| // just the 2 element has to be incremented
            _c = _mm256_and_si256(_mask2, _one);

            // n + c
            _n = _mm256_add_epi64(_n, _c);

            // if ((a * a + b * b) < 4.0 && n < nMaxIterations) goto repeat
            if (_mm256_movemask_pd(_mm256_castsi256_pd(_mask2)) > 0)
                goto repeat;

            #if defined _WIN32
                pFractalIterations[(x * nScreenHeightSize) + y + 0] = int(_n.m256i_i64[3]);
                pFractalIterations[(x * nScreenHeightSize) + y + 1] = int(_n.m256i_i64[1]);
                pFractalIterations[(x * nScreenHeightSize) + y + 2] = int(_n.m256i_i64[2]);
                pFractalIterations[(x * nScreenHeightSize) + y + 3] = int(_n.m256i_i64[0]);
            #endif
            #ifndef _WIN32
                pFractalIterations[(x * nScreenHeightSize) + y + 0] = int(_n[3]);
                pFractalIterations[(x * nScreenHeightSize) + y + 1] = int(_n[1]);
                pFractalIterations[(x * nScreenHeightSize) + y + 2] = int(_n[2]);
                pFractalIterations[(x * nScreenHeightSize) + y + 3] = int(_n[0]);
            #endif
            _y_pos = _mm256_add_pd(_y_pos, _y_jump);
        }
        // x_pos += x_scale;
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
    int nColorMode = 0;
    int nFracMode = 0;

    int* pFractalIterations;

    // MPI Coord stuff
    int nMyRank, nNodesSize;

    double** pNodesParam;

// UI variables
protected:

    bool toggleBenchmark = true;
    bool toggleHelp = true;
    bool toggleScreenShotView = false;
    bool toggleScreenShot = false;
    int nScreenShotCount = 0;

    std::string calcName;
    olc::Pixel color;
    olc::Pixel fracModeColor[3];
    olc::Pixel fracColorCol[3];
    olc::Pixel fracCol;

    olc::Pixel colorIterations = olc::WHITE;
    olc::Pixel zoomColor[2];
    olc::Pixel mouseColor;

public:
	bool OnUserCreate() override
	{
        char* ssc = getenv("SCREENSHOTCOUNT");
        if(ssc)
            nScreenShotCount = atoi(ssc);

		nWidth = ScreenWidth();
        nHeight = ScreenHeight();
        pFractalIterations = new int[ScreenWidth() * ScreenHeight()]{ 0 };

        // MPI
        if(!UseMPI)
        {
            nNodesSize = 1;
            nMyRank = 0;
        }
        else
        {
            nNodesSize = MPI::COMM_WORLD.Get_size();
            std::cout << "Nodes size: " << nNodesSize << "\n";
            nMyRank = MPI::COMM_WORLD.Get_rank();
        }

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

        mouseColor = olc::WHITE;
        if(GetMouse(0).bHeld)
        {
            vOffset -= (vMouse - vStartPan) / vScale;
            vStartPan = vMouse;
            mouseColor = olc::RED;
        }

        olc::vd2d vMouseBeforeZoom;
        ScreenToWorld(vMouse, vMouseBeforeZoom);

        zoomColor[0] = olc::WHITE;
        zoomColor[1] = olc::WHITE;
        if (GetKey(olc::Key::E).bHeld) { vScale *= 1.1; zoomColor[0] = olc::RED; }
		if (GetKey(olc::Key::Q).bHeld) { vScale *= 0.9; zoomColor[1] = olc::RED; }
		
		olc::vd2d vMouseAfterZoom;
		ScreenToWorld(vMouse, vMouseAfterZoom);
		vOffset += (vMouseBeforeZoom - vMouseAfterZoom);

        // Now we have a smaller screen, and want to map to the world coord
        
        olc::vi2d pixel_tl = {  0,  0 };
        olc::vi2d pixel_br = {  nWidth,  nHeight };

        olc::vd2d frac_tl = { -2.0, -1.0 };
        olc::vd2d frac_br = {  1.0,  1.0 };

        // Then in the limites we now have the cartesian coords we want to draw
        // cartesian plane starting at top-left to bottom-right

        olc::vd2d frac_real;
        olc::vd2d frac_imag;

        ScreenToFrac(pixel_tl, pixel_br, frac_tl, frac_br, frac_real, frac_imag);

        //Toggles
        if(GetKey(olc::Key::H).bPressed) { toggleHelp = !toggleHelp; }
        if(GetKey(olc::Key::I).bPressed) { toggleBenchmark = !toggleBenchmark; }
        if(GetKey(olc::Key::S).bPressed) { toggleScreenShotView = !toggleScreenShotView; }
        if(GetKey(olc::Key::P).bPressed) { toggleScreenShot = !toggleScreenShot; }
        

        // Calculation Option
        if (GetKey(olc::Key::K1).bPressed) { nFracMode = 0; }
        if (GetKey(olc::Key::K2).bPressed) { nFracMode = 1; }
        if (GetKey(olc::Key::K3).bPressed) { nFracMode = 2; }
        // Color Option
        if (GetKey(olc::Key::F1).bPressed) nColorMode = 0;
        if (GetKey(olc::Key::F2).bPressed) nColorMode = 1;
        if (GetKey(olc::Key::F3).bPressed) nColorMode = 2;
        // Modify the max iteration on the fly
        if (GetKey(olc::UP).bPressed)   { nMaxIteration += 32; colorIterations.g -= 16; colorIterations.b -= 16; }
		if (GetKey(olc::DOWN).bPressed) { nMaxIteration -= 32; colorIterations.g += 16; colorIterations.b += 16; }
        if (nMaxIteration < 32) { nMaxIteration = 32; colorIterations.g = 255; colorIterations.b = 255; }

        // Divide Fractal
        DivideFractal(pNodesParam, pixel_tl, pixel_br, frac_real, frac_imag, nMaxIteration, nNodesSize);

        // Cont the time with chrono clock
        auto tStart = std::chrono::high_resolution_clock::now();

        for(int i = 0; i < nNodesSize - 1; i++)
        {
            MPI::COMM_WORLD.Send((void*)pNodesParam[i], 10, MPI::DOUBLE, i+1, 0);
        }
        
        switch (nFracMode)
        {
            case 0:
                CreateFractalSequential({(int)pNodesParam[nNodesSize-1][0], (int)pNodesParam[nNodesSize-1][1]}, 
                            {(int)pNodesParam[nNodesSize-1][2], (int)pNodesParam[nNodesSize-1][3]}, 
                            {pNodesParam[nNodesSize-1][4], pNodesParam[nNodesSize-1][5]}, 
                            {pNodesParam[nNodesSize-1][6], pNodesParam[nNodesSize-1][7]}, 
                            pFractalIterations, (unsigned int)pNodesParam[nNodesSize-1][8]); break;
            case 1:
                CreateFractalParallel({(int)pNodesParam[nNodesSize-1][0], (int)pNodesParam[nNodesSize-1][1]}, 
                            {(int)pNodesParam[nNodesSize-1][2], (int)pNodesParam[nNodesSize-1][3]}, 
                            {pNodesParam[nNodesSize-1][4], pNodesParam[nNodesSize-1][5]}, 
                            {pNodesParam[nNodesSize-1][6], pNodesParam[nNodesSize-1][7]}, 
                            pFractalIterations, (unsigned int)pNodesParam[nNodesSize-1][8]); break;
            case 2:
                CreateFractalParallelAVX({(int)pNodesParam[nNodesSize-1][0], (int)pNodesParam[nNodesSize-1][1]}, 
                            {(int)pNodesParam[nNodesSize-1][2], (int)pNodesParam[nNodesSize-1][3]}, 
                            {pNodesParam[nNodesSize-1][4], pNodesParam[nNodesSize-1][5]}, 
                            {pNodesParam[nNodesSize-1][6], pNodesParam[nNodesSize-1][7]}, 
                            pFractalIterations, (unsigned int)pNodesParam[nNodesSize-1][8]); break;

        }

        for(int i = 0; i < nNodesSize - 1; i++)
        {
            MPI::COMM_WORLD.Recv((void*)(pFractalIterations + ((int)pNodesParam[i][0] * (int)pNodesParam[i][3])), ((ScreenWidth()*ScreenHeight()) / nNodesSize), MPI::INT, i + 1, MPI::ANY_TAG);
        }

        auto tEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> fTime = tEnd - tStart;

        switch (nColorMode)
        {
            case 0: 
                // Render result to screen
                if(toggleScreenShot)
                {
                    std::ofstream screenShot("screen_shot_" + std::to_string(nScreenShotCount++) +".ppm");
                    screenShot << "P3\n" << std::to_string(ScreenHeight()) + " " + std::to_string(ScreenWidth()) + "\n" << "256\n";
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            // Coloring Algorithm - Picked from https://solarianprogrammer.com/2013/02/28/mandelbrot-set-cpp-11/
                            double t = (double)n/(double)nMaxIteration;
                            // Use smooth polynomials for r, g, b
                            int rr = (int)(9*(1-t)*t*t*t*255);
                            int rg = (int)(15*(1-t)*(1-t)*t*t*255);
                            int rb =  (int)(8.5*(1-t)*(1-t)*(1-t)*t*255);
                            screenShot << rr << " " << rg << " " << rb << std::endl;
                        }
                    }
                    screenShot.close();
                    toggleScreenShot = false;
                }
                else
                {         
                    #pragma omp parallel for
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            // Coloring Algorithm - Picked from https://solarianprogrammer.com/2013/02/28/mandelbrot-set-cpp-11/
                            double t = (double)n/(double)nMaxIteration;
                            // Use smooth polynomials for r, g, b
                            int rr = (int)(9*(1-t)*t*t*t*255);
                            int rg = (int)(15*(1-t)*(1-t)*t*t*255);
                            int rb =  (int)(8.5*(1-t)*(1-t)*(1-t)*t*255);
                            Draw(x, y, olc::Pixel(rr, rg, rb, 255));
                        }
                    }
                }
                break;

            case 1:
                if(toggleScreenShot)
                {
                    std::ofstream screenShot("screen_shot_" + std::to_string(nScreenShotCount++) +".ppm");
                    screenShot << "P3\n" << std::to_string(ScreenHeight()) + " " + std::to_string(ScreenWidth()) + "\n" << "256\n";
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            if (n == nMaxIteration)
                                screenShot << 0 << " " << 0 << " " << 0 << std::endl;
                            else
                            {
                                n = n % 255;
                                screenShot << n << " " << n << " " << n << std::endl;
                            }
                        }
                    }
                    screenShot.close();
                    toggleScreenShot = false;
                }
                else
                {
                    #pragma omp parallel for
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            if (n == nMaxIteration)
                                Draw(x, y, olc::Pixel(0,0,0));
                            else
                                Draw(x, y, olc::Pixel(n,n,n));
                        }
                    }
                } 
                break;
            case 2:
                if(toggleScreenShot)
                {
                    std::ofstream screenShot("screen_shot_" + std::to_string(nScreenShotCount++) +".ppm");
                    screenShot << "P3\n" << std::to_string(ScreenHeight()) + " " + std::to_string(ScreenWidth()) + "\n" << "256\n";
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            screenShot << (n*n)%255 << " " << (n)%255 << " " << (n*3)%255 << std::endl;
                        }
                    }
                    screenShot.close();
                    toggleScreenShot = false;
                }
                else
                {
                    #pragma omp parallel for
                    for (int x = 0; x < ScreenWidth(); x++)
                    {
                        for (int y = 0; y < ScreenHeight(); y++)
                        {
                            int n = pFractalIterations[x * ScreenHeight() + y];
                            Draw(x, y, olc::Pixel(n*n, n, n*3, 255));
                        }
                    } 
                }
                break;
        }

        
        if(nFracMode == 0)
        {
            calcName = "Sequencial";
            color = olc::CYAN;
            fracModeColor[0] = olc::CYAN;
            fracModeColor[1] = olc::WHITE;
            fracModeColor[2] = olc::WHITE;
        }
        else if (nFracMode == 1)
        {
            calcName = "Paralelo";
            color = olc::YELLOW;
            fracModeColor[0] = olc::WHITE;
            fracModeColor[1] = olc::YELLOW;
            fracModeColor[2] = olc::WHITE;
        }
        else if (nFracMode == 2)
        {
            calcName = "Paralelo c/ Ins. Vetoriais.";
            color = olc::GREEN;
            fracModeColor[0] = olc::WHITE;
            fracModeColor[1] = olc::WHITE;
            fracModeColor[2] = olc::GREEN;
        }

        if(nColorMode == 0)
        {
            fracCol = fracColorCol[0] = olc::Pixel(126, 156, 247);
            fracColorCol[1] = olc::WHITE;
            fracColorCol[2] = olc::WHITE;
        }
        else if(nColorMode == 1)
        {
            fracColorCol[0] = olc::WHITE;
            fracCol = fracColorCol[1] = olc::DARK_GREY;
            fracColorCol[2] = olc::WHITE;
        }
        else if(nColorMode == 2)
        {
            fracColorCol[0] = olc::WHITE;
            fracColorCol[1] = olc::WHITE;
            fracCol = fracColorCol[2] = olc::RED;
        }

        constexpr int uiDist = 25;

        if(toggleBenchmark)
        {
            DrawString(0, 10, calcName, color, 2);
            DrawString(0, 35,  "Modelo de programacao: " + std::to_string(nFracMode + 1) + "/ 3", color, 2);
            DrawString(0, 60,  "Carga Computacional (iteracoes): " + std::to_string(nMaxIteration), colorIterations, 2);
            DrawString(0, 85,  "Coloracao: F" + std::to_string(nColorMode + 1) + "/ F3", fracCol, 2);
            DrawString(0, 110, "Tempo decorrido: " + std::to_string(fTime.count()) + "s", color, 3);
            DrawString(0, 145, "FPS: " + std::to_string(1/fTime.count()) + "", color, 3);
        }
        
        if(toggleHelp)
        {
            DrawString(700, 10 + uiDist *  0, "Controles:" , olc::Pixel(0, 255, 47), 2);
            DrawString(700, 10 + uiDist *  1, "Mover:" , olc::Pixel(2, 189, 36), 2);
            DrawString(700, 10 + uiDist *  2, "Segurar e Arrastar com mouse" , mouseColor, 2);
            DrawString(700, 10 + uiDist *  3, "Zoom:" , olc::Pixel(2, 189, 36), 2);
            DrawString(700, 10 + uiDist *  4, "Tecla: E - Zoom In" , zoomColor[0], 2);
            DrawString(700, 10 + uiDist *  5, "Tecla: Q - Zoom Out" , zoomColor[1], 2);
            DrawString(700, 10 + uiDist *  6, "Mudar carga computacional (iteracoes):" , olc::Pixel(2, 189, 36), 2);
            DrawString(700, 10 + uiDist *  7, "Tecla: UP " , olc::WHITE, 2);
            DrawString(700, 10 + uiDist *  8, "Tecla: DOWN " , olc::WHITE, 2);
            DrawString(700, 10 + uiDist *  9, "Mudar o modelo de programacao:" , olc::Pixel(2, 189, 36), 2);
            DrawString(700, 10 + uiDist * 10, "Tecla: 1 - Sequencial" , fracModeColor[0], 2);
            DrawString(700, 10 + uiDist * 11, "Tecla: 2 - Paralelo" , fracModeColor[1], 2);
            DrawString(700, 10 + uiDist * 12, "Tecla: 3 - Paralelo c/ Instruc. Vet." , fracModeColor[2], 2);
            DrawString(700, 10 + uiDist * 13, "Mudar o modelo de coloracao:" , olc::Pixel(2, 189, 36), 2);
            DrawString(700, 10 + uiDist * 14, "Tecla: F1 - Azul" , fracColorCol[0], 2);
            DrawString(700, 10 + uiDist * 15, "Tecla: F2 - Cinza" , fracColorCol[1], 2);
            DrawString(700, 10 + uiDist * 16, "Tecla: F3 - Vermelho." , fracColorCol[2], 2);
        }

        if(!toggleScreenShotView)
        {
            FillRect({0, ScreenHeight() - 85}, {370, 85}, olc::VERY_DARK_BLUE);
            DrawString(10, ScreenHeight() - 75, "I - Mostrar Benchmarks", olc::WHITE, 2);
            DrawString(10, ScreenHeight() - 50, "H - Controles e ajuda", olc::WHITE, 2);
            DrawString(10, ScreenHeight() - 25, "P - Print Screen", olc::WHITE, 2);
        }
		else
            toggleBenchmark = toggleHelp = false;

        DrawString(ScreenWidth() - 660, ScreenHeight() - 30, "Mandel2Us - github.com/lucaszm7/Mandel2Us", olc::Pixel(255,255,255,123), 2);
        
        return true;
	}

    ~MandelbrotFractal()
    {
        delete[] pFractalIterations;
        for(int i = 0; i < nNodesSize; ++i)
            delete[] pNodesParam[i];
        delete[] pNodesParam;
        setenv("SCREENSHOTCOUNT", std::to_string(nScreenShotCount).c_str(), 1);
    }

// Pan and Zoom Created with help of tutorials on channel @OneLoneCoder
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
    int nScreenWidth  = 1200;
    int nScreenHeight = 800;
    int nMyRank, nNodesSize;

    if(!UseMPI)
    {
        MPI::Init(argc, argv);
        MandelbrotFractal demo;
        if (demo.Construct(nScreenWidth, nScreenHeight, 1, 1, false, false))
            demo.Start();
        MPI::Finalize();
    }

    else
    {
        MPI::Init(argc, argv);

        nNodesSize = MPI::COMM_WORLD.Get_size();
        nMyRank = MPI::COMM_WORLD.Get_rank();
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
                CreateFractalParallelAVX({(int)pParam[0], (int)pParam[1]}, {(int)pParam[2], (int)pParam[3]}, 
                            {pParam[4], pParam[5]}, {pParam[6], pParam[7]}, 
                            pFractalIterations, (int)pParam[8], (int)pParam[3]);
                
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
    }
    return 0;
}
