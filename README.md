# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

## Build: Linux
 1. Clone o repositório
   ```
   git clone https://github.com/lucaszm7/Mandel2Us.git
   ```
 2. Instale as dependências
   ```
   make install
   ```
 3. Compile o Programa 
   ```
   make gen
   ```
 4. Rode :)
   ```
   make run
   ```
   
 ## Build Windows
  - Com Msys2 + mingw-64
   ```
   g++ -fopenmp Application.cpp -luser32 -lgdi32 -lopengl32 -lgdiplus -lShlwapi -ldwmapi -lstdc++fs -static -std=c++17 -O3 -mavx2  -o app
   ```
  - Com MSVC (Alpha)
   ```
   cl /EHsc /openmp /O2 /Ot /std:c++17 /arch:AVX2 Application.cpp
   ```
 
---

![Mandelbrot_Fractal__bigbig_adobe](https://user-images.githubusercontent.com/42661760/171764189-d58f25b9-5090-47b2-baf3-dd0992efab3b.gif)
- Video Showcase: https://www.youtube.com/watch?v=9-DVTdkkEjQ

## Features!

### Clusterization with MPI!
> Dynamically divide the fractal in N nodes, beeing N <= Screen Width
> Return the value from each nodes to the master node

### Parallelization with OpenMP!
> Dynamic divide the portion of the fractal in each node, using MPI

### Multi-PLataform Rendering!
> Using olcPixelGameEngine to Draw the fractal in a window!

### Handle pan and Zoom!
> Implemented pan and zoom in the fractal, supporting before features!

### Different coloring algorithms:
![screen_shot_0](https://user-images.githubusercontent.com/42661760/177466041-be9d8a8d-4bec-4a08-9712-f9325af655a4.png)
> <img width="599" alt="_mandel2us_print_mark_2" src="https://user-images.githubusercontent.com/42661760/177465943-2001bb61-cc1e-4e23-86af-1db0c15bf034.png">
![screen_shot_1_](https://user-images.githubusercontent.com/42661760/177465956-ebf08e27-6543-49bf-9ab3-b616c5a6d45c.png)
![screen_shot_2](https://user-images.githubusercontent.com/42661760/177465979-8228379d-23e2-43e4-a7d4-bf00e6678661.png)
![screen_shot_06](https://user-images.githubusercontent.com/42661760/177466011-29821b44-ccf8-41a8-838c-b391c32974a6.png)

## Compilar Manualmente
> mpic++ -fopenmp  Application.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17 -o app

## Rodar Manualmente
> mpirun --hostfile hosts ./app

## TODO:
- CUDA(?)

## DONE:
- UI
- Intrinsec Functions
- MakeFile
- Separate regions of Mandelbrot Between computers
- Already paralelize with OpenMP, but has to really get better.
