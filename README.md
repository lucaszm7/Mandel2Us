# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

> ![Mandelbrot_Fractal__bigbig_adobe](https://user-images.githubusercontent.com/42661760/171764189-d58f25b9-5090-47b2-baf3-dd0992efab3b.gif)
> [![Video!](https://img.youtube.com/vi/9-DVTdkkEjQ/0.jpg)](https://www.youtube.com/watch?v=YOUTUBE_VIDEO_ID_HERE)

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
> <img width="399" alt="image" src="https://user-images.githubusercontent.com/42661760/174455667-a2f3f7c1-568b-43c5-bef0-83225e35b7d7.png">
> <img width="401" alt="image" src="https://user-images.githubusercontent.com/42661760/174455671-e1e14f78-8cbb-4af5-b087-adaeeb5d8b4b.png">
> <img width="400" alt="image" src="https://user-images.githubusercontent.com/42661760/174455681-3162f7fc-3bd3-48ed-b245-c1dc7d482dd9.png">
> <img width="399" alt="image" src="https://user-images.githubusercontent.com/42661760/174455767-94b58546-99e3-4f7e-83fb-34a4e0be868d.png">


## Build para Linux
> sudo apt-get install build-essential libglu1-mesa-dev libpng-dev

## Compilar com Makefile
Compilar && rodar
> make

## Rodar
> make run

## Compilar Manualmente
> mpic++ -fopenmp  Application.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17 -o app

## Rodar Manualmente
> mpirun --hostfile hosts ./app

## TODO:
- Intrinsec Functions
- CUDA(?)

## DONE:
- MakeFile
- Separate regions of Mandelbrot Between computers
- Already paralelize with OpenMP, but has to really get better.
