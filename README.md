# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

## Build para Linux
> sudo apt-get install build-essential libglu1-mesa-dev libpng-dev

## Compilar com Makefile
Compilar && rodar
> make

## Rodar
> make run

## Compilar
> mpic++ -fopenmp  Application.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17 -o app

## Rodar
> mpirun --hostfile hosts ./app

## TODO:
- Separate regions of Mandelbrot Between computers
- Already paralelize with OpenMP, but has to really get better.

## DONE:
- MakeFile
