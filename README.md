# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

![mandelbrot_compimido_265_medio_AdobeCreativeCloudExpress](https://user-images.githubusercontent.com/42661760/171763161-525cf337-887e-4b70-b207-8ca7c120337d.gif)

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
