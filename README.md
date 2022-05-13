# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

## Build para Linux
> sudo apt-get install build-essential libglu1-mesa-dev libpng-dev

## Compilar
> mpic++ -o app Application.cpp -lX11 -lGL -lpthread -lpng -lstdc++fs -std=c++17

## Rodar
> mpirun --hostfile hosts ./app
