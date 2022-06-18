# Fractal de Mandelbrot construido em C++ com MPI e OpenMP, e a parte gráfica com olcPixelGameEngine

![Mandelbrot_Fractal__bigbig_adobe](https://user-images.githubusercontent.com/42661760/171764189-d58f25b9-5090-47b2-baf3-dd0992efab3b.gif)

## Features!

### Different coloring algorithms:
<img width="399" alt="image" src="https://user-images.githubusercontent.com/42661760/174455667-a2f3f7c1-568b-43c5-bef0-83225e35b7d7.png">
<img width="401" alt="image" src="https://user-images.githubusercontent.com/42661760/174455671-e1e14f78-8cbb-4af5-b087-adaeeb5d8b4b.png">
<img width="400" alt="image" src="https://user-images.githubusercontent.com/42661760/174455681-3162f7fc-3bd3-48ed-b245-c1dc7d482dd9.png">
<img width="399" alt="image" src="https://user-images.githubusercontent.com/42661760/174455767-94b58546-99e3-4f7e-83fb-34a4e0be868d.png">


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
