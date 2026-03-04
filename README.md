# KMC3D

Kinetic Monte Carlo - 3D grid with directions over a sphere

Julia script used to predict the charge mobility of organic semiconductors by means of Kinetic Monte Carlo simulations

---

The code allows for the definition of a 3D simulation grid on the fly by just defining (in order): 

- The number of **neighboring sites** for each site of the grid (*n*)
- The number of individual **simulation walks** to be performed (the final result would be the average of all of them)
- The direction of the **applied electric field** in form of a vector: X Y Z components
- The **number of hops** to be performed (*nhops*) along each simulation walk
- The **temperature (T)** at which the simulation is to be carried out
- The hole (or electron for electron transport) **reorganization energy** (lambda)
- The **distance** (d) between the grid sites (resulting from the modeling of a relistic thin film morphology)
- The statistical **electronic coupling** (tab) and **site energy differences** (eab) distributions

An example input would look like this: 

```
6 
5000
10000 0.0 0.0
2000
298.15
0.191
12.37
l -0.0105 0.5243 10.0
n 4.2958 195.5866 600.0
```

where the electric field must be provided in V/cm, the temperature in K, the hole reorganization energy in eV, and the distance in angstroms. 

For the electronic coupling distribution, the code accepts either "l" or "d" standing for lorentzian or double exponential distributions. The first two items correspond to the center of the distributions and the scale factor, respectively. The last item corresponds to the maximum value to be clamped (positive and negative).

For the energy difference distribution, only normal distribution "n" are accepted. The first two items correspond to the average and the standard deviation. The last item corresponds to the maximum value to be clamped (positive and negative). 

---

For executing the code, the following command might be used: 

```
julia --threads 8 KMC.jl IN > OUT
```

where IN is the input file and OUT is the output file. If no output file is supplied, the code will print the average hole mobility, together with the average number of hops, to the system standard output. 

The output mobility is in cm^2 / (V s)
