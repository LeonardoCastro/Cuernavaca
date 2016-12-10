# Cuernavaca Highway Simulations

This repository contains the numerical simulations done to study the phase simulations of the Cuernavaca Bypass, in central Mexico.

The aim of this study is to know if the public work to expand the bypass is going to bring benefits to the population living in the surrounding areas.

S1 and S2 are the two different senses found in the highway. S1 is the North-South direction and S2 the South-North Direction.

The suffix CR correspond to the state of the highway before the expansion. The latter consider the creation of two different highways, the first having the same on and off-ramps and three lanes instead of two. Because of the ramps, this first system is called the Local one. The second highway consists of two lanes and the elimination of all on and off-ramps. Because of these facts, we call it the Express-Pass.

The sripts are written in [Julia](http://julialang.org/), v0.5. A module called (Highway)[https://github.com/LeonardoCastro/Highway] is used to create the different highway and to consider the dynamics of the vehicles in it.
