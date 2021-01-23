# Lotka-Volterra-modelisation
Small project for CentraleSupelec to model a continuous state system, the Lotka-Voltera equations in 2 dimensions with and without competition.

## Usage
Run the file plot.py
You can change the equation (competition, predaction or the corrected Lotka-Voltera equations) by (de)commenting the lines in the function lokta_volterra.

Both overflow.py and equation.py are my own helper modules used respectively for managing overflow problem with divergence of equations to prevent overflow error and to be able to use some numerical method (such as Euler explicit, implicit and Runge Kutta 4)


Have fun!
