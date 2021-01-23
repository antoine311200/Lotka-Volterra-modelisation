import matplotlib.pyplot as plt
import numpy as np

from equation import *


def corrected_equation(r):
    return lambda t, X: np.array([X[0]*(r[0][0]*(1-X[0]/100)-r[0][1]*X[1]), X[1]*(r[1][1]*X[0]-r[1][0])])


def equation_predation(r):
    return lambda t, X: np.array([X[0]*(r[0][0]-r[0][1]*X[1]), X[1]*(r[1][1]*X[0]-r[1][0])])


def equation_intern_competition(pm):
    return lambda t, X: np.array([X[0]*(pm['a']-pm['b']*X[1]-pm['e']*X[0]), X[1]*(-pm['c']+pm['d']*X[0]-pm['f']*X[1])])


space = np.linspace(1, 10, 4)
initial_populations = [[space[i], space[i]] for i in range(space.shape[0])]
initial_population = [4, 4]
time = np.linspace(0, 20, 2000)


def lotka_volterra(initial_population, time):
    r = [[4, 3], [4, 2]]
    pm = {'a': 4, 'b': 3, 'c': 4, 'd': 2, 'e': 0.4, 'f': 0.5}
    # Pour intern-competition, utiliser pm. Pour les autres, utiliser r.
    run = np.transpose(MultiCoupled.RungeKutta(
        equation_predation(r), initial_population, 0, 20, 1999)[1])
    """
    run = np.transpose(MultiCoupled.RungeKutta(
        equation_intern_competition(pm), initial_population, 0, 20, 1999)[1])
    """
    """
    run = np.transpose(MultiCoupled.RungeKutta(
        corrected_equation(r), initial_population, 0, 20, 1999)[1])
    """
    return run


def show(initial_population, time, gradient=1):
    value = lotka_volterra(initial_population, time)
    # Change the body of the Lotka-Volterra function to change the equation
    plt.subplot(221)
    plt.plot(time, value[0], color=(1, gradient, 0),
             label=str(initial_population[0]))
    plt.xlabel('Time')
    plt.ylabel('Preys')
    plt.legend(loc='upper right', title='Initial Preys')
    plt.subplot(223)
    plt.plot(time, value[1], color=(0, gradient, 1),
             label=str(initial_population[1]))
    plt.xlabel('Time')
    plt.ylabel('Predators')
    plt.legend(loc='upper right', title='Initial Predators')
    plt.subplot(222)
    plt.plot(value[0], value[1], color=(1, gradient, 0))
    plt.xlabel('Preys')
    plt.ylabel('Predators')
    plt.get_current_fig_manager().resize(1400, 900)
    return plt


def shows(initial_populations, time):
    for i in range(len(initial_populations)):
        gradient = i/len(initial_populations)
        print(initial_populations[i])
        show(initial_populations[i], time, gradient=gradient)
        print(initial_populations[i], "shown")
    plt.get_current_fig_manager().resize(1400, 900)
    return plt


shows(np.sort(initial_populations, axis=0), time)
plt.show()
