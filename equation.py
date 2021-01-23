# Goal : implement the RungeKutta algorithm (Euler algorithm does not converge)
# Library

from math import *
from overflow import *


# Coupled equations

class Coupled:

    @staticmethod
    def Euler_modified(f, g, a, b, x0, y0, n, *X1, **Y1):
        ecart = 0

        h = (b-a)/n
        t, x, y = a, x0, y0
        T, X, Y = [a], [x0], [y0]

        for i in range(n):
            x, y = x+h*f(t+h/2, x+h/2*f(t, x, y), y+h/2*g(t, x, y)
                         ), y+h*g(t+h/2, x+h/2*f(t, x, y), y+h/2*g(t, x, y))
            X.append(x)
            Y.append(y)
            t += h
            T.append(t)

            if(X1 and Y1):
                ecart += (X1[i] - x)**2 + (Y1[i] - y)**2
        if X1 and Y1:
            return T, X, Y, ecart
        return T, X, Y

    @staticmethod
    def Heun(f, g, a, b, x0, y0, n, *X1, **Y1):
        ecart = 0

        h = (b-a)/n
        t, x, y = a, x0, y0
        T, X, Y = [a], [x0], [y0]

        for i in range(n):
            x, y = x+h/2*(f(t, x, y)+f(t+h, x+h*f(t, x, y), y+h*g(t, x, y))), y + \
                h/2*(g(t, x, y)+g(t+h, x+h*f(t, x, y), y+h*g(t, x, y)))
            X.append(x)
            Y.append(y)
            t += h
            T.append(t)

            if(X1 and Y1):
                ecart += (X1[i] - x)**2 + (Y1[i] - y)**2

        if X1 and Y1:
            return T, X, Y, ecart
        return T, X, Y

    @staticmethod
    def RungeKutta(f, g, a, b, x0, y0, n, X1=[], Y1=[]):
        ecart = 0

        h = (b-a)/n
        t, x, y = a, x0, y0
        T, X, Y = [a], [x0], [y0]

        for i in range(n):
            k1x = f(t, x, y)
            k1y = g(t, x, y)

            k2x = f(t+h/2, x+h/2*k1x, y+h/2*k1y)
            k2y = g(t+h/2, x+h/2*k1x, y+h/2*k1y)

            k3x = f(t+h/2, x+h/2*k2x, y+h/2*k2y)
            k3y = g(t+h/2, x+h/2*k2x, y+h/2*k2y)

            k4x = f(t+h, x+h*k3x, y+h*k3y)
            k4y = g(t+h, x+h*k3x, y+h*k3y)

            x, y = x+h/6*(k1x+2*k2x+2*k3x+k4x), y+h/6*(k1y+2*k2y+2*k3y+k4y)

            X.append(x)
            Y.append(y)
            t += h
            T.append(t)

            if len(X1) != 0 and len(Y1) != 0:
                ecart += (X1[i] - x)**2 + (Y1[i] - y)**2
        if len(X1) != 0 and len(Y1) != 0:
            return T, X, Y, ecart
        return T, X, Y

# Multi-coupled equations


class MultiCoupled:

    @staticmethod
    def Euler_modified(F, X0, a, b, n, X1=[]):
        ecart = 0

        h = (b-a)/n
        t, X = a, X0
        T, Xn = [a], [X0]

        for i in range(n):
            X = X+h*F(t+h/2, X+h/2*F(t, X))
            t += h

            Xn.append(X)
            T.append(t)

            if len(X1) != 0:
                for j in range(len(X0)):
                    ecart += (X[j][i] - X[j])**2
        if len(X1) != 0:
            return T, Xn, ecart
        return T, Xn

    @staticmethod
    def RungeKutta(F, X0, a, b, n, X1=[], tp="only"):
        ecart = 0
        ecarts = [0 for i in range(len(X0))]

        h = (b-a)/n
        t, X = a, X0
        T, Xn = [a], [X0]

        for i in range(n):
            if len(X1) != 0:
                if len(X) == 2:
                    ecart += (X1[i][0] - X[0])**2 + (X1[i][1] - X[1])**2
                elif tp == "only":
                    for j in range(len(X0)):
                        ecart += (X1[i][j] - X[j])**2
                elif tp == "multi":
                    for j in range(len(X0)):
                        ecarts[j] += (X1[i][j] - X[j])**2

            k1 = F(t, X)
            k2 = F(t+h/2, X+h/2*k1)
            k3 = F(t+h/2, X+h/2*k2)
            k4 = F(t+h, X+h*k3)

            X = overflow(X+h/6*(k1+2*k2+2*k3+k4))

            Xn.append(X)
            t += h
            T.append(t)

        for j in range(len(X0)):
            ecarts[j] = sqrt(ecarts[j]/len(X0))/sqrt(len(X0))

        if len(X1) != 0:
            if tp == "cost":
                return T, Xn, sqrt(ecart/len(X0))/sqrt(len(X0))
            elif tp == "only":
                return T, Xn, ecart/2
            elif tp == "multi":
                return T, Xn, ecarts
        return T, Xn
