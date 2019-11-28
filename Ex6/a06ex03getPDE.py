#!/usr/bin/env python3

# ------------------------------------------------------------------------\
# Assignment 4, Exercise 1c                                               |
#                                                             submitted by|
#                                                                         |
#                        Kagan Atci | 338131 | Physical Engineering, M.Sc.|
#                     Navneet Singh | 380443 | Scientific Computing, M.Sc.|
#                   Riccardo Parise | 412524 | Scientific Computing, M.Sc.|
#        Daniel V. Herrmannsdoerfer | 412543 | Scientific Computing, M.Sc.|
#                                                                         |
#                                                          in  Python 3.5 |
# ------------------------------------------------------------------------/


def a06ex03getPDE(k, f, g, N):
    print("a06func")
    return


def main():

    def f(x): return x ** 2 + 3 * x - 4
    def k(x): return 2 * x
    def g(x): return 0
    n = 10

    a06ex03getPDE(k, f, g, n)


if __name__ == '__main__':
    main()
