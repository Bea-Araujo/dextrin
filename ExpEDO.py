#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np


def dfB(y, x):
    global K
    global Eex
    global Een
    G = y
    dG = []

    Ki_1ex = K[0]
    Ki_1en = K[13]

    Ki_2en = K[14]
    ks_2ex = K[1]
    Km_2ex = K[2]

    ks_3ex = K[3]
    Km_3ex = K[4]

    ks_4ex = K[5]
    Km_4ex = K[6]
    ks_4en = K[15]
    Km_4en = K[16]

    ks_5ex = K[7]
    Km_5ex = K[8]
    ks_5en = K[17]
    Km_5en = K[18]

    ks_6ex = K[9]
    Km_6ex = K[10]
    ks_6en = K[19]
    Km_6en = K[20]

    ks_7ex = K[11]
    Km_7ex = K[12]
    ks_7en = K[21]
    Km_7en = K[22]

    # numpy.append(G, S0[0])
    # numpy.append(G, S0[1])
    # numpy.append(G, S0[2])
    # numpy.append(G, S0[3])
    # numpy.append(G, S0[4])
    # numpy.append(G, S0[5])
    # numpy.append(G, S0[6])
    dG1 = (ks_7ex * Eex * G[6] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (ks_6ex * Eex * G[5] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (ks_5ex * Eex * G[4] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (
                  ks_4ex * Eex * G[3] / (
                  1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[
              5] / Km_6ex + G[6] / Km_7ex)) + (ks_3ex * Eex * G[2] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (ks_2ex * Eex * G[1] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex))
    dG2 = (ks_3ex * Eex * G[2] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (ks_5en * Een * G[4] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)) + \
          2 * (ks_4en * Een * G[3] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)) - (
                  ks_2ex * Eex * G[1] / (
                  1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[
              5] / Km_6ex + G[6] / Km_7ex))
    dG3 = (ks_4ex * Eex * G[3] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) + (ks_7en * Een * G[6] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[
        6] / Km_7en)) + 2 * (ks_6en * Een * G[5] / (1 + G[0] / Ki_1en + G[1] /
                                                    Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[
                                                        6] / Km_7en)) + (ks_5en * Een * G[4] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)) - (
                  ks_3ex * Eex * G[2] / (
                  1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[
              5] / Km_6ex + G[6] / Km_7ex))
    dG4 = (ks_7en * Een * G[6] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)) + (
                  ks_5ex * Eex * G[4] / (
                  1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[
              5] / Km_6ex + G[6] / Km_7ex)) - (
                  (ks_4ex * Eex * G[3] / (
                          1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[
                      5] / Km_6ex + G[6] / Km_7ex)) + (ks_4en * Een * G[3] / (
                  1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[
              6] / Km_7en)))
    dG5 = (ks_6ex * Eex * G[5] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) - ((ks_5ex * Eex * G[4] / (1 + G[0] / Ki_1ex + G[1] / Km_2ex +
                                                 G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
                                                     6] / Km_7ex)) + (ks_5en * Een * G[4] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)))
    dG6 = (ks_7ex * Eex * G[6] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] / Km_7ex)) - ((ks_6ex * Eex * G[5] / (1 + G[0] / Ki_1ex + G[1] / Km_2ex +
                                                 G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
                                                     6] / Km_7ex)) + (ks_6en * Een * G[5] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en)))

    dG7 = -(ks_7ex * Eex * G[6] / (
            1 + G[0] / Ki_1ex + G[1] / Km_2ex + G[2] / Km_3ex + G[3] / Km_4ex + G[4] / Km_5ex + G[5] / Km_6ex + G[
        6] /
            Km_7ex)) - (ks_7en * Een * G[6] / (
            1 + G[0] / Ki_1en + G[1] / Ki_2en + G[3] / Km_4en + G[4] / Km_5en + G[5] / Km_6en + G[6] / Km_7en))

    dG = [dG1, dG2, dG3, dG4, dG5, dG6, dG7]
    # dG.append(0)
    # dG.append(0)
    # dG.append(0)
    # dG.append(0)
    # dG.append(0)
    # dG.append(0)
    # dG.append(dG7)

    return dG


# from dfEDO.py import dfB

K1 = [1, 1, 4.04, 4.7, 183.33, 0.54, 519.04, 0.21,
      313.88, 0.36, 480.76, 0.26, 563.63, 0.22]
K2 = [1, 0.17, 0.89, 0.17, 0.89, 0.17, 0.89, 0.17, 0.89]
K = K1 + K2

S0 = [0, 0, 0, 0, 0, 0, 10]
# print(len(S0))
# Exoenzima
Eex = 1
# Endoenzima
Een = 1

# Eixo "x", tempo de reação em segundos (s):
t_grafico = np.linspace(0, 400)

# O Sistema de Equações Diferenciais (EDO):
# odefun = @(t,G) dfEDO(G,K,Eex,Een);

# Solução do sistema de EDO:
# [X, Y] = ode45(odefun, t_grafico, S0);
y = odeint(dfB, S0, t_grafico)
# print(y)
plt.figure(1)
plt.plot(t_grafico, y)

# plt.plot(t_grafico, y[1:])

print("Cabo")
# print(y)
# print(y[1:])
# print(t_grafico)
y2 = []
i = 0
while i < len(y):
    temp = y[i][1:]
    print(temp)
    y2.append(temp)
    i += 1

print(y2)
plt.figure(2)
plt.plot(t_grafico, y2)
plt.show()
