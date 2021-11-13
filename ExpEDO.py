#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plotter
from scipy.integrate import odeint
import numpy as np
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

window = tk.Tk()
window.title("Dextroit")
window.configure(background="#15224F",width=500,height=400)
##construção de entries e labels de concentração

variaveis = [0,0,0,0,0,0,0,0,0]
identidade_entries = []
contador = 0


#Interface de entrada de valores
def salvar_valores():
    global variaveis
    concentracoes_inicias = []
    contador = 0
    while contador < len(variaveis):
        try:
            concentracoes_inicias.append(float(variaveis[contador].get()))
        except:
            concentracoes_inicias.append(0.0)
        contador += 1
    print(concentracoes_inicias)
    S0 = concentracoes_inicias[:7]
    global Eex
    global Een
    Eex = concentracoes_inicias[7]
    Een = concentracoes_inicias[8]
    parametros(S0)
    return

#Parâmetros para equações

# from dfEDO.py import dfB
def parametros(S0):
    global K
    global Eex
    global Een
    K1 = [1, 1, 4.04, 4.7, 183.33, 0.54, 519.04, 0.21,313.88, 0.36, 480.76, 0.26, 563.63, 0.22]
    K2 = [1, 0.17, 0.89, 0.17, 0.89, 0.17, 0.89, 0.17, 0.89]
    K = K1 + K2

    # print(len(S0))
    # Exoenzima
    #Eex = 1
    # Endoenzima
    #Een = 1

    # Eixo "x", tempo de reação em segundos (s):
    t_grafico = np.linspace(0, 400)

    # O Sistema de Equações Diferenciais (EDO):
    # odefun = @(t,G) dfEDO(G,K,Eex,Een);

    # Solução do sistema de EDO:
    # [X, Y] = ode45(odefun, t_grafico, S0);
    y = odeint(dfB, S0, t_grafico)
    # print(y)

    figure, eixo = plotter.subplots(1, 2)
    plotter.subplots_adjust(hspace=2)
    eixo[0].set_title("G1 visível")
    eixo[0].legend(loc="upper left")

    eixo[0].plot(t_grafico, y, label = ["G1","G2","G3","G4","G5","G6","G7"])
    eixo[0].set_xlabel('Tempo')

    eixo[0].set_ylabel('Concentração')
    # plt.figure("G1 visível")
    # plt.plot(t_grafico, y)
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
    #plt.figure("G1 Invisível")
    #plt.plot(t_grafico, y2)
    #plt.show()
    eixo[1].set_title("G1 invisível")

    eixo[1].plot(t_grafico, y2, label = ["G2","G3","G4","G5","G6","G7"])

    eixo[1].set_xlabel('Tempo')

    eixo[1].set_ylabel('Concentração')

    eixo[1].legend(loc="upper right")

    plotter.show()
    return


##Equações diferenciais
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


#Construção de interface
while contador < len(variaveis):
    variaveis[contador] = tk.StringVar(window)
    concentracao = tk.Entry(textvariable = variaveis[contador])  # botao de iniciar ativ
    concentracao.config(width=15)
    identidade_entries.append(concentracao)
    if contador < 7:
        textos = tk.Label(text="Concentração G" + str(contador+1), fg="white", bg="#1B1F49", width=35, height=3)
    elif contador == 7:
        textos = tk.Label(text="Concentração Exoenzima", fg="white", bg="#1B1F49", width=35, height=3)
    elif contador == 8:
        textos = tk.Label(text="Concentração Endoenzima", fg="white", bg="#1B1F49", width=35, height=3)
    textos.grid(row=contador+1,column=1)
    concentracao.grid(row=contador + 1, column=2)
    contador += 1

botao_enviar = tk.Button(window, width=15, text="Enviar!", bg="#aff587", command=salvar_valores)
botao_enviar.grid(row=10, column=1)

window.mainloop()
