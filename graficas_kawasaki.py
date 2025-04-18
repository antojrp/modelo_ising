import numpy as np
import matplotlib.pyplot as plt
import os

# Parámetros del sistema
L = 50  # Tamaño de la red
T_i = 0.5  # Temperatura inicial
T_f = 3.5  # Temperatura final
n = 7  # Número de temperaturas

# Directorio donde se encuentran los archivos de energía
directorio = "resultados/Kawasaki/"

# Crear un arreglo para las temperaturas
temperaturas = np.linspace(T_i, T_f, n)

# Graficar la energía en función del tiempo para cada temperatura
plt.figure(figsize=(10, 6))

for T in temperaturas:
    # Usar str(T) para evitar ceros innecesarios a la derecha
    archivo_energia = os.path.join(directorio, f"kawasaki_energia_L_{L}_T_{str(T)}.dat")
    if os.path.exists(archivo_energia):
        # Leer las dos columnas: tiempo y energía
        datos = np.loadtxt(archivo_energia, delimiter=",")
        tiempo = datos[:, 0]  # Primera columna: tiempo
        energia = datos[:, 1]  # Segunda columna: energía
        plt.plot(tiempo, energia, label=f"T = {T}")  # Graficar energía vs tiempo
    else:
        print(f"Archivo no encontrado: {archivo_energia}")

# Configurar la gráfica
plt.xlabel("Tiempo (pasos)")
plt.ylabel("Energía")
plt.title("Energía vs Tiempo para diferentes temperaturas (Kawasaki)")
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()