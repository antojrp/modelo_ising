import numpy as np
import matplotlib.pyplot as plt
import os

L = 40  # Tamaño del retículo (LxL)

# Obtener el directorio del script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Construir rutas absolutas para los archivos
magnetizacion_file = "resultados/H=0/ising_magnetizacion_L_"+str(L)+".dat"
susceptibilidad_file = "resultados/H=0/ising_susceptibilidad_L_"+str(L)+".dat"

temperaturas, medias_magnetizacion, varianzas_magnetizacion = np.loadtxt(magnetizacion_file, delimiter=",", unpack=True)
temperaturas_susc, susceptibilidades = np.loadtxt(susceptibilidad_file, delimiter=",", unpack=True)

# Gráfica de magnetización vs temperatura con barras de error
plt.figure(figsize=(8, 6))
plt.errorbar(
    temperaturas,
    medias_magnetizacion,
    yerr=np.sqrt(varianzas_magnetizacion),
    fmt="o",
    label="Magnetización",
    capsize=5,
)
plt.xlabel("Temperatura (T)")
plt.ylabel("Magnetización promedio")
plt.title("Magnetización vs Temperatura")
plt.legend()
plt.grid()
plt.show()

# Gráfica de susceptibilidad vs temperatura
plt.figure(figsize=(8, 6))
plt.plot(temperaturas_susc, susceptibilidades, "-o", label="Susceptibilidad magnética")
plt.xlabel("Temperatura (T)")
plt.ylabel("Susceptibilidad magnética")
plt.title("Susceptibilidad vs Temperatura")
plt.legend()
plt.grid()
plt.show()

# Definir la temperatura crítica del modelo de Ising
T_c = 2 / np.log(1 + np.sqrt(2))  # Aproximadamente 2.269

#Gráfica de magnetización*L**b vs (1 - T/T_c) * L
plt.figure(figsize=(8, 6))
plt.errorbar(
    (1 - temperaturas / T_c) * L,
    medias_magnetizacion * L**(1/8),
    yerr=np.sqrt(varianzas_magnetizacion),
    fmt="o",
    label="Magnetización",
    capsize=5,
)
plt.xlabel("(1 - T/T_c) * L")   
plt.ylabel("Magnetización promedio")
plt.title("Magnetización vs (1 - T/T_c) * L")
plt.legend()
plt.grid()
plt.show()

# Gráfica de susceptibilidad vs (1 - T/T_c) * L
plt.figure(figsize=(8, 6))
plt.plot((1 - temperaturas / T_c) * L, susceptibilidades, "-o", label="Susceptibilidad magnética")
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Susceptibilidad magnética")
plt.title("Susceptibilidad vs (1 - T/T_c) * L")
plt.legend()
plt.grid()
plt.show()

