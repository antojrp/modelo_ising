import numpy as np
import matplotlib.pyplot as plt

# Cargar resultados finales por temperatura
# Formato por línea: T, <|m|>, Var(m), error, mc*tau, c
data = np.loadtxt("results.dat")

T = data[:, 0]
magnet = data[:, 1]
variance = data[:, 2]
error = data[:, 3]
tau = data[:, 4]
corr = data[:, 5]

# Cargar magnetización cruda (paso a paso)
magnetization_time_series = np.loadtxt("magnetization.dat")

# --- Gráfica de magnetización media vs temperatura ---
plt.figure(figsize=(8,6))
plt.errorbar(T, magnet, yerr=error, fmt='o-', capsize=3, label=r'$\langle |m| \rangle$')
plt.xlabel("Temperatura T")
plt.ylabel("Magnetización media")
plt.title("Magnetización absoluta promedio vs Temperatura")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("magnet_vs_temp.png")
plt.show()

# --- Gráfica de la serie temporal de magnetización ---
plt.figure(figsize=(10,4))
plt.plot(magnetization_time_series, lw=0.5)
plt.xlabel("Paso de medida")
plt.ylabel("Magnetización |m|")
plt.title("Evolución temporal de la magnetización")
plt.grid(True)
plt.tight_layout()
plt.savefig("magnet_time_series.png")
plt.show()

# --- Gráfica de la autocorrelación y tiempo de correlación ---
plt.figure(figsize=(8,6))
plt.plot(T, tau, 's-', label="mc·τ (tiempo de correlación)")
plt.xlabel("Temperatura T")
plt.ylabel("mc·τ")
plt.title("Tiempo de correlación vs Temperatura")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("tau_vs_temp.png")
plt.show()
