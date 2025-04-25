import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

# Parámetros del sistema

L = 100  # Tamaño de la red
T = 0.5  # Temperatura específica a graficar
directorio = "resultados/Kawasaki/"  # Directorio donde se encuentran los archivos
colors = cm.get_cmap('Dark2').colors
# Construir la ruta del archivo de clases
archivo_clases = os.path.join(directorio, f"kawasaki_clases_L_{L}_T_{T}_p=0.5.dat")

# Verificar si el archivo existe
if not os.path.exists(archivo_clases):
    print(f"Error: El archivo {archivo_clases} no existe.")
    exit(1)

# Leer los datos del archivo
try:
    datos = np.loadtxt(archivo_clases, delimiter=",", skiprows=1)  # Saltar la cabecera
    pasos = datos[:, 0]  # Primera columna: pasos
    clases = datos[:, 1:]  # Resto de las columnas: número de espines en cada clase
except Exception as e:
    print(f"Error al leer el archivo {archivo_clases}: {e}")
    exit(1)

# Configurar el tamaño de las fuentes
plt.rcParams.update({'font.size': 16})

# Graficar el número de espines en cada clase a lo largo del tiempo
plt.figure(figsize=(10, 6))
for i in range(5):  # Iterar sobre las columnas (clases)
    plt.plot(pasos, clases[:, i] + clases[:, 9 - i], color=colors[i], label=f"{4 - i} vecinos")

# Configurar la gráfica
plt.xlabel("Paso")
plt.ylabel("Número de espines")
plt.title(f"Evolución del número de espines\ncon vecinos con el mismo signo (L={L}, T={T}, p=0.5)")
plt.legend()
plt.grid()
plt.tight_layout()

# Guardar la gráfica como PDF en la misma ubicación que el archivo de datos
nombre_pdf = os.path.join(directorio, f"kawasaki_clases_L_{L}_T_{T}_p=0.5.pdf")
plt.savefig(nombre_pdf)
print(f"Gráfica guardada como {nombre_pdf}")

# Mostrar la gráfica
plt.show()