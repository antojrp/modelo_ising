import numpy as np
import matplotlib.pyplot as plt
import os

# Tamaños del retículo (LxL)
L_values = [20,40,80,120,160]

# Obtener el directorio del script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Inicializar listas para almacenar datos de diferentes tamaños
data = {}

# Leer los datos para cada tamaño de red
for L in L_values:
    resultados_file = f"resultados/H=0_peque/ising_resultados_L_{L}.dat"

    # Verificar si el archivo existe
    if not os.path.exists(resultados_file):
        print(f"Advertencia: El archivo {resultados_file} no existe. Se omitirá este tamaño.")
        continue

    # Cargar los datos del archivo
    try:
        # Leer las columnas: T, Magnetización, Error, Susceptibilidad, Cumulante
        temperaturas, medias_magnetizacion, errores, susceptibilidades, momentos = np.loadtxt(
            resultados_file, delimiter=",", skiprows=1, unpack=True
        )
        data[L] = {
            "temperaturas": temperaturas,
            "medias_magnetizacion": medias_magnetizacion,
            "errores": errores,
            "susceptibilidades": susceptibilidades,
            "momentos": momentos,
        }
    except Exception as e:
        print(f"Error al cargar los datos del archivo {resultados_file}: {e}")
        continue

# Verificar si se cargaron datos
if not data:
    print("Error: No se cargaron datos de ningún archivo.")
    exit(1)

# Gráfica de magnetización vs temperatura con barras de error para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.errorbar(
        dataset["temperaturas"],
        dataset["medias_magnetizacion"],
        yerr=dataset["errores"],
        fmt="o",
        label=f"L = {L}",
        capsize=5,
    )
plt.xlabel("Temperatura (T)")
plt.ylabel("Magnetización promedio")
plt.title("Magnetización vs Temperatura")
plt.legend()
plt.grid()
plt.show()

# Gráfica de susceptibilidad vs temperatura para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        dataset["temperaturas"],
        dataset["susceptibilidades"],
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("Temperatura (T)")
plt.ylabel("Susceptibilidad magnética")
plt.title("Susceptibilidad vs Temperatura")
plt.legend()
plt.grid()
plt.show()

# Gráfica de momento vs temperatura para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        dataset["temperaturas"],
        dataset["momentos"],
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("Temperatura (T)")
plt.ylabel("Cumulante cuarto orden")
plt.title("Cumulante cuarto orden vs Temperatura")
plt.legend()
plt.grid()
plt.show()

# Definir la temperatura crítica del modelo de Ising
T_c = 2 / np.log(1 + np.sqrt(2))  # Aproximadamente 2.269

# Gráfica de magnetización*L**b vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.errorbar(
        (1 - dataset["temperaturas"] / T_c) * L,
        dataset["medias_magnetizacion"] * L**(1/8),
        yerr=dataset["errores"],
        fmt="o",
        label=f"L = {L}",
        capsize=5,
    )
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Magnetización promedio")
plt.title("Magnetización vs (1 - T/T_c) * L")
plt.legend()
plt.grid()
plt.show()

# Gráfica de susceptibilidad vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        (1 - dataset["temperaturas"] / T_c) * L,
        dataset["susceptibilidades"],
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Susceptibilidad magnética")
plt.title("Susceptibilidad vs (1 - T/T_c) * L")
plt.legend()
plt.grid()
plt.show()

# Gráfica de momento vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        (1 - dataset["temperaturas"] / T_c) * L,
        dataset["momentos"],
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Cumulante cuarto orden")
plt.title("Cumulante cuarto orden vs (1 - T/T_c) * L")
plt.legend()
plt.grid()
plt.show()
