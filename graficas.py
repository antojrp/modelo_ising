import numpy as np
import matplotlib.pyplot as plt
import os

# Tamaños del retículo (LxL)
L_values = [40, 80, 120, 160, 200]
# Parámetro crítico del modelo de Ising
nu = 0.8 
gamma = 7/4
beta = 1/8
# Obtener el directorio del script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Inicializar listas para almacenar datos de diferentes tamaños
data = {}

plt.rcParams.update({'font.size': 16})  # Tamaño general de las fuentes

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

# Definir la temperatura crítica del modelo de Ising
T_c = 2 / np.log(1 + np.sqrt(2))  # Aproximadamente 2.269

# Gráfica de magnetización vs temperatura con barras de error para diferentes tamaños
plt.figure(figsize=(8, 6))
T_teo = np.linspace(0.01, T_c, 500)
m_teo = (1 - 1 / np.sinh(2 / T_teo)**4)**(1 / 8)
m_teo[-1]=0
plt.plot(T_teo, m_teo, 'k-', linewidth=2, label="Solución \n de Onsager")   
for L, dataset in data.items():
    plt.errorbar(
        dataset["temperaturas"],
        dataset["medias_magnetizacion"],
        yerr=dataset["errores"],
        fmt="o",
        label=f"L = {L}",
        capsize=5,
    )
plt.axvline(T_c, color="gray", linestyle="--", label="T_c")  # Línea vertical en T_c
plt.xlabel("Temperatura (T)")
plt.ylabel("m")
plt.title("Magnetización frente a Temperatura", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_magnetizacion_vs_T.pdf")
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
plt.axvline(T_c, color="gray", linestyle="--", label="T_c")  # Línea vertical en T_c
plt.xlabel("Temperatura (T)")
plt.ylabel("$\chi_T$")
plt.title("Susceptibilidad frente a Temperatura", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_susceptibilidad_vs_T.pdf")
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
plt.axvline(T_c, color="gray", linestyle="--", label="T_c")  # Línea vertical en T_c
plt.xlabel("Temperatura (T)")
plt.ylabel("$U_4$")
plt.title("Cumulante frente a Temperatura", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_momento_vs_T.pdf")
plt.show()

# Gráfica de magnetización*L**b vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.errorbar(
        (1 - dataset["temperaturas"] / T_c) * L ** (1/nu),
        dataset["medias_magnetizacion"] * L**(beta/nu),
        yerr=dataset["errores"],
        fmt="o",
        label=f"L = {L}",
        capsize=5,
    )
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Magnetización promedio")
plt.title("Magnetización vs (1 - T/T_c) * L", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_magnetizacion_vs_Tc.pdf")
plt.show()

# Gráfica de susceptibilidad vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        (1 - dataset["temperaturas"] / T_c) * L ** (1/nu),
        dataset["susceptibilidades"] * L**(-(gamma)/nu),
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("(1 - T/T_c) * L")
plt.ylabel("Susceptibilidad magnética")
plt.title("Susceptibilidad vs (1 - T/T_c) * L", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_susceptibilidad_vs_Tc.pdf")
plt.show()

# Gráfica de momento vs (1 - T/T_c) * L para diferentes tamaños
plt.figure(figsize=(8, 6))
for L, dataset in data.items():
    plt.plot(
        (1 - dataset["temperaturas"] / T_c) * L**(1/nu),
        dataset["momentos"],
        "-o",
        label=f"L = {L}",
    )
plt.xlabel("(1 - T/T$_c$)L$^{1/\\nu}$")
plt.ylabel("Cumulante cuarto orden")
plt.title("Cumulante cuarto orden reescalado ($\\nu=0.8$)", fontsize=16)  # Reducir tamaño del título
plt.legend()
plt.grid()
plt.savefig("resultados/H=0_peque/ising_momento_vs_Tc.pdf")
plt.show()

