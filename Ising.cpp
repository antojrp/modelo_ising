# include <iostream>
# include <random> // Librería para números aleatorios
# include <time.h>
# include <fstream>
# include <math.h>

using namespace std;
const int L=40; // Tamaño de la red
const int n=79; // Número de temperaturas
const int M0=1000; // Número de pasos de termalización
const int M = 8192; // Número de pasos de medida por temperatura


// Inicializar el generador de números aleatorios
std::mt19937 rng(time(nullptr)); // Semilla basada en el tiempo

void generar_datos(int s[L][L], int N)
{
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Distribución uniforme continua entre 0.0 y 1.0
    float p = 0.5; // Probabilidad de asignar -1
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            if (dist(rng) < p) // Ajustar la probabilidad para asignar -1
            {
                s[i][j] = -1;
            }
            else // Caso contrario, asignar 1
            {
                s[i][j] = 1;
            }
        }
    }
}

double variacion_energia(int s[L][L], int n, int m, int N, double J, double H)
{
    return  2*J*s[n][m]*(s[n][(m+1)%N]+s[n][(m+N-1)%N]+s[(n+1)%N][m]+s[(n+N-1)%N][m]) + 2 * H * s[n][m];
}

void precalcular_factores(double h[5][2] , double J, double H, double T) {
    for (int j = -4, idx = 0; j <= 4; j += 2, idx++) {
        for (int i = -1, idy = 0; i <= 1; i += 2, idy++) {
            h[idx][idy] = exp(-2*J*j / T - 2*H*i/T); // Precálculo de los factores de Boltzmann-Gibbs
        }
    }
}

double magnetizacion(int s[L][L], int N)
{
    double M = 0;
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            M += s[i][j];
        }
    }
    return abs(M / (N*N)); // Magnetización promedio
}

void ising(int s[L][L], double T, int N, double J, double H, double h[5][2]) {
    std::uniform_int_distribution<int> dist_pos(0, N - 1); // Distribución para posiciones
    std::uniform_real_distribution<double> dist_prob(0.0, 1.0); // Distribución para probabilidades

    int n = dist_pos(rng);
    int m = dist_pos(rng);

    // Calcular la suma de los vecinos
    int suma_vecinos = s[n][(m + 1) % N] + s[n][(m + N - 1) % N] + s[(n + 1) % N][m] + s[(n + N - 1) % N][m];
    int idx = (s[n][m]*suma_vecinos + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
    int idy = (s[n][m] + 1) / 2; // Mapear el espín a índices (0, 1) 
    double p;


    p=min(1.0, h[idx][idy]); // Asegurarse de que p no exceda 1.0
    double e = dist_prob(rng);

    if (e < p) {
        s[n][m] = -s[n][m]; // Cambiar el estado del espín
    }
}

void crear_fichero(string nombre)
{
    ofstream salida(nombre);
    salida << "";
    salida.close();
}

void escribir_datos(int s[L][L], int N, string direccion)
{
    ofstream salida(direccion, ios::app);
    for (int i = 0; i < N; i++) {           // Recorre filas
        for (int j = 0; j < N; j++) {       // Recorre columnas
            salida << s[i][j];
            if (j < N - 1) salida << ",";   // Añade coma si no es la última columna
        }
        salida << "\n";
    }
    salida << "\n";
    salida.close();
}

string format_double(double value) {
    string str = to_string(value);
    str.erase(str.find_last_not_of('0') + 1, string::npos); // Elimina ceros a la derecha
    if (str.back() == '.') {
        str.pop_back(); // Elimina el punto decimal si es el último carácter
    }
    return str;
}

// Función para calcular la susceptibilidad magnética
double susceptibilidad_magnetica(double magnetizaciones[], int nsim, double T) {
    double M_promedio = 0;
    double M2_promedio = 0;

    for (int i = 0; i < nsim; i++) {
        M_promedio += magnetizaciones[i];
        M2_promedio += magnetizaciones[i] * magnetizaciones[i];
    }

    M_promedio /= nsim;
    M2_promedio /= nsim;

    return (M2_promedio - M_promedio * M_promedio) / T;
}

// Función para escribir resultados en un archivo
void escribir_resultados(const string& direccion, double temperaturas[], double medias[], double varianzas[], int n) {
    ofstream salida(direccion);
    for (int i = 0; i < n; i++) {
        salida << temperaturas[i] << "," << medias[i] << "," << varianzas[i] << "\n";
    }
    salida.close();
}

double energia_total(int s[L][L], int L, double J, double H) {
    double E = 0;
    for (int j = 0; j < L; j++) {
        for (int i = 0; i < L; i++) {
            E -= J * s[i][j] * (s[i][(j + 1) % L] + s[(i + 1) % L][j]); // Interacción con los vecinos
            E -= H * s[i][j];
        }
    }
    return E;
}

int main()
{
    int s[L][L];
    double T_i = 0.1; // Temperatura inicial
    double T_f = 4; // Temperatura final

    double J = 1; // Interacción entre espines
    double H = 0; // Campo magnético externo



    string direccion_magnetizacion = "resultados/H=" + format_double(H) + "/ising_magnetizacion_L_" + to_string(L) + ".dat";
    string direccion_susceptibilidad = "resultados/H=" + format_double(H) + "/ising_susceptibilidad_L_" + to_string(L) + ".dat";
    string direccion_energia = "resultados/H=" + format_double(H) + "/ising_energia_L_" + to_string(L) + ".dat";
    ofstream salida_energia(direccion_energia); // Crear archivo para la energía

    double temperaturas[n];
    double medias_magnetizacion[n];
    double varianzas_magnetizacion[n];
    double susceptibilidades[n];
    generar_datos(s, L);
    for (int j = 0; j < n; j++) {
        double T = T_f - j * (T_f - T_i) / (n-1); // Temperatura decreciente
        double h[5][2];
        precalcular_factores(h, J, H, T); // Precalcular factores de Boltzmann-Gibbs
        temperaturas[j] = T;
        string direccion = "resultados/H=" + format_double(H) + "/ising_data_L_" + to_string(L) + "_T_" + format_double(T) + ".dat";
        crear_fichero(direccion);


        cout << "Procesando temperatura " << j + 1 << " de " << n << " (T = " << T << ")" << endl;

        double magnetizaciones[M]; // Arreglo para almacenar las magnetizaciones finales
        escribir_datos(s, L, direccion);
        // Termalización
        for (int i = 0; i < M0 * L * L; i++) {
            ising(s, T, L, J, H, h);
            if (i % (L * L) == 0) {
                escribir_datos(s, L, direccion); // Escribir cada L^2 pasos
            }
        }
        // Medición
        for (int i = 0; i < M; i++) {
            for(int j = 0; j < L * L; j++) {
                ising(s, T, L, J, H, h); // Actualizar el sistema
            }
            magnetizaciones[i] = magnetizacion(s, L); // Calcular la magnetización
            salida_energia << energia_total(s, L, J, H) << "\n"; // Escribir la energía en el archivo
            escribir_datos(s, L, direccion); // Escribir cada paso de medida
        }

        // Calcular la media y la varianza de la magnetización
        double suma = 0, suma_cuadrados = 0;
        for (int sim = 0; sim < M; sim++) {
            suma += magnetizaciones[sim];
            suma_cuadrados += magnetizaciones[sim] * magnetizaciones[sim];
        }
        double media = suma / M;
        double varianza = (suma_cuadrados / M) - (media * media);

        medias_magnetizacion[j] = media;
        varianzas_magnetizacion[j] = varianza;

        // Calcular la susceptibilidad magnética
        susceptibilidades[j] = susceptibilidad_magnetica(magnetizaciones, M, T);

        cout << "  Temperatura " << T << " completada. Media magnetización: " << media << ", Varianza: " << varianza << endl;
    }

    // Escribir resultados en los archivos
    escribir_resultados(direccion_magnetizacion, temperaturas, medias_magnetizacion, varianzas_magnetizacion, n);

    ofstream salida_susceptibilidad(direccion_susceptibilidad);
    for (int i = 0; i < n; i++) {
        salida_susceptibilidad << temperaturas[i] << "," << susceptibilidades[i] << "\n";
    }
    salida_susceptibilidad.close();

    salida_energia.close(); // Cerrar el archivo de energía

    return 0;
}