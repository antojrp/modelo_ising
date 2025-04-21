# include <iostream>
# include <cstdlib> // Para rand() y srand()
# include <ctime>   // Para time()
# include <fstream>
# include <cmath>

using namespace std;
const int L=40; // Tamaño de la red
const int n=60; // Número de temperaturas
const int M0=1000; // Número de pasos de termalización
const int M = 50000; // Número de pasos de medida por temperatura

// Inicializar el generador de números aleatorios
void inicializar_aleatorio() {
    srand(time(nullptr)); // Semilla basada en el tiempo
}

void generar_datos(int s[L][L], int N)
{
    float p = 0.5; // Probabilidad de asignar -1
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < N; i++)
        {
            if ((rand() / (double)RAND_MAX) < p) // Ajustar la probabilidad para asignar -1
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

void ising(int s[L][L], int N, double h[5][2]) {
    int n = rand() % N; // Generar posición aleatoria
    int m = rand() % N;

    // Calcular la suma de los vecinos
    int suma_vecinos = s[n][(m + 1) % N] + s[n][(m + N - 1) % N] + s[(n + 1) % N][m] + s[(n + N - 1) % N][m];
    int idx = (s[n][m] * suma_vecinos + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
    int idy = (s[n][m] + 1) / 2; // Mapear el espín a índices (0, 1) 
    double p = h[idx][idy];
    double e = rand() / (double)RAND_MAX; // Generar probabilidad aleatoria

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
    inicializar_aleatorio(); // Inicializar el generador de números aleatorios

    int s[L][L];
    double T_i = 2.26; // Temperatura inicial
    double T_f = 2.32; // Temperatura final

    double J = 1; // Interacción entre espines
    double H = 0; // Campo magnético externo
    double rm = 0, rm2 = 0, rm4 = 0, rm0 = 0, rm1 = 0, c = 0, tau=0;

    string direccion_magnetizacion = "resultados/H=" + format_double(H) + "/ising_magnetizacion_L_" + to_string(L) + ".dat";
    string direccion_susceptibilidad = "resultados/H=" + format_double(H) + "/ising_susceptibilidad_L_" + to_string(L) + ".dat";
    string direccion_momento = "resultados/H=" + format_double(H) + "/ising_momento_L_" + to_string(L) + ".dat";
    string direccion_energia = "resultados/H=" + format_double(H) + "/ising_energia_L_" + to_string(L) + ".dat";
    ofstream salida_energia(direccion_energia); // Crear archivo para la energía

    double temperaturas[n];
    double medias_magnetizacion[n];
    double errores[n];
    double susceptibilidades[n];
    double momentos[n];

    generar_datos(s, L);
    
    for (int j = 0; j < n; j++) {
        double T = T_f - j * (T_f - T_i) / (n-1); // Temperatura decreciente
        double h[5][2];
        precalcular_factores(h, J, H, T); 
        temperaturas[j] = T;
        string direccion = "resultados/H=" + format_double(H) + "/ising_data_L_" + to_string(L) + "_T_" + format_double(T) + ".dat";
        // crear_fichero(direccion);


        cout << "Procesando temperatura " << j + 1 << " de " << n << " (T = " << T << ")" << endl;

        // escribir_datos(s, L, direccion);

        // Termalización
        for (int i = 0; i < M0 * L * L; i++) {
            ising(s, L, h);
            if (i % (L * L) == 0) {
                // escribir_datos(s, L, direccion); // Escribir cada L^2 pasos
            }
        }

        // Medición
        rm=0.0;
        rm2=0.0;
        rm4=0.0;
        rm1 = magnetizacion(s, L);
        c=0.0;
        tau=0.0;
        for (int i = 0; i < M; i++) {
            for(int j = 0; j < 10 * L * L; j++) {
                ising(s, L, h); // Actualizar el sistema
            }
            rm0 = magnetizacion(s, L); // Magnetización instantánea
            rm += rm0;
            rm2 += rm0 * rm0;
            rm4 += rm0 * rm0 * rm0 * rm0; 
            c += rm0 * rm1; 
            rm1 = rm0;

            salida_energia << energia_total(s, L, J, H) << "\n"; // Escribir la energía en el archivo
            // escribir_datos(s, L, direccion); // Escribir cada paso de medida
        }
        // escribir_datos(s, L, direccion);
        cout << "  Rm " << rm << " Rm2: "<< rm2 << endl;
        rm = rm / M;
        rm2 =rm2/M;
        rm4= 1 - (rm4/M) / (3 * rm2 * rm2);
        rm2 = rm2 - rm * rm; 
        c = (c / M - rm * rm) / rm2;
        cout << "  Correlación: " << c << " Rm2: "<< rm2 << endl;
        if (c != 1.0) {
            tau = c / (1.0 - c); // Tiempo de autocorrelación
        }

        double error = sqrt(rm2 * (2.0 * tau + 1.0) / M);

        medias_magnetizacion[j] = rm;
        errores[j] = error; // Guardamos el error, no la varianza
        susceptibilidades[j] = L*L*(rm2) / T;
        momentos[j] = rm4; // Momento de la magnetización

        cout << "  Temperatura " << T << " completada. Media magnetización: " << medias_magnetizacion[j] << ", error: " << error << endl;
    }

    // Escribir resultados en los archivos
    escribir_resultados(direccion_magnetizacion, temperaturas, medias_magnetizacion, errores, n);

    ofstream salida_susceptibilidad(direccion_susceptibilidad);
    for (int i = 0; i < n; i++) {
        salida_susceptibilidad << temperaturas[i] << "," << susceptibilidades[i] << "\n";
    }

    ofstream salida_momento(direccion_momento);
    for (int i = 0; i < n; i++) {
        salida_momento << temperaturas[i] << "," << momentos[i] << "\n";
    }


    salida_susceptibilidad.close();
    salida_momento.close(); 

    salida_energia.close(); // Cerrar el archivo de energía

    return 0;
}