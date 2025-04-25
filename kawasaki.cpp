# include <iostream>
# include <random> // Librería para números aleatorios
# include <time.h>
# include <fstream>
# include <math.h>

using namespace std;
const int L=100; // Tamaño de la red
const int n=17; // Número de temperaturas
const int M0=12000; // Número de pasos de termalización
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

void clase(int s[L][L], int c[L][L], int i, int j) // Calcular la clase de los espines
{
    int vecinos_1 = (s[i][(j+1)%L]+s[i][(j+L-1)%L]+s[(i+1)%L][j]+s[(i+L-1)%L][j]+4)/2; // Vecinos con espín 1
    c[i][j] = 5*(1.5-s[i][j]/2.0)-vecinos_1; // Guardar la clase de espín
}

void cambiar_espin(int s[L][L], int c[L][L], int i, int j){
    int vecinos_contrario = ((-s[i][j])*(s[i][(j+1)%L]+s[i][(j+L-1)%L]+s[(i+1)%L][j]+s[(i+L-1)%L][j])+4)/2; // Vecinos con espín contrario
    std::uniform_int_distribution<int> dist(0.0, vecinos_contrario-1); // Distribución uniforme continua entre 0.0 y 1.0
    int random_value = dist(rng); // Generar un número aleatorio entre 0 y 1
    int contador=0;
    int k, l;
    int dx[] = {-1, 1, 0, 0};
    int dy[] = {0, 0, -1, 1};
    for(int d=0; d<4; d++)
    {
        k = dx[d];
        l = dy[d];
        if (s[i][j]!=s[(i+k+L)%L][(j+l+L)%L])
        {
            if(contador==random_value)
            {
                int aux=s[i][j]; // Guardar el espín
                s[i][j] = s[(i+k+L)%L][(j+l+L)%L]; // Cambiar el espín
                s[(i+k+L)%L][(j+l+L)%L] = aux; // Cambiar el espín
                break;
            }
            contador = contador +1;
        }
        
        
        
    }
    int m,o;
    for(int d=0; d<4; d++)
    {
        o = dx[d];
        m = dy[d];
            clase(s, c, (i+o+L)%L, (j+m+L)%L); 
            clase(s, c, (i+k+o+L)%L, (j+l+m+L)%L);
        
    }
}

double ising_kawasaki(int s[L][L], double h[5][2], int c[L][L]) {
    double suma=0;
    int contador[11]={0,0,0,0,0,0,0,0,0,0,0}; // Contador de clases de espín
    for (int j = 0; j < L; j++)
    {
        for (int i = 0; i < L; i++)
        {
            for (int k=1; k<=2; k++)
            {
                if ((s[(i + (k/2))%L][(j + (k%2))%L ]!= s[i][j])){
                    int suma_vecinos = s[i][(j + 1) % L] + s[i][(j + L - 1) % L] + s[(i + 1) % L][j] + s[(i + L - 1) % L][j];
                    int idx = (s[i][j]*suma_vecinos + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
                    int idy = (s[i][j] + 1) / 2; // Mapear el espín a índices (0, 1) 


                    int suma_vecinos2 = s[(i + (k/2))%L][(j + (k%2) + 1) % L] + s[(i + (k/2))%L][(j + (k%2) + L - 1) % L] + s[(i + (k/2) + 1) % L][(j + (k%2))%L ] + s[(i + (k/2) + L - 1) % L][(j + (k%2))%L ];
                    int idx2 = (s[(i + (k/2))%L][(j + (k%2))%L ]*suma_vecinos2 + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
                    int idy2 = (s[(i + (k/2))%L][(j + (k%2))%L ] + 1) / 2; // Mapear el espín a índices (0, 1)

                    suma += h[idx][idy]*h[idx2][idy2]/(h[4-idx][idy]*h[4-idx2][idy2]);
                }
            }
        }
    }



    std::uniform_real_distribution<double> dist(0.0, suma); // Distribución uniforme continua entre 0.0 y 1.0
    double random_value = dist(rng); // Generar un número aleatorio entre 0 y 1
    double cumulative_probability = 0.0;

    for (int j = 0; j < L; j++)
    {
        for (int i = 0; i < L; i++)
        {
            for (int k=1; k<=2; k++)
            {
                if ((s[(i + (k/2))%L][(j + (k%2))%L ]!= s[i][j])){
                    int suma_vecinos = s[i][(j + 1) % L] + s[i][(j + L - 1) % L] + s[(i + 1) % L][j] + s[(i + L - 1) % L][j];
                    int idx = (s[i][j]*suma_vecinos + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
                    int idy = (s[i][j] + 1) / 2; // Mapear el espín a índices (0, 1) 


                    int suma_vecinos2 = s[(i + (k/2))%L][(j + (k%2) + 1) % L] + s[(i + (k/2))%L][(j + (k%2) + L - 1) % L] + s[(i + (k/2) + 1) % L][(j + (k%2))%L ] + s[(i + (k/2) + L - 1) % L][(j + (k%2))%L ];
                    int idx2 = (s[(i + (k/2))%L][(j + (k%2))%L ]*suma_vecinos2 + 4) / 2; // Mapear delta_E (-4, -2, 0, 2, 4) a índices (0, 1, 2, 3, 4)
                    int idy2 = (s[(i + (k/2))%L][(j + (k%2))%L ] + 1) / 2; // Mapear el espín a índices (0, 1)


                    cumulative_probability +=  h[idx][idy]*h[idx2][idy2]/(h[4-idx][idy]*h[4-idx2][idy2]);

                    if (random_value < cumulative_probability) {
                        int aux = s[(i + (k/2))%L][(j + (k%2))%L ]; // Guardar el espín
                        s[(i + (k/2))%L][(j + (k%2))%L ] = s[i][j]; // Cambiar el espín
                        s[i][j] = aux; // Cambiar el espín
                        std::uniform_real_distribution<double> dist3(0.0, 1.0); // Distribución uniforme continua entre 0 y 1
                        double random_value3 = dist3(rng); // Generar un número aleatorio entre 0 y 1
                        return -L*L*log(random_value3)/suma;
                    }
                }
            }
        }
    }

    std::uniform_real_distribution<double> dist3(0.0, 1.0); // Distribución uniforme continua entre 0 y 1
    double random_value3 = dist3(rng); // Generar un número aleatorio entre 0 y 1
    return -L*L*log(random_value3)/suma;
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
    int s[L][L], c[L][L]; // Matriz de espines y matriz de clases
    double T_i = 2.2; // Temperatura inicial
    double T_f = 2.36; // Temperatura final
    double dt, t=0; // Tiempo y delta tiempo
    double J = 1; // Interacción entre espines
    double H = 0; // Campo magnético externo


    double temperaturas[n];
    double medias_magnetizacion[n];
    double varianzas_magnetizacion[n];
    double susceptibilidades[n];
    generar_datos(s, L);

    for (int j = 0; j < n; j++) {
        t=0; // Reiniciar el tiempo para cada temperatura
        double T = T_f - j * (T_f - T_i) / (n - 1); // Temperatura decreciente
        double h[5][2];
        precalcular_factores(h, J, H, T); // Precalcular factores de Boltzmann-Gibbs
        temperaturas[j] = T;
        generar_datos(s, L);
        // Crear un archivo específico para la energía de esta temperatura
        string direccion_energia = "resultados/Kawasaki/kawasaki_energia_L_" + to_string(L) + "_T_" + format_double(T) + ".dat";
        ofstream salida_energia(direccion_energia);

        string direccion = "resultados/Kawasaki/ising_data_L_" + to_string(L) + "_T_" + format_double(T) + ".dat";
        crear_fichero(direccion);

        cout << "Procesando temperatura " << j + 1 << " de " << n << " (T = " << T << ")" << endl;

        escribir_datos(s, L, direccion);

        for (int j = 0; j < L; j++) {
            for (int k = 0; k < L; k++) {
                clase(s, c, j, k); // Calcular la clase de los espines
            }
        }

        // Termalización
        for (int i = 0; i < M0; i++) {
            dt=ising_kawasaki(s, h, c); // Actualizar el sistema
            cout << "Paso de termalización: " << i << endl;
            // if (static_cast<int>((t) / (L * L)) != static_cast<int>(t / (L * L))) {
            if (i % (20) == 0) { // Cada L^2 pasos
                double energia = energia_total(s, L, J, H); // Calcular la energía total
                 // salida_energia << static_cast<int>((t + dt) / (L * L)) << "," << energia << "\n"; // Escribir la energía en el archivo
                salida_energia << i / (L * L) << "," << energia << "\n"; // Escribir la energía en el archivo
                escribir_datos(s, L, direccion); // Escribir cada L^2 pasos
            }
            t=t+dt; // Actualizar el tiempo

        }

        salida_energia.close(); // Cerrar el archivo de energía para esta temperatura
    }

    return 0;
}