#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <numeric>
#include <algorithm>

// Parámetros del sistema
const int L = 40;           // Longitud del lado de la red (LxL)
const int N = L * L;        // Número total de espines
const int M = 8192;         // Número de pasos de medida por temperatura
const int M0 = 1000;        // Número de pasos de termalización
const int mc = 1;           // Número de pasos de Monte Carlo entre medidas

// Generadores de números aleatorios
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> dist(0.0, 1.0);
std::uniform_int_distribution<int> int_dist(0, N - 1);

// Devuelve un número aleatorio uniforme en [0,1)
double ran_u() {
    return dist(rng);
}

// Devuelve un índice aleatorio entre 0 y N-1
int i_ran(int N) {
    return int_dist(rng);
}

// Función que calcula los vecinos de cada sitio bajo condiciones periódicas
void neighbors(std::vector<int> &n1, std::vector<int> &n2, std::vector<int> &n3, std::vector<int> &n4, int L) {
    for (int iy = 1; iy <= L; ++iy) {
        for (int ix = 1; ix <= L; ++ix) {
            int i = (iy - 1) * L + (ix - 1); // Índice lineal

            // Vecino a la derecha
            int ix1 = ix + 1;
            if (ix1 == L + 1) ix1 = 1;
            n1[i] = (iy - 1) * L + (ix1 - 1);

            // Vecino arriba
            int iy2 = iy + 1;
            if (iy2 == L + 1) iy2 = 1;
            n2[i] = (iy2 - 1) * L + (ix - 1);

            // Vecino a la izquierda
            int ix3 = ix - 1;
            if (ix3 == 0) ix3 = L;
            n3[i] = (iy - 1) * L + (ix3 - 1);

            // Vecino abajo
            int iy4 = iy - 1;
            if (iy4 == 0) iy4 = L;
            n4[i] = (iy4 - 1) * L + (ix - 1);
        }
    }
}

int main() {
    // Arreglos de espines y vecinos
    std::vector<int> s(N), n1(N), n2(N), n3(N), n4(N);
    std::vector<double> h(9, 1.0); // Tabla de factores de Boltzmann (índices de -4 a +4)

    // Archivos de salida
    std::ofstream out_magnet("magnetization.dat"); // Magnetización en cada paso
    std::ofstream out_results("results.dat");      // Resultados finales por temperatura

    // Calcular los vecinos para cada sitio
    neighbors(n1, n2, n3, n4, L);

    // Configuración inicial aleatoria de espines
    for (int i = 0; i < N; ++i)
        s[i] = (ran_u() < 0.5) ? 1 : -1;

    // Bucle sobre la temperatura (de 4.0 a 0.1 en pasos de -0.1)
    for (double T = 2.36; T >= 2.2; T -= 0.01) {
        // Tabla de factores de Boltzmann-Gibbs para ΔE = 2J*s*Σvecinos
        for (int j = -4; j <= 4; j += 2)
            h[j + 4] = std::min(1.0, std::exp(-2.0 * j / T));

        // Termalización: se actualiza el sistema para alcanzar estado estacionario
        for (int ij = 0; ij < M0 * N; ++ij) {
            int i = i_ran(N); // Elegir espín aleatorio
            int ib = s[i] * (s[n1[i]] + s[n2[i]] + s[n3[i]] + s[n4[i]]); // Energía local
            if (ran_u() < h[ib + 4]) s[i] = -s[i]; // Regla de Metropolis
        }

        // Inicialización de observables
        double rm = 0.0, rm2 = 0.0, c = 0.0;
        double rm1 = std::abs(std::accumulate(s.begin(), s.end(), 0.0)) / N;

        // Bucle de medición
        for (int im = 0; im < M; ++im) {
            // M pasos de Monte Carlo por espín
            for (int ij = 0; ij < mc * N; ++ij) {
                int i = i_ran(N);
                int ib = s[i] * (s[n1[i]] + s[n2[i]] + s[n3[i]] + s[n4[i]]);
                if (ran_u() < h[ib + 4]) s[i] = -s[i];
            }

            // Medición de magnetización absoluta normalizada
            double rm0 = std::abs(std::accumulate(s.begin(), s.end(), 0.0)) / N;
            out_magnet << rm0 << "\n"; // Guardar magnetización

            // Acumulación de observables
            rm += rm0;
            rm2 += rm0 * rm0;
            c += rm0 * rm1;
            rm1 = rm0;
        }

        // Cálculo de valores promedios
        rm /= M;
        rm2 = rm2 / M - rm * rm;
        c = (c / M - rm * rm) / rm2;

        // Tiempo de correlación integrado
        double tau = (c != 1.0) ? c / (1.0 - c) : 0.0;

        // Error estadístico estimado
        double error = std::sqrt(rm2 * (2 * tau + 1) / M);

        // Guardar resultados finales para esta T
        out_results << T << " " << rm << " " << rm2 << " " << error << " " << mc * tau << " " << c << "\n";
    }

    // Cierre de archivos
    out_magnet.close();
    out_results.close();
    return 0;
}
