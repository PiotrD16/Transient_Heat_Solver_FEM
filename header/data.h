#pragma once
#include <filesystem>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <limits>

#include "gaussLegendre.h"

#define ELEMENT_SIZE 4

using std::vector;
// dane wejsciowe

struct GlobalData {
    double simulationTime;       // całkowity czas symulacji [s]
    double simulationStepTime;   // krok czasowy delta t [s]
    double conductivity;         // przewodność cieplna k [W/(m*K)]
    double alfa;                 // współczynnik konwekcji alfa [W/(m^2*K)]
    double tOt;                  // temperatura otoczenia [K]
    double initialTemp;          // temperatura początkowa w całej siatce [K]
    double density;              // gęstość materiału [kg/m^3]
    double specificHeat;         // ciepło właściwe Cp [J/(kg*K)]
    int nN;                      // liczba węzłów globalnych
    int nE;                      // liczba elementów skończonych
    int npc;                     // liczba punktów całkowania w jednym kierunku
};

// wspolrzedne kazdego node'a

struct Node {
    double x, y;   // współrzędne węzła
    int BC;        // 1 = warunek brzegowy, 0 = brak

    Node(): x(0), y(0), BC(0) {};

    [[nodiscard]] double get_x () const{return x;}
    [[nodiscard]] double get_y () const{return y;}
    [[nodiscard]] double get_BC() const{return BC;}
};

// struktura Element zawiera numery punktow w rogach kazdego elementu
// czytane przeciwnie do ruchu wskazowek zegara

struct Jakobian {
    double J[2][2]{};      // macierz Jacobiego ∂(x,y)/∂(ksi,eta)
    double J1[2][2]{};     // macierz odwrotna do J → ∂(ksi,eta)/∂(x,y)
    double detJ;           // wyznacznik |J| = lokalna skalowana powierzchnia
    double dN_dx[4];       // pochodne funkcji kształtu w globalnym x
    double dN_dy[4];       // pochodne funkcji kształtu w globalnym y

    Jakobian() {
        J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        J1[0][0] = J1[0][1] = J1[1][0] = J1[1][1] = 0.0;
        detJ = 0.0;
    }
};

struct ElemUniv {
    vector<vector<double>> dN_dE;
    vector<vector<double>> dN_dn;

    explicit ElemUniv(const int npc) {
        dN_dE.resize(npc*npc, std::vector<double>(4));
        dN_dn.resize(npc*npc, std::vector<double>(4));
    }

    void calculateDerivatives(const Integral& integral);
    void displayDerivatives() const;
};

struct Surface {
    double ksi[2];
    double eta[2];
    double N[2][4]; // 2 punkty * 4 funkcje ksztaltu
};

struct Element {
    vector<int>element_id;  // [ELEMENT_SIZE]
    vector<Jakobian>Jakobis;// [npc]
    double HBC [ELEMENT_SIZE][ELEMENT_SIZE];
    Surface surfaces[ELEMENT_SIZE];
    double P[4];
    double C[4][4];
    int npc;

    Element() = default;

    explicit Element(const int npc): npc(npc) {
        element_id.resize(ELEMENT_SIZE);
        Jakobis.resize(npc*npc);

        auto nodes = Integral().getSchema2N().getNodes();

        // ŚCIANA 0 : (ksi, eta) = (ksi, -1)
        for (int pc=0; pc<npc; pc++){
            surfaces[0].ksi[pc] = nodes[pc];
            surfaces[0].eta[pc] = -1;
        }
        // ŚCIANA 1 : (ksi, eta) = (+1, eta)
        for (int pc=0; pc<npc; pc++){
            surfaces[1].ksi[pc] = 1;
            surfaces[1].eta[pc] = nodes[pc];
        }
        // ŚCIANA 2 : (ksi, eta) = (ksi, +1)
        for (int pc=0; pc<npc; pc++){
            surfaces[2].ksi[pc] = nodes[npc-1-pc]; // odwrócony kierunek
            surfaces[2].eta[pc] = 1;
        }
        // ŚCIANA 3 : (ksi, eta) = (-1, eta)
        for (int pc=0; pc<npc; pc++){
            surfaces[3].ksi[pc] = -1;
            surfaces[3].eta[pc] = nodes[npc-1-pc];
        }
    }

    // obliczanie jakobianu dla punktow calkowania
    void calculateJakobian(const vector<Node>& nodes, const ElemUniv& elemUniv);
    // obliczanie wyznacznikow dla jakobianow
    void calculateDetJ();
    // obliczanie odwrotnego jakobianu i liczymy pochodne dN/dX i dN/dY
    void calculateInverseJ(const ElemUniv& elemUniv);
    void displayJakobianData(int elementIndex) const;
    // obliczanie macierzy H - nie jest zapisywana bo potem jest potrzebna tylko do przeksztalcen
    void calculateMatrixH(
        const double& k,
        vector<vector<double>>& HGlobal,
        vector<vector<double>>& CGlobal,
        double density,
        double specificHeat);
    void calculateHbc(const std::vector<Node>& nodes,double alpha,double t_ot,std::vector<std::vector<double>>& HGlobal,std::vector<double>& PGlobal);
    void displayHbcAndP(int elementIndex) const;
};

struct SystemOfEquation {
    vector<double> t;     // wartosci temperatur
    vector<vector<double>> H;   // wektor H globalny
    vector<double> P;     // wektor P globalny
    vector<vector<double>> C;

    SystemOfEquation(const int nN, vector<double> P, vector<vector<double>> H) {
        t.reserve(nN);
        H.reserve(nN);
        P.reserve(nN);
        C.reserve(nN);

        this->P = P;
        this->H = H;
    }
    SystemOfEquation(const int nN, vector<double> P, vector<vector<double>> H, vector<vector<double>> C) {
        t.reserve(nN);
        H.reserve(nN);
        P.reserve(nN);
        C.reserve(nN);

        this->P = P;
        this->H = H;
        this->C = C;
    }

    // trzeba dodac jeszcze metode ktora bedzie obliczac wartosci temperatury (vector t)
    void solve();
    void solve2(const double stepTime); // tutaj rozwiazujemy uklad niestacjonarny
    [[nodiscard]] vector<double> getT () const {return t;}
};

struct Grid {
    vector<Element> element;
    vector<Node> nodes;

    Grid(const int nN,const int nE, const int npc) {
        element.reserve(nE);
        nodes.resize(nN);
        for (int i = 0; i < nE; i++)
            element.emplace_back(npc);
    }
};

bool readData(const std::string& filename, GlobalData &globalData, Grid* &grid);

void displayData(const GlobalData &globalData, Grid* &grid);