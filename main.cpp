#include <filesystem>
#include <iostream>
#include <fstream>
#include "header/data.h"
#include "header/gaussLegendre.h"

double f1(double x) {
    return  5*x*x + 3*x + 6;
}

double f2(double x, double y) {
    return 5*x*x*y*y + 3*x*y + 6;
}

int main() {

    GlobalData globalData{};
    Grid* grid = nullptr;

    std::string file1 = "Test1_4_4.txt";
    std::string file2 = "Test2_4_4_MixGrid.txt";
    std::string file3 = "Test3_31_31_kwadrat.txt";

    bool success = readData(file1, globalData, grid);
    if (!success) {
        std::cerr << "Error opening file " << std::endl;
        return -1;
    }
    displayData(globalData, grid);

    std::vector<std::vector<double>> HGlobal(globalData.nN, std::vector<double>(globalData.nN, 0.0));
    std::vector<double> PGlobal(globalData.nN, 0.0);
    ElemUniv elem_univ(globalData.npc);
    const Integral integral;
    elem_univ.calculateDerivatives(integral);
    elem_univ.displayDerivatives();

    // elem_univ zawiera juz obliczone pochodne czastkowe
    vector<vector<double>> CGlobal(globalData.nN, vector<double>(globalData.nN, 0.0));

    for (int e = 0; e < globalData.nE; ++e) {
        auto& elem = grid->element[e];
        elem.calculateJakobian(grid->nodes, elem_univ);
        elem.calculateDetJ();
        elem.calculateInverseJ(elem_univ);
        elem.calculateMatrixH(globalData.conductivity, HGlobal, CGlobal, globalData.density, globalData.specificHeat);

        elem.calculateHbc(grid->nodes, globalData.alfa, globalData.tOt, HGlobal, PGlobal);
        elem.displayHbcAndP(e);
    }

    // wyswietlenie macierzy H globalnej
    std::cout << "\nMacierz globalna H:\n";
    for (int i = 0; i < globalData.nN; ++i) {
        for (int j = 0; j < globalData.nN; ++j) {
            std::cout << std::setw(6) << std::setprecision(2) <<  HGlobal[i][j] << " ";
        }
        std::cout << "\n";
    }

    // wyswietlanie wketora P globalnego:
    std::cout << "\nWektor P globalny: \n";
    for (int i = 0; i < globalData.nN; i++) {
        std::cout << std::setw(10) << std::setprecision(5) << PGlobal[i] << " ";
    }

    // Do tego momentu wszystko liczone jest poprawnie: teraz wyznaczam wartosci temperatur po uprzednim
    // policzeniu H globalnego i P globalnego
    SystemOfEquation equation(globalData.nN, PGlobal, HGlobal);
    equation.solve();
    vector<double> solvedTemperatures = equation.getT();

    std::cout<<"\n";
    std::cout << "\nWartosci temperatur dla rozwiazania stacjonarnego: \n";
    for (auto i: solvedTemperatures) {
        std::cout << std::setw(10) << i;
    }
    std::cout << "\n";

    // teraz rozwiazujemy zadanie niestacjonarne: do równania dochodzi macierz C globalna
    // SystemOfEquation equation2(globalData.nN, PGlobal, HGlobal, CGlobal);
    // equation2.solve2(globalData.simulationStepTime);

    std::cout << "\nMacierz C globalna: \n";
    for (int i = 0; i < globalData.nN; ++i) {
        for (int j = 0; j < globalData.nN; ++j) {
            std::cout << std::setw(6) << std::setprecision(5) <<  CGlobal[i][j] << " ";
        }
        std::cout << "\n";
    }

    // teraz zalatwiam rozwiazanie niestacjonarne: macierze H, C i P sa policzone, pozostaje
    // rozwiazac ostatni uklad rownan, gdzie t0 w pierwszej iteracji jest rowne t1, tzn.
    // temperatura obliczona w poprzedniej iteracji rowna jest temperaturze poczatkowej w kolejnej iteracji

    SystemOfEquation equation2(globalData.nN, PGlobal, HGlobal, CGlobal);

    equation2.t.resize(globalData.nN);
    std::fill(equation2.t.begin(), equation2.t.end(), globalData.initialTemp);

    double currentTime = 0.0;
    while (currentTime < globalData.simulationTime) {
        currentTime += globalData.simulationStepTime;

        // Rozwiązanie dla bieżącego kroku (t0 -> t1)
        // Wynik zostanie nadpisany w equation.t
        equation2.solve2(globalData.simulationStepTime);

        // Wyświetlenie wyników min/max dla tego kroku
        double minT = equation2.t[0];
        double maxT = equation2.t[0];
        for (double val : equation2.t) {
            if (val < minT) minT = val;
            if (val > maxT) maxT = val;
        }

        std::cout << "Time: " << currentTime << " Min: " << minT << " Max: " << maxT << std::endl;
    }

    delete grid;

    // std::cout << "Calculating 1st function integral -> " << integral.calculateIntegralOneDim(f1, 2) << "\n";
    // std::cout << "Calculating 2nd function integral -> " << integral.calulateIntegralTwoDim(f2,4) << "\n";

    return 0;
}