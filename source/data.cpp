//
// Created by piotr on 23.10.2025.
//

#include "../header/data.h"

using std::cout, std::endl;

bool readData(const std::string& filename, GlobalData &globalData, Grid* &grid) {
    std::ifstream file(filename);
    if ( !file.is_open()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return false;
    }
    std::string temp;

    file >> temp >> globalData.simulationTime;
    file >> temp >> globalData.simulationStepTime;
    file >> temp >> globalData.conductivity;
    file >> temp >> globalData.alfa;
    file >> temp >> globalData.tOt;
    file >> temp >> globalData.initialTemp;
    file >> temp >> globalData.density;
    file >> temp >> globalData.specificHeat;
    file >> temp >> temp >> globalData.nN;
    file >> temp >> temp >> globalData.nE;

    globalData.npc = 2;

    grid = new Grid(globalData.nN, globalData.nE, globalData.npc);


    std::string line;
    std::getline(file, line);            // consumes end of previous line
    std::getline(file, line);

    // Wczytywanie node'ow
    for (int i = 0; i < globalData.nN; i++) {
        int temp_int;
        file >> temp_int;
        file.ignore(1, ' ');
        file >> grid->nodes[i].x;
        file.ignore(1, ' ');
        file >> grid->nodes[i].y;
    }

    file >> temp >> temp;

    // Wczytywanie elementow
    for (int i = 0; i < globalData.nE; i++) {
        int temp_int;
        file >> temp_int;
        file.ignore(1, ',');
        for (int j = 0; j < ELEMENT_SIZE; j++) {
            file >> grid->element[i].element_id[j];
            if (j < ELEMENT_SIZE - 1) file.ignore(1, ',');
        }
    }

    std::string bcLabel;
    file >> bcLabel;

    if (bcLabel == "*BC") {

        // Pobierz resztę linii (usunięcie \n po *BC)
        std::string bcLine;
        std::getline(file, bcLine);

        // Teraz pobierz właściwą linię z numerami
        std::getline(file, bcLine);

        std::stringstream ss(bcLine);
        int bcNode;

        while (ss >> bcNode) {
            grid->nodes[bcNode - 1].BC = 1;
            if (ss.peek() == ',') ss.ignore();
        }
    }

        file.close();

        return true;
}

    void displayData(const GlobalData &globalData, Grid* &grid) {
        std::cout << "Dane globalne:\n";
        std::cout << "SimulationTime:     " << globalData.simulationTime     << "\n";
        std::cout << "SimulationStepTime: " << globalData.simulationStepTime << "\n";
        std::cout << "Conductivity:       " << globalData.conductivity       << "\n";
        std::cout << "Alfa:               " << globalData.alfa               << "\n";
        std::cout << "Tot (T otoczenia):  " << globalData.tOt                << "\n";
        std::cout << "InitialTemp:        " << globalData.initialTemp        << "\n";
        std::cout << "Density:            " << globalData.density            << "\n";
        std::cout << "SpecificHeat:       " << globalData.specificHeat       << "\n";
        std::cout << "Nodes number:       " << globalData.nN                 << "\n";
        std::cout << "Elements number:    " << globalData.nE                 << "\n";
        std::cout << "Punkty calkowania:  " << globalData.npc                 << "\n";

        std::cout << "Wspolrzedne wezlow:\n";
        for (int i = 0; i < globalData.nN; i++) {
            std::cout << "Node " << i + 1 << ": x=" << grid->nodes[i].x
                      << ", y=" << grid->nodes[i].y << " BC=" << grid->nodes[i].BC << "\n";
        }

        std::cout << "Elementy siatki:\n";
        for (int i = 0; i < globalData.nE; i++) {
            std::cout << "Element " << i + 1 << ": ";
            for (int j = 0; j < ELEMENT_SIZE; j++) {
                std::cout << grid->element[i].element_id[j];
                if (j < ELEMENT_SIZE - 1) std::cout << ", ";
            }
            std::cout << "\n";
        }
    }

    void ElemUniv::calculateDerivatives(const Integral& integral) {
        const auto schema = integral.getSchema2N();
        int totalPoints = dN_dE.size();
        const auto nodes = schema.getNodes();

        vector<std::pair<double,double>> gaussPoints(totalPoints);
        int index = 0;
        for (double ksi : nodes) {
            for (double eta : nodes) {
                if (index >= totalPoints) break;
                gaussPoints[index++] = {ksi, eta};
            }
        }

        for (int i = 0; i < totalPoints; i++) {

            double ksi = gaussPoints[i].first;
            double eta = gaussPoints[i].second;

            dN_dE[i][0] = -0.25 * (1 - eta); // dN1/dE
            dN_dE[i][1] =  0.25 * (1 - eta); // dN2/dE
            dN_dE[i][2] =  0.25 * (1 + eta); // dN3/dE
            dN_dE[i][3] = -0.25 * (1 + eta); // dN4/dE

            dN_dn[i][0] = -0.25 * (1 - ksi); // dN1/dN
            dN_dn[i][1] = -0.25 * (1 + ksi); // dN2/dN
            dN_dn[i][2] =  0.25 * (1 + ksi); // dN3/dN
            dN_dn[i][3] =  0.25 * (1 - ksi); // dN4/dN
        }
    }

    void ElemUniv::displayDerivatives() const {
        const int npc = dN_dE.size();

        cout << "Pochodne dN/dE:\n";
        for (int i = 0; i < npc; i++) {
            cout << "Node " << i + 1 << ": ";
            for (int j = 0; j < 4; j++) {
                cout << dN_dE[i][j] << " ";
            }
            cout << endl;
        }

        cout << "\nPochodne dN/dN:\n";
        for (int i = 0; i < npc; i++) {
            cout << "Node " << i + 1 << ": ";
            for (int j = 0; j < 4; j++) {
                cout << dN_dn[i][j] << " ";
            }
            cout << endl;
        }
    }

    void Element::calculateJakobian(const vector<Node>& nodes, const ElemUniv& elemUniv) {
        const int totalPoints = npc * npc;

        // macierz jakobiego powstaje w kazdym punkcie calkowania
        for (int i = 0; i < totalPoints; ++i) {

            Jakobis[i].J[0][0] = Jakobis[i].J[0][1] = Jakobis[i].J[1][0] = Jakobis[i].J[1][1] = 0.0;

            // iterujemy po 4 węzłach elementu
            for (int j = 0; j < ELEMENT_SIZE; ++j) {
                int nid = element_id[j] - 1; // numeracja w pliku zaczyna się od 1

                // dN/dksi * x_j
                Jakobis[i].J[0][0] += elemUniv.dN_dE[i][j] * nodes[nid].x;

                // dN/dksi * y_j
                Jakobis[i].J[0][1] += elemUniv.dN_dE[i][j] * nodes[nid].y;

                // dN/deta * x
                Jakobis[i].J[1][0] += elemUniv.dN_dn[i][j] * nodes[nid].x;

                // dN/deta * y
                Jakobis[i].J[1][1] += elemUniv.dN_dn[i][j] * nodes[nid].y;
            }
        }
    }

    void Element::calculateDetJ() {
        int totalPoints = npc * npc;
        for (int i = 0; i < totalPoints; ++i) {
            double det = Jakobis[i].J[0][0] * Jakobis[i].J[1][1] - Jakobis[i].J[0][1] * Jakobis[i].J[1][0];
            Jakobis[i].detJ = det;
        }
    }

    void Element::calculateInverseJ(const ElemUniv& elemUniv) {
        int totalPoints = npc * npc;
        for (int pc = 0; pc < totalPoints; ++pc) {
            double det = Jakobis[pc].detJ;
            if (fabs(det) < 1e-14) {
                std::cerr << "Warning: detJ ~ 0 in element, pc=" << pc << "\n";
                continue;
            }
            // odwrotna macierz J - rowniez liczona w kazdym punkcie calkowania
            // J^{-1} = (1/det) * [ J11  -J01 ; -J10  J00 ]
            Jakobis[pc].J1[0][0] =  Jakobis[pc].J[1][1] / det;   // dξ/dx
            Jakobis[pc].J1[0][1] = -Jakobis[pc].J[0][1] / det;   // dξ/dy
            Jakobis[pc].J1[1][0] = -Jakobis[pc].J[1][0] / det;   // dη/dx
            Jakobis[pc].J1[1][1] =  Jakobis[pc].J[0][0] / det;   // dη/dy

            // obliczamy globalne pochodne: dN/dx = dN/dksi * dksi/dx + dN/deta * deta/dx
            for (int j = 0; j < ELEMENT_SIZE; ++j) {

                Jakobis[pc].dN_dx[j] = Jakobis[pc].J1[0][0] * elemUniv.dN_dE[pc][j]
                                     + Jakobis[pc].J1[1][0] * elemUniv.dN_dn[pc][j]; // dN/dx = dN/dξ * dξ/dx + dN/dη * dη/dx
                Jakobis[pc].dN_dy[j] = Jakobis[pc].J1[0][1] * elemUniv.dN_dE[pc][j]
                                     + Jakobis[pc].J1[1][1] * elemUniv.dN_dn[pc][j]; // dN/dy = dN/dξ * dξ/dy + dN/dη * dη/dy
            }
        }
    }


    void Element::displayJakobianData(int elementIndex) const {
        cout << "Obliczenia macierzy Jakobiego dla elementu nr " << elementIndex + 1 << "\n";

        for (int pc = 0; pc < npc * npc; ++pc) {
            cout << "Macierz Jakobiego dla " << pc + 1 << " punktu calkowania\n";
            cout << Jakobis[pc].J[0][0] << " " << Jakobis[pc].J[0][1] << "\n";
            cout << Jakobis[pc].J[1][0] << " " << Jakobis[pc].J[1][1] << "\n";

            // Wyznacznik
            double det = Jakobis[pc].detJ;
            cout << "DetJ = " << det << " dla " << pc + 1 << " punktu calkowania\n";

            // Macierz odwrotna
            cout << "Macierz Jakobiego odwrotna dla " << pc + 1 << " punktu calkowania\n";
            cout << Jakobis[pc].J1[0][0] << " " << Jakobis[pc].J1[0][1] << "\n";
            cout << Jakobis[pc].J1[1][0] << " " << Jakobis[pc].J1[1][1] << "\n\n";

            cout << "Pochodne dN/dx i dN/dy:\n";
            for (int j = 0; j < ELEMENT_SIZE; ++j) {
                cout << "  N" << j + 1 << ": "
                     << "dN/dx = " << Jakobis[pc].dN_dx[j]
                     << ", dN/dy = " << Jakobis[pc].dN_dy[j] << "\n";
            }
            cout << endl;
        }
    }

    // kazdy element bedzie mial po 4 macierze H liczone w punktach calkowania i u nas beda mialy wymiar ELEMENT_SIZE * ELEMENT_SIZE
    void Element::calculateMatrixH(
        const double& k,
        vector<vector<double>>& HGlobal,
        vector<vector<double>>& CGlobal,
        double density,
        double specificHeat)
{
    Integral integral;
    auto weights = integral.getSchema2N().getWeights();
    auto nodes = integral.getSchema2N().getNodes();
    int totalPoints = npc * npc;

    double HL[ELEMENT_SIZE][ELEMENT_SIZE] = {0.0};
    double CL[ELEMENT_SIZE][ELEMENT_SIZE] = {0.0};   // <--- DODANA MACIERZ C

    std::vector<std::vector<std::vector<double>>> Hpc(
        totalPoints,
        std::vector<std::vector<double>>(ELEMENT_SIZE, std::vector<double>(ELEMENT_SIZE, 0.0))
    );

    // =================== OBRÓBKA H ===================
    for (int pc = 0; pc < totalPoints; ++pc) {
        double detJ = Jakobis[pc].detJ;

        for (int i = 0; i < ELEMENT_SIZE; i++) {
            for (int j = 0; j < ELEMENT_SIZE; j++) {

                Hpc[pc][i][j] =
                    k * (
                    Jakobis[pc].dN_dx[i] * Jakobis[pc].dN_dx[j] +
                    Jakobis[pc].dN_dy[i] * Jakobis[pc].dN_dy[j]
                ) * detJ;
            }
        }
    }

    // =================== OBRÓBKA C ===================
    for (int pc = 0; pc < totalPoints; ++pc) {

        double ksi = nodes[pc % npc];
        double eta = nodes[pc / npc];

        double detJ = Jakobis[pc].detJ;

        // FUNKCJE KSZTAŁTU N1..N4
        double N[4];
        N[0] = 0.25 * (1 - ksi) * (1 - eta);
        N[1] = 0.25 * (1 + ksi) * (1 - eta);
        N[2] = 0.25 * (1 + ksi) * (1 + eta);
        N[3] = 0.25 * (1 - ksi) * (1 + eta);

        double w = weights[pc % npc] * weights[pc / npc];

        for (int i = 0; i < ELEMENT_SIZE; i++) {
            for (int j = 0; j < ELEMENT_SIZE; j++) {

                CL[i][j] +=
                    density * specificHeat *
                    N[i] * N[j] *
                    detJ *
                    w;
            }
        }
    }

    // =================== SUMOWANIE H ===================
    for (int pc = 0; pc < totalPoints; ++pc) {

        int ksiIndex = pc % npc;
        int etaIndex = pc / npc;

        double w = weights[ksiIndex] * weights[etaIndex];

        for (int i = 0; i < ELEMENT_SIZE; ++i)
            for (int j = 0; j < ELEMENT_SIZE; ++j)
                HL[i][j] += Hpc[pc][i][j] * w;
    }

    // =================== WYPISANIE MACIERZY H ===================
    std::cout << "\nMacierz H dla elementu:\n";
    for (int i = 0; i < ELEMENT_SIZE; i++) {
        for (int j = 0; j < ELEMENT_SIZE; j++)
            std::cout << HL[i][j] << "\t";
        std::cout << "\n";
    }

    // =================== AGREGACJA H + C ===================
    for (int i = 0; i < ELEMENT_SIZE; ++i) {
        int id_i = element_id[i] - 1;

        for (int j = 0; j < ELEMENT_SIZE; ++j) {
            int id_j = element_id[j] - 1;

            HGlobal[id_i][id_j] += HL[i][j];
            CGlobal[id_i][id_j] += CL[i][j];
        }
    }
}


void Element::calculateHbc(const std::vector<Node>& nodes,double alpha,double t_ot,std::vector<std::vector<double>>& HGlobal,std::vector<double>& PGlobal) {

    const Integral integral;
    const auto gn = integral.getSchema2N().getNodes();   // gn od nodes bo byl jakis konflikt nazw
    const auto weights = integral.getSchema2N().getWeights();

    // wyzerowanie tablic zeby nie bylo tam smieci
    for (int i = 0; i < ELEMENT_SIZE; ++i) {
        P[i] = 0.0;
        for (int j = 0; j < ELEMENT_SIZE; ++j)
            HBC[i][j] = 0.0;
    }

    // iterowanie po 4 bokach elementu:
    for (int side = 0; side < 4; ++side) {
        int id1 = element_id[side] - 1;
        int id2 = element_id[(side + 1) % 4] - 1;

        // bok należy traktować jako brzegowy tylko gdy oba węzły mają BC==1
        if (!(nodes[id1].BC == 1 && nodes[id2].BC == 1))
            continue;

        // długość boku (J dla transformacji 1D) -> długość/2
        double dx = nodes[id2].x - nodes[id1].x;
        double dy = nodes[id2].y - nodes[id1].y;
        double detJ_edge = std::sqrt(dx*dx + dy*dy) / 2.0;

        // Dla każdego punktu Gaussa na krawędzi
        for (int g = 0; g < 2; ++g) {
            double ksi = gn[g];
            double weight = weights[g];

            // funkcje kształtu 1D na krawędzi (liniowe)
            double Nedge[2] = { 0.5 * (1.0 - ksi), 0.5 * (1.0 + ksi) };

            // odwzorowuje na funkcje 4-nodowego elementu (reszta = 0)
            double N[4] = {0.0, 0.0, 0.0, 0.0};
            N[side] = Nedge[0];
            N[(side + 1) % 4] = Nedge[1];

            // obliczanie lokalnego Hbc i P
            for (int i = 0; i < ELEMENT_SIZE; ++i) {
                for (int j = 0; j < ELEMENT_SIZE; ++j) {
                    HBC[i][j] += alpha * N[i] * N[j] * detJ_edge * weight;
                }
                P[i] += alpha * N[i] * t_ot * detJ_edge * weight ;
            }
        }
    }

    // Agregacja wektora P (wektor globalny jest w mainie i jest przekazywany do tej struktury):
    for (int i = 0; i < ELEMENT_SIZE; ++i) {
        int I = element_id[i] - 1;
        if (I < 0) continue;
        PGlobal[I] += P[i];

        // dodawanie macierzy HBC do H globalnej
        for(int j=0;j<4;j++){
            int J = element_id[j] - 1;
            HGlobal[I][J] += HBC[i][j];
        }
    }
}

void Element::displayHbcAndP(int elementIndex) const {
    std::cout << "\n===== ELEMENT " << elementIndex + 1 << " =====\n";

    std::cout << "Macierz HBC dla elementu " << elementIndex + 1 << ":\n";
    for (int i = 0; i < ELEMENT_SIZE; ++i) {
        for (int j = 0; j < ELEMENT_SIZE; ++j)
            std::cout << std::setw(12) << HBC[i][j] << " ";
        std::cout << "\n";
    }

    std::cout << "Wektor P dla elementu " << elementIndex + 1 << ":\n";
    for (int i = 0; i < ELEMENT_SIZE; ++i)
        std::cout << "P[" << i + 1 << "] = " << P[i] << "\n";
}

void SystemOfEquation::solve() {
    int n = H.size();

    // Kopia macierzy H do macierzy A (żeby nie niszczyć oryginału, jeśli potrzebny później)
    vector<vector<double>> A = H;
    vector<double> b(n);

    // Przekształcenie równania: H * t + P = 0  =>  H * t = -P
    // Tworzymy wektor wyrazów wolnych b
    for (int i = 0; i < n; i++) {
        b[i] = P[i];
    }

    // --- Eliminacja Gaussa z Częściowym Wyborem Elementu Głównego (Partial Pivoting) ---
    for (int k = 0; k < n - 1; k++) {

        // 1. KROK: Wybór pivota (szukamy największej wartości w kolumnie k poniżej przekątnej)
        int maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[maxRow][k])) {
                maxRow = i;
            }
        }

        // Zamiana wierszy (macierzy A i wektora b), jeśli znaleziono lepszy pivot
        if (maxRow != k) {
            std::swap(A[k], A[maxRow]);
            std::swap(b[k], b[maxRow]);
        }

        // Sprawdzenie czy macierz nie jest osobliwa
        if (std::abs(A[k][k]) < 1e-14) {
            std::cerr << "Blad: Macierz osobliwa lub zbyt bliska zeru na przekatnej w wierszu " << k << "\n";
            return; // Można rzucić wyjątek
        }

        // 2. KROK: Eliminacja wierszy poniżej k
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];

            // Zerujemy element pod przekątną (dla porządku, matematycznie jest 0)
            A[i][k] = 0.0;

            // Odejmujemy przeskalowany wiersz k od wiersza i
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            // To samo dla wektora wyników
            b[i] -= factor * b[k];
        }
    }

    // --- Podstawianie Wsteczne (Back Substitution) ---
    t.resize(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        // Suma znanych już wartości (na prawo od przekątnej)
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * t[j];
        }
        // Wyliczenie niewiadomej t[i]
        t[i] = (b[i] - sum) / A[i][i];
    }
}

void SystemOfEquation::solve2(const double stepTime) {
    int n = H.size();

    // Upewnij się, że wektor t ma odpowiedni rozmiar (jest to t0 dla tego kroku).
    // Jeśli to pierwsza iteracja, t musi być wcześniej zainicjowane temperaturą początkową.
    if (t.size() != n) {
        t.resize(n, 0.0); // Lub inna bezpieczna inicjalizacja, jeśli t było puste
    }

    // Lokalne kopie macierzy efektywnej A i wektora b
    // A reprezentuje lewą stronę równania: [H] + [C]/dT
    vector<vector<double>> A(n, vector<double>(n));
    // b reprezentuje prawą stronę równania: {P} + ([C]/dT)*{t0}
    vector<double> b(n);

    // 1. Budowa macierzy efektywnej A oraz wektora prawej strony b
    for (int i = 0; i < n; ++i) {
        double CdT_times_t0_i = 0.0; // Składnik ([C]/dT * t0) dla wiersza i

        for (int j = 0; j < n; ++j) {
            double value_CdT = C[i][j] / stepTime;

            // Element macierzy efektywnej A = H + C/dT
            A[i][j] = H[i][j] + value_CdT;

            // Obliczanie części wektora obciążeń pochodzącej od pojemności cieplnej
            // Mnożymy wiersz macierzy C/dT przez wektor temperatur z poprzedniego kroku (t[j])
            CdT_times_t0_i += value_CdT * t[j];
        }

        // Konstrukcja wektora wyrazów wolnych b
        // b = P + (C/dT)*t0
        // Zakładamy, że P jest dodatnie (ciepło wchodzące), zgodnie z poprawką w poprzedniej części rozmowy
        b[i] = P[i] + CdT_times_t0_i;
    }

    // 2. Rozwiązanie układu równań A * t_nowe = b metodą Eliminacji Gaussa
    // (Kod analogiczny do metody solve(), ale operujący na macierzy A i wektorze b)

    for (int k = 0; k < n - 1; k++) {
        // Wybór pivota (Partial Pivoting)
        int maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[maxRow][k])) {
                maxRow = i;
            }
        }

        if (maxRow != k) {
            std::swap(A[k], A[maxRow]);
            std::swap(b[k], b[maxRow]);
        }

        if (std::abs(A[k][k]) < 1e-14) continue; // Zabezpieczenie przed dzieleniem przez 0

        // Eliminacja
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            A[i][k] = 0.0;
            for (int j = k + 1; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // 3. Podstawianie wsteczne (Back Substitution)
    // Wynik zapisujemy bezpośrednio do wektora t, aktualizując stan układu na nowy krok czasowy
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * t[j];
        }
        t[i] = (b[i] - sum) / A[i][i];
    }
}