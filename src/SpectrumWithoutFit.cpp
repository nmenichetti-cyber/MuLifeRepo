#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>

//I dati da usare sono nella cartella WavedumpFinal, non la carico su GitHub perché è troppo pesante

struct EventWave {
    int eventNumber = -1;
    std::vector<unsigned short> ch;
};

struct MultiChannelEvent {
    int eventNumber = -1;
    std::vector<double> sum;
};

bool CreateEvent(std::ifstream& in, EventWave& ev){

    ev.ch.clear();
    std::string line;

    // Numero associazione al numero dell'evento, se il file non è importato da false
    while (std::getline(in, line)) {
        if (line.find("Event Number") != std::string::npos) {
            ev.eventNumber = std::stoi(line.substr(line.find(':') + 1));
            break;
        }
    }
    if (!in) return false;

    // Legge i campioni finché non trova il prossimo evento
    std::streampos lastPos;
    while (true) {
        lastPos = in.tellg();

        //Se non trovi una riga testuale, fermati;

        if (!std::getline(in, line)) break;

        //Se trovi l'header, fermati;

        if (line.find("Event Number") != std::string::npos) {
            in.seekg(lastPos);  
            break;
        }

        //Se è vuota, fermati;

        if (line.empty()) continue;

        //Se c'è un numero, lo salvi nell'attributo ch;

        unsigned int adc;
        std::stringstream ss(line);
        if (ss >> adc)
            ev.ch.push_back(static_cast<unsigned short>(adc));
    }
    return true;
}


void AnalyzeEvent(const std::vector<double>& wave, TH1D* h, double thr = 200)
{
    int N = wave.size();

    static std::vector<double> dwave;
    static std::vector<double> dwave_smooth;

    dwave.clear();
    dwave.reserve(N);

    // Derivata discreta
    for (int j = 1; j < N - 1; ++j)
        dwave.push_back((wave[j + 1] - wave[j - 1]) / 2.0);

    dwave_smooth.assign(dwave.size(), 0.0);

    // Smoothing semplice, ripulisce da picchi spuri in derivata
    for (size_t i = 1; i + 1 < dwave.size(); ++i) {
        double s = dwave[i - 1] + dwave[i] + dwave[i + 1];
        dwave_smooth[i] = std::abs(s / 3.0);
    }

    // Primo fronte, cerco massimo assoluto
    auto it1 = std::max_element(dwave_smooth.begin(), dwave_smooth.end());
    if (*it1 < thr) return;

    int bl_ind = std::distance(dwave_smooth.begin(), it1);

    // Secondo fronte, cerco un massimo assoluto successivo
    int minDist = 20;
    if (bl_ind + minDist >= (int)dwave_smooth.size()) return;

    auto it2 = std::max_element(dwave_smooth.begin() + bl_ind + minDist, dwave_smooth.end());
    if (*it2 < thr) return;

    double c = 0.0564706;

    h->Fill((*it2) * c);
}

void WaveAnalysis(const char* f0name, const char* f1name, const char* f2name, const char* f3name, double thr = 200)
{
    std::ifstream f0(f0name), f1(f1name), f2(f2name), f3(f3name);

    //Controllo di import
    if (!f0 || !f1 || !f2 || !f3) {
        std::cerr << "Errore apertura file\n";
        return;
    }

    //Creo Evento e istogramma

    EventWave w0, w1, w2, w3;

    TH1D* h = new TH1D("hAmp2", "Electron energy", 100, 0, 70);

    int nev = 0;

    while ( CreateEvent(f0, w0) && CreateEvent(f1, w1) && CreateEvent(f2, w2) && CreateEvent(f3, w3) ) {

        size_t N = std::min({ w0.ch.size(), w1.ch.size(), w2.ch.size(), w3.ch.size() });

        static std::vector<double> sum;
        sum.resize(N);

        for (size_t i = 0; i < N; ++i) {sum[i] = w0.ch[i] + w1.ch[i] + w2.ch[i] + w3.ch[i];}

        AnalyzeEvent(sum, h, thr);
        ++nev;
    }

    std::cout << "Eventi analizzati: " << nev << "\n";

    TCanvas* c = new TCanvas("cAmp2", "Electron energy spectrum", 800, 600);
    h->GetXaxis()->SetTitle("Electron energy [MeV]");
    h->GetYaxis()->SetTitle("Counts");
    h->Draw();
}
