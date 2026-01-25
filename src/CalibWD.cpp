#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>

#include <TROOT.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>

//I dati da usare sono nella cartella WavedumpDec10, non la carico su GitHub perché è troppo pesante

struct EventWave {
    int eventNumber = -1;
    std::vector<unsigned short> ch;
};

struct MultiChannelEvent {
    int eventNumber = -1;
    std::vector<double> sum;
};

std::map<int, EventWave> ReadWaveFile(const std::string &fname)
{
    std::ifstream in(fname);
    if (!in.is_open()) {
        std::cerr << "Impossibile aprire " << fname << std::endl;
        return {};
    }

    std::map<int, EventWave> events;
    std::string line;
    EventWave *currentWave = nullptr;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        if (line.find("Event Number") != std::string::npos) {
            size_t pos = line.find(':');
            if (pos != std::string::npos) {
                int ev = std::stoi(line.substr(pos + 1));
                EventWave &ew = events[ev];
                ew.eventNumber = ev;
                ew.ch.clear();
                currentWave = &ew;
            }
            continue;
        }

        if (!currentWave) continue;

        if ((line[0] >= '0' && line[0] <= '9') || line[0] == '-') {
            unsigned int adc;
            std::stringstream ss(line);
            ss >> adc;
            currentWave->ch.push_back(static_cast<unsigned short>(adc));
        }
    }

    return events;
}

void AnalyzeStep(const std::map<int, MultiChannelEvent> &merged, double thr = 400)
{
    if (merged.empty()) return;

        //Inizializzo istogramma

        TH1D *h = new TH1D("hAmp2", "Cosmic ray step amp", 100, 0, 2500);
        TGraph g;

        //Loop for prima importo il vettore delle somme contenuto in emrged (un MultiChannelEvent) e poi il numero di punti

        for (const auto &kv : merged) {
            const auto &wave = kv.second.sum;
            int N = wave.size();
            if (N < 50) continue;

            int j = 1;

            //Faccio la derivata discreta

            std::vector<double> dwave;

            while (j < N-2) {

                double d = (wave[j+1]-wave[j-1])/2;
                dwave.push_back(d);
                j++;

            }

            //Medio per evitare picchi spuri nella derivata

            std::vector<double> dwave_smooth(dwave.size(),0);
            
            for (size_t i=3; i < dwave.size()-3;i++){
                double sum = 0;
                for (int k=-1; k<=1; k++) {sum += dwave[i+k];};
                dwave_smooth[i] = std::abs(sum/3);
            }

            auto it_1 = std::max_element(dwave_smooth.begin(), dwave_smooth.end());

            double max_1 = *it_1;

            if (max_1 < thr) {continue;}

            h->Fill(max_1);
        }

    TCanvas *c = new TCanvas("cAmp2","AmpSpetrumCR",800,600);
    h->GetXaxis()->SetTitle("Amplitude [a.u.]");
    h->GetYaxis()->SetTitle("Counts [pure]");
    h->Draw();
    int maxBin = h->GetMaximumBin();  
    double xMax = h->GetBinCenter(maxBin); 
    cout << 60/xMax << endl;
}


void WaveAnalysis(const char* f0, const char* f1, const char* f2, const char* f3, double thr = 400)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    auto e0 = ReadWaveFile(f0);
    auto e1 = ReadWaveFile(f1);
    auto e2 = ReadWaveFile(f2);
    auto e3 = ReadWaveFile(f3);

    std::vector<double> calib = {13.7, 10.7, 10.4, 10.5}; 

    std::map<int, MultiChannelEvent> merged;

    for (const auto &kv : e0) {
        int ev = kv.first;
        if (!e1.count(ev) || !e2.count(ev) || !e3.count(ev)) continue;

        const auto &w0 = e0[ev].ch;
        const auto &w1 = e1[ev].ch;
        const auto &w2 = e2[ev].ch;
        const auto &w3 = e3[ev].ch;

        size_t N = std::min({w0.size(), w1.size(), w2.size(), w3.size()});
        MultiChannelEvent m;
        m.eventNumber = ev;
        m.sum.resize(N);

        for (size_t i=0; i<N; i++)
            m.sum[i] = w0[i] + w1[i] + w2[i] + w3[i];

        merged[ev] = std::move(m);
    }

    std::cout << "Eventi combinati: " << merged.size() << "\n";

    AnalyzeStep(merged, thr);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo totale: " 
              << std::chrono::duration<double>(t1-t0).count() 
              << " s\n";
}