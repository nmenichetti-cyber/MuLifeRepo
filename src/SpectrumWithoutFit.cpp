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

void AnalyzeEnergies(const std::map<int, MultiChannelEvent> &merged, double thr = 1)
{
    if (merged.empty()) return;

    //Inizializzo istogramma

    TH1D *h = new TH1D("hAmp2", "Electron energy", 100, 0, 150);
    TGraph g;

    //Loop for prima importo il vettore delle somme contenuto in emrged (un MultiChannelEvent) e poi il numero di punti

    int counter = 0;

    for (const auto &kv : merged) {
        const auto &wave = kv.second.sum;
        int N = wave.size();
        if (N < 50) continue;

        int j = 1;

        //Faccio la derivata discreta

        vector<double> dwave;

        while (j < N-2) {

            double d = (wave[j+1]-wave[j-1])/2;
            dwave.push_back(d);
            j++;

        }

        //Medio per evitare picchi spuri nella derivata

        vector<double> dwave_smooth(dwave.size(),0);
        
        for (size_t i=3; i<dwave.size()-3;i++){
            double sum = 0;
            for (int k=-1; k<=1; k++) {sum += dwave[i+k];};
            dwave_smooth[i] = sum/3;
        }

        //Ora cerco i massimi, e salvo gli indici in un vettore

        vector<int> pks; 

        for (size_t l=1; l<dwave.size()-2; l++){
            if (std::abs(dwave_smooth[l])>thr && 
                std::abs(dwave_smooth[l]) > std::abs(dwave_smooth[l-1]) && 
                std::abs(dwave_smooth[l]) > abs(dwave_smooth[l+1])){pks.push_back(l+1);}
        }
        
        if (pks.size() >= 2){
            h->Fill(abs(wave[pks[1]]-wave[pks[0]]));
            counter++;
        }

                // --- Canvas per waveform originale ---
        //std::string canvasNameWave = "cWave_" + std::to_string(kv.first);
        //std::string canvasTitleWave = "Waveform Event " + std::to_string(kv.first);
        //TCanvas *cWave = new TCanvas(canvasNameWave.c_str(), canvasTitleWave.c_str(), 800, 400);

        //TGraph *grWave = new TGraph(wave.size(), wave.data());
        //grWave->SetLineColor(kBlack);
        //grWave->SetTitle(("Waveform Event " + std::to_string(kv.first)).c_str());
        //grWave->GetXaxis()->SetTitle("Campione");
        //grWave->GetYaxis()->SetTitle("ADC / MeV");
        //grWave->Draw("AL");

        // --- Canvas separato per derivata smussata ---
        //std::string canvasNameDeriv = "cDeriv_" + std::to_string(kv.first);
        //std::string canvasTitleDeriv = "Derivata Event " + std::to_string(kv.first);
        //TCanvas *cDeriv = new TCanvas(canvasNameDeriv.c_str(), canvasTitleDeriv.c_str(), 800, 400);

        //TGraph *grD = new TGraph(dwave_smooth.size(), dwave_smooth.data());
        //grD->SetLineColor(kRed);
        //grD->SetTitle(("Derivata Smussata Event " + std::to_string(kv.first)).c_str());
        //grD->GetXaxis()->SetTitle("Campione");
        //grD->GetYaxis()->SetTitle("d(ADC)/dC");
        //grD->Draw("AL");

        // --- opzionale: picchi sulla derivata ---
        //for (int idx : pks) {
        //    TLine *l = new TLine(idx, 0, idx, dwave_smooth[idx]);
        //    l->SetLineColor(kBlue);
        //    l->SetLineStyle(2);
        //    l->Draw();
        //}

        // Aggiorna canvas
        //cWave->Update();
        //cDeriv->Update();

        // --- Pausa interattiva ---
        //std::cout << "Premi Invio per vedere il prossimo evento...\n";
        //std::cin.get();

        // --- Chiudi canvas prima di passare al prossimo evento ---
        //cWave->Close();
        //cDeriv->Close();

        // --- opzionale: libera memoria per i grafici ---
        //delete grWave;
        //delete grD;
}

TCanvas *c = new TCanvas("cAmp2","Electron energy spectrum",800,600);
h->GetXaxis()->SetTitle("Electron energy [MeV]");
h->GetYaxis()->SetTitle("Counts [pure]");
h->Draw();
cout << counter << endl;

}

void WaveAnalysis(const char* f0, const char* f1, const char* f2, const char* f3, double thr = 0.4)
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
            m.sum[i] = (calib[0]*w0[i] + calib[1]*w1[i] + calib[2]*w2[i] + calib[3]*w3[i]) * 1e-3;

        merged[ev] = std::move(m);
    }

    std::cout << "Eventi combinati: " << merged.size() << "\n";

    AnalyzeEnergies(merged, thr);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo totale: " 
              << std::chrono::duration<double>(t1-t0).count() 
              << " s\n";
}