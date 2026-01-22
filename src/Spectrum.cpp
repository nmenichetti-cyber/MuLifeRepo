// ============================================================================
//  WaveAnalysis_Optimized.C  (VERSIONE OTTIMIZZATA PER PILE-UP)
// ============================================================================

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

// ============================================================================
//  STRUTTURE DATI
// ============================================================================

struct EventWave {
    int eventNumber = -1;
    std::vector<unsigned short> ch;
};

struct MultiChannelEvent {
    int eventNumber = -1;
    std::vector<double> sum;  // ora double per unità fisiche
};

// ============================================================================
//  LETTURA FILE WAVEFORMS
// ============================================================================

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

// ============================================================================
//  ANALISI PILE-UP OTTIMIZZATA
// ============================================================================

void AnalyzePileupAmplitude(const std::map<int, MultiChannelEvent> &merged,
                            double minStepAmp = 0.2,  // in unità fisiche
                            double chiSquareThreshold = 0.7)
{
    if (merged.empty()) return;

    TH1D *hAmp2 = new TH1D("hAmp2", "Ampiezza secondo gradino; |amp2| [MeV]; Events",
                           100, 0, 100);

    TGraph g;
    TF1 fSingle("fSingle", "[0] + [1]/(1 + exp(-(x-[2])/[3]))");
    TF1 fDouble("fDouble",
                "[0] + [1]/(1 + exp(-(x-[2])/[3])) + [4]/(1 + exp(-(x-[5])/[6]))");

    for (const auto &kv : merged) {
        const auto &wave = kv.second.sum;
        int N = wave.size();
        if (N < 50) continue;

        // ---------- BASELINE ----------
        int nBase = std::min(50, N/4);
        double baseline = 0;
        for (int i=0; i<nBase; i++) baseline += wave[i];
        baseline /= nBase;

        int minVal = *std::min_element(wave.begin(), wave.end());
        int maxVal = *std::max_element(wave.begin(), wave.end());
        double ampInit = (std::abs(maxVal - baseline) > std::abs(minVal - baseline)) ?
                          maxVal - baseline : minVal - baseline;

        if (std::abs(ampInit) < minStepAmp) continue;

        // ---------- MIDPOINT ----------
        double mid = baseline + 0.5*ampInit;
        int iMid = 0;
        for (int i=0;i<N;i++) {
            if ((ampInit>0 && wave[i]>=mid) || (ampInit<0 && wave[i]<=mid)) { iMid=i; break; }
        }

        int iStart = std::max(0, iMid-30);
        int iEnd   = std::min(N-1, iMid+70);
        int nPts   = iEnd-iStart+1;

        // ---------- TGRAPH ----------
        g.Set(nPts);
        for (int i=iStart;i<=iEnd;i++) g.SetPoint(i-iStart, i-iStart, wave[i]);

        // ---------- FIT SINGOLO ----------
        fSingle.SetRange(0, nPts-1);
        fSingle.SetParameters(baseline, ampInit, iMid-iStart, 5.0);
        fSingle.SetParLimits(1, -2*std::abs(ampInit), 2*std::abs(ampInit));
        fSingle.SetParLimits(3, 0.5, 50.0);

        if(g.Fit(&fSingle, "QNR0")!=0) continue;
        double chi2s = fSingle.GetChisquare()/std::max(1,fSingle.GetNDF());
        if(chi2s<1.5) continue;

        // ---------- FIT DOPPIO ----------
        fDouble.SetRange(0,nPts-1);
        fDouble.SetParameters(
            fSingle.GetParameter(0),
            fSingle.GetParameter(1),
            fSingle.GetParameter(2),
            fSingle.GetParameter(3),
            0.5*ampInit,
            fSingle.GetParameter(2)+15,
            5.0
        );
        fDouble.SetParLimits(4, -2*std::abs(ampInit), 2*std::abs(ampInit));
        fDouble.SetParLimits(5, 5, nPts-1);
        fDouble.SetParLimits(6, 0.5, 50.0);

        if(g.Fit(&fDouble,"QNR0")!=0) continue;

        double chi2d = fDouble.GetChisquare()/std::max(1,fDouble.GetNDF());
        double amp2 = std::abs(fDouble.GetParameter(4));
        double sep  = fDouble.GetParameter(5) - fDouble.GetParameter(2);

        if(chi2d/chi2s < chiSquareThreshold && amp2>0.1*std::abs(ampInit) && sep>5 && sep<100) {
            hAmp2->Fill(amp2);
        }
    }

    TCanvas *c = new TCanvas("cAmp2","Electron energy spectrum",800,600);
    hAmp2->Draw();
    std::cout << "Mean = " << hAmp2->GetMean() << " RMS = " << hAmp2->GetRMS() << "\n";
}

// ============================================================================
//  MAIN ROOT FUNCTION
// ============================================================================

void WaveAnalysis(const char* f0,
                  const char* f1,
                  const char* f2,
                  const char* f3,
                  double minStepAmp = 0.2,       // MeV
                  double chiSquareThreshold = 0.7)
{
    auto t0 = std::chrono::high_resolution_clock::now();

    auto e0 = ReadWaveFile(f0);
    auto e1 = ReadWaveFile(f1);
    auto e2 = ReadWaveFile(f2);
    auto e3 = ReadWaveFile(f3);

    std::vector<double> calib = {13.7, 10.7, 10.4, 10.5}; // costanti canali

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
            m.sum[i] = (calib[0]*w0[i] + calib[1]*w1[i] + calib[2]*w2[i] + calib[3]*w3[i]) * 1e-3; // somma calibrata in MeV

        merged[ev] = std::move(m);
    }

    std::cout << "Eventi combinati: " << merged.size() << "\n";

    AnalyzePileupAmplitude(merged, minStepAmp, chiSquareThreshold);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo totale: " 
              << std::chrono::duration<double>(t1-t0).count() 
              << " s\n";
}

