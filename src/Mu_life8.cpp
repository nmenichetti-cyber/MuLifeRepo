/*****************************************************************************************
 *  Mu_life_new_v8.C
 *
 *  Macro ROOT per misura del lifetime del muone da stream FIFO (CH,CT) con parole reset.
 *
 *  Novità (v8):
 *   - Estrazione dt (START->STOP) UNA SOLA VOLTA (senza tagliare su tmin/tmax) e riuso.
 *   - Modalità SCAN automatica:
 *        * loop su configurazioni (tmin,tmax,nbins,modello)
 *        * fit Poisson likelihood (L) + integral per bin (I) nel range (R)
 *        * metriche: deviance/ndf, RMS(r_D), NLL (min FCN)
 *        * CSV: tmin,tmax,nbins,model,tau,etau,dev_over_ndf,rms_rD,nll
 *        * TGraphErrors e ROOT file di output + PNG riassuntivi
 *   - Test quantitativo bias START:
 *        * stima rate START (±sqrt(N)/T)
 *        * test Poisson su conteggi in finestre (index of dispersion -> p-value)
 *        * stima impatto su tau con modello di "competizione" (START come processo Poisson)
 *
 *  Uso consigliato (ROOT interattivo, robusto):
 *     root -l
 *     .L Mu_life_new_v8.C+
 *     Mu_life_new("FIFOread_Take5.txt", 300, 0.8, 20.0);   // analisi singola (come prima)
 *     Mu_life_scan("FIFOread_Take5.txt");                  // scan automatico
 *
 *  Nota:
 *   - NLL qui è salvato come MinFcnValue() dal TFitResult (FCN minimo).
 *     In fit "L" ROOT minimizza una funzione di likelihood (con costanti possibili),
 *     quindi NLL è un diagnostico comparativo, non un assoluto "statistico puro".
 *****************************************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>
#include <limits>
#include <sstream>
#include <iomanip>

// ROOT
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TPaveText.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TMath.h"

// =====================================================================
//                    COSTANTI HARDWARE / DECODIFICA
// =====================================================================

const unsigned int BIT_START = 1u;         // 1
const unsigned int BIT_STOP  = 1u << 1;    // 2  (STOP generale)
const unsigned int BIT_B8    = 1u << 2;    // 4
const unsigned int BIT_B9    = 1u << 3;    // 8
const unsigned int BIT_B10   = 1u << 4;    // 16
const unsigned int BIT_B11   = 1u << 5;    // 32

const unsigned int STOP_GENERIC_MASK = (BIT_STOP | BIT_B8 | BIT_B9 | BIT_B10 | BIT_B11);
const unsigned int BLOCK_MASK        = (BIT_B8 | BIT_B9 | BIT_B10 | BIT_B11);

const unsigned int RESET_FLAG        = (1u << 31);
const unsigned int COUNTER_MASK      = 0x3FFFFFFF;   // 30 bit bassi

// Tick e reset in microsecondi
const double tick_us    = 0.005;                             // 5 ns
const double reset_t_us = (double)(1ULL << 30) * tick_us;    // offset per ogni reset

// Parametri logici richiesti
const int    EARLY_STOP_MAX_TICKS  = 10;     // stop "immediato" entro 10 eventi dopo lo start
const double FINAL_STOP_MAX_US     = 20.0;   // stop fisico entro 20 µs dallo start
const int    EARLY_BLOCK_WINDOW    = 2;      // ±2 eventi per stimare i blocchi dello stop "immediato"
const int    FINAL_BLOCK_WINDOW    = 3;      // ±3 eventi per stimare i blocchi dello stop finale

// =====================================================================
//                      STRUTTURA EVENTO
// =====================================================================

struct Event {
    std::size_t  index;      // indice della riga nel file originale
    double       t_us;       // tempo assoluto [µs]
    unsigned int ch;         // channel word "piena"
    bool         isStart;
    bool         isStop;
    unsigned int stopMask;   // ch & STOP_GENERIC_MASK

    Event(std::size_t i, double t, unsigned int c)
        : index(i),
          t_us(t),
          ch(c),
          isStart((c & BIT_START) != 0u),
          isStop((c & STOP_GENERIC_MASK) != 0u),
          stopMask(c & STOP_GENERIC_MASK) {}
};

// =====================================================================
//                        FUNZIONI DI SUPPORTO
// =====================================================================

inline bool IsResetWord(unsigned int ch) {
    return (ch & RESET_FLAG) != 0u;
}

// OR dei bit dei blocchi (8,9,10,11) in una finestra di eventi
unsigned int CollectBlockMask(const std::vector<Event>& evs,
                              int centerIndex,
                              int halfWindow)
{
    unsigned int mask = 0u;
    int iMin = std::max(0, centerIndex - halfWindow);
    int iMax = std::min((int)evs.size() - 1, centerIndex + halfWindow);
    for (int i = iMin; i <= iMax; ++i) {
        mask |= (evs[i].ch & BLOCK_MASK);
    }
    return mask;
}

// ---------------------------------------------------------------------
// Deviance residuals per fit Poisson ("L") con opzione "I" coerente.
// ---------------------------------------------------------------------
TH1D* MakeDevianceResidualHist(const TH1* h, TF1* f,
                               double fitMin, double fitMax,
                               const char* name="hDevRes",
                               const char* title="Deviance residuals; t_{decay} [#mu s]; r_{D}")
{
    TH1D* hr = (TH1D*)h->Clone(name);
    hr->Reset("ICES");
    hr->SetTitle(title);

    const int nb = h->GetNbinsX();

    for (int i = 1; i <= nb; ++i) {
        const double x0 = h->GetXaxis()->GetBinLowEdge(i);
        const double x1 = h->GetXaxis()->GetBinUpEdge(i);

        if (x1 <= fitMin || x0 >= fitMax) continue;

        const double Ni = h->GetBinContent(i);
        const double binw = (x1 - x0);

        double mu = f->Integral(x0, x1) / binw;
        if (mu <= 0.0) mu = 1e-12;

        double di = 0.0;
        if (Ni > 0.0) {
            di = 2.0 * (Ni * std::log(Ni / mu) - (Ni - mu));
        } else {
            di = 2.0 * mu;
        }
        if (di < 0.0) di = 0.0;

        const double sgn = (Ni >= mu) ? 1.0 : -1.0;
        const double ri  = sgn * std::sqrt(di);

        hr->SetBinContent(i, ri);
        hr->SetBinError(i, 1.0);
    }

    return hr;
}

struct ResidQuality {
    int    nBinsUsed   = 0;
    double mean        = 0.0;
    double rms         = 0.0;
    int    nAbsGt3     = 0;
    int    nAbsGt5     = 0;
    double fracAbsGt3  = 0.0;
    double fracAbsGt5  = 0.0;
    double sumR2       = 0.0; // ~ devianza Poisson nel range
    int    ndf         = 0;
    double r2OverNdf   = 0.0;
};

ResidQuality ComputeResidualQuality(const TH1* hres, double fitMin, double fitMax, int nPar=0)
{
    ResidQuality q;
    const int nb = hres->GetNbinsX();

    double sum = 0.0;
    int n = 0;
    for (int i = 1; i <= nb; ++i) {
        const double x0 = hres->GetXaxis()->GetBinLowEdge(i);
        const double x1 = hres->GetXaxis()->GetBinUpEdge(i);
        if (x1 <= fitMin || x0 >= fitMax) continue;

        sum += hres->GetBinContent(i);
        n++;
    }

    q.nBinsUsed = n;
    if (n == 0) return q;

    q.mean = sum / (double)n;

    double sumsq = 0.0;
    int nGt3 = 0;
    int nGt5 = 0;

    for (int i = 1; i <= nb; ++i) {
        const double x0 = hres->GetXaxis()->GetBinLowEdge(i);
        const double x1 = hres->GetXaxis()->GetBinUpEdge(i);
        if (x1 <= fitMin || x0 >= fitMax) continue;

        const double r  = hres->GetBinContent(i);
        const double dr = r - q.mean;

        sumsq += dr * dr;

        if (std::fabs(r) > 3.0) nGt3++;
        if (std::fabs(r) > 5.0) nGt5++;

        q.sumR2 += r * r;
    }

    q.rms = std::sqrt(sumsq / (double)n);

    q.nAbsGt3 = nGt3;
    q.nAbsGt5 = nGt5;
    q.fracAbsGt3 = (double)nGt3 / (double)n;
    q.fracAbsGt5 = (double)nGt5 / (double)n;

    q.ndf = n - nPar;
    if (q.ndf < 1) q.ndf = 1;
    q.r2OverNdf = (q.ndf > 0) ? (q.sumR2 / (double)q.ndf) : 0.0;

    return q;
}

void PrintResidualQuality(const TH1* hres, double fitMin, double fitMax, int nPar, const char* label="")
{
    ResidQuality q = ComputeResidualQuality(hres, fitMin, fitMax, nPar);

    std::cout << "\n============= DEVIANCE RESIDUAL QUALITY " << label << " =============\n";
    std::cout << "Range:      [" << fitMin << ", " << fitMax << "] us\n";
    std::cout << "Bins used:  " << q.nBinsUsed << "\n";
    std::cout << "<r_D>:      " << q.mean << "\n";
    std::cout << "RMS(r_D):   " << q.rms  << "   (ideal ~ 1)\n";
    std::cout << "|r_D|>3:    " << q.nAbsGt3 << "  (" << 100.0*q.fracAbsGt3 << " %)\n";
    std::cout << "|r_D|>5:    " << q.nAbsGt5 << "  (" << 100.0*q.fracAbsGt5 << " %)\n";
    std::cout << "Sum(r_D^2): " << q.sumR2 << "   (~ Poisson deviance nel range)\n";
    std::cout << "ndf:        " << q.ndf << "   (bins used - nPar)\n";
    std::cout << "dev/ndf:    " << q.r2OverNdf << "\n";
    std::cout << "===================================================================\n";
}

// ---------------------------------------------------------------------
// Utility: istogramma 1D distribuzione residui (opzionale; qui lo lasciamo)
// ---------------------------------------------------------------------
TH1D* MakeResidualDistribution(const TH1* hres, double fitMin, double fitMax,
                               int nBins=60, double xmin=-6.0, double xmax=6.0,
                               const char* name="hDevResDist",
                               const char* title="Deviance residual distribution;r_{D};counts")
{
    TH1D* hdist = new TH1D(name, title, nBins, xmin, xmax);

    const int nb = hres->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
        const double x0 = hres->GetXaxis()->GetBinLowEdge(i);
        const double x1 = hres->GetXaxis()->GetBinUpEdge(i);
        if (x1 <= fitMin || x0 >= fitMax) continue;

        hdist->Fill(hres->GetBinContent(i));
    }

    return hdist;
}

// =====================================================================
//                 I/O: helper per prefisso output
// =====================================================================
std::string BasenameNoExt(const std::string& path)
{
    std::string s = path;
    // strip directory
    const std::size_t pos = s.find_last_of("/\\");
    if (pos != std::string::npos) s = s.substr(pos + 1);

    // strip extension
    const std::size_t dot = s.find_last_of('.');
    if (dot != std::string::npos) s = s.substr(0, dot);

    if (s.empty()) s = "mu_life";
    return s;
}

// =====================================================================
//                      ESTRAZIONE EVENTI + PAIRING
// =====================================================================

struct PairingStats {
    long long nStartSeen                 = 0;
    long long nStartDiscardNewStartEarly = 0;
    long long nStartDiscardNewStartLate  = 0;
    long long nStartNoEarlyStop          = 0;
    long long nStartNoFinalStop          = 0;
    long long nPairsBuilt                = 0;  // dt costruiti (prima di qualunque taglio tmin/tmax)
};

struct DataBundle {
    std::vector<Event>         events;
    std::vector<double>        startTimes_us;
    std::vector<double>        dt_all_us;
    std::vector<unsigned int>  stopBlocks_all;
    double                     liveTime_us = 0.0;
    long long                  nRowsRead   = 0;
    long long                  nResetsSeen = 0;
    PairingStats               pstats;
};

bool LoadEventsFromFile(const char* filename, DataBundle& db)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire il file " << filename << "\n";
        return false;
    }

    std::vector<unsigned int> CH;
    std::vector<unsigned int> CT;
    CH.reserve(200000);
    CT.reserve(200000);

    unsigned int ch_tmp = 0;
    unsigned int ct_tmp = 0;
    while (fin >> ch_tmp >> ct_tmp) {
        CH.push_back(ch_tmp);
        CT.push_back(ct_tmp);
    }
    fin.close();

    db.nRowsRead = (long long)CH.size();

    if (CH.empty() || CH.size() != CT.size()) {
        std::cerr << "[ERRORE] File vuoto o colonne di lunghezza diversa.\n";
        return false;
    }

    db.events.clear();
    db.events.reserve(CH.size());

    long long n_reset = -1;     // parte da -1, così il primo reset → offset 0
    bool seenFirstReset = false;

    for (std::size_t i = 0; i < CH.size(); ++i) {
        unsigned int ch = CH[i];
        unsigned int ct = CT[i];

        if (IsResetWord(ch)) {
            seenFirstReset = true;
            n_reset += 1;
            db.nResetsSeen++;
            continue;
        }

        if (!seenFirstReset) {
            continue; // buffer prima del primo reset
        }

        unsigned int ctr = (ct & COUNTER_MASK);
        double t_us_abs = (double)ctr * tick_us + (double)n_reset * reset_t_us;

        // Consideriamo solo eventi con almeno un bit significativo (START o STOP generico)
        unsigned int mask = ch & (BIT_START | STOP_GENERIC_MASK);
        if (mask == 0u) continue;

        db.events.emplace_back(i, t_us_abs, ch);
    }

    if (db.events.empty()) {
        std::cerr << "[ERRORE] Nessun evento utile dopo il primo reset.\n";
        return false;
    }

    // Livetime (approssimazione: primo->ultimo evento utile)
    db.liveTime_us = db.events.back().t_us - db.events.front().t_us;
    if (db.liveTime_us <= 0) db.liveTime_us = 0.0;

    // Estrai tempi START
    db.startTimes_us.clear();
    db.startTimes_us.reserve(db.events.size()/5);
    for (const auto& ev : db.events) {
        if (ev.isStart) db.startTimes_us.push_back(ev.t_us);
    }

    return true;
}

void BuildDecayPairs(const std::vector<Event>& events,
                     std::vector<double>& dt_all_us,
                     std::vector<unsigned int>& stopBlocks_all,
                     PairingStats& pst,
                     double maxDtUs = FINAL_STOP_MAX_US)
{
    dt_all_us.clear();
    stopBlocks_all.clear();
    pst = PairingStats{};

    dt_all_us.reserve(20000);
    stopBlocks_all.reserve(20000);

    const std::size_t N = events.size();
    std::size_t i = 0;

    while (i < N) {
        const Event& evStart = events[i];

        if (!evStart.isStart) {
            ++i;
            continue;
        }

        pst.nStartSeen++;

        const std::size_t idxStart = i;
        const double tStart = evStart.t_us;

        bool discardThisStart = false;

        // 1) Cerca stop "immediato" entro EARLY_STOP_MAX_TICKS eventi
        bool foundEarlyStop = false;
        std::size_t idxEarlyStop = 0;

        for (std::size_t j = idxStart + 1;
             j < N && j <= idxStart + (std::size_t)EARLY_STOP_MAX_TICKS;
             ++j) {

            const Event& ev2 = events[j];

            // Se nel mezzo appare un nuovo START → scartiamo quello vecchio (bias potenziale)
            if (ev2.isStart) {
                pst.nStartDiscardNewStartEarly++;
                i = j;
                discardThisStart = true;
                break;
            }

            if ((ev2.stopMask & STOP_GENERIC_MASK) != 0u) {
                foundEarlyStop = true;
                idxEarlyStop = j;
                break;
            }
        }

        if (discardThisStart) continue;

        if (!foundEarlyStop) {
            pst.nStartNoEarlyStop++;
            ++i;
            continue;
        }

        // 2) Cerca lo STOP FINALE (BIT_STOP) entro maxDtUs
        bool foundFinalStop = false;
        std::size_t idxFinalStop = 0;

        for (std::size_t j = idxEarlyStop + 1; j < N; ++j) {
            const Event& ev2 = events[j];

            const double dt_tmp = ev2.t_us - tStart;
            if (dt_tmp > maxDtUs) break;

            if (ev2.isStart) {
                pst.nStartDiscardNewStartLate++;
                i = j;
                discardThisStart = true;
                break;
            }

            if ((ev2.ch & BIT_STOP) != 0u) {
                foundFinalStop = true;
                idxFinalStop = j;
                break;
            }
        }

        if (discardThisStart) continue;

        if (!foundFinalStop) {
            pst.nStartNoFinalStop++;
            ++i;
            continue;
        }

        // 3) Coppia valida: salva dt (senza tagli su tmin/tmax) + maschera blocchi
        const Event& evStop = events[idxFinalStop];
        const double dt = evStop.t_us - tStart;

        if (dt >= 0.0 && dt <= maxDtUs) {
            const unsigned int finalBlockMask =
                CollectBlockMask(events, (int)idxFinalStop, FINAL_BLOCK_WINDOW);

            dt_all_us.push_back(dt);
            stopBlocks_all.push_back(finalBlockMask);
            pst.nPairsBuilt++;
        }

        i = idxFinalStop + 1;
    }
}

// =====================================================================
//                         BIAS START: ANALISI
// =====================================================================

struct StartProcessStats {
    long long nStart = 0;
    double liveTime_us = 0.0;

    double rStart_per_us  = 0.0;
    double erStart_per_us = 0.0;

    // Poisson window test (index of dispersion)
    double window_us = 1000.0;
    int    nWindows  = 0;
    double meanCount = 0.0;
    double varCount  = 0.0;
    double dispersionIndex = 0.0; // var/mean (Poisson ideal ~1)
    double chi2_disp = 0.0;
    int    ndf_disp  = 0;
    double pval_disp = 0.0;

    // Inter-arrival stats
    double meanDelta_us = 0.0;
    double rmsDelta_us  = 0.0;
};

StartProcessStats AnalyzeStartProcess(const std::vector<double>& startTimes_us,
                                     double liveTime_us,
                                     double window_us = 1.3e7)
{
    StartProcessStats s;
    s.nStart = (long long)startTimes_us.size();
    s.liveTime_us = liveTime_us;
    s.window_us = window_us;

    if (liveTime_us > 0.0 && s.nStart > 0) {
        s.rStart_per_us  = (double)s.nStart / liveTime_us;
        s.erStart_per_us = std::sqrt((double)s.nStart) / liveTime_us; // Poisson
    }

    // Inter-arrival
    if (startTimes_us.size() >= 2) {
        std::vector<double> d;
        d.reserve(startTimes_us.size()-1);
        for (std::size_t i = 0; i+1 < startTimes_us.size(); ++i) {
            double dt = startTimes_us[i+1] - startTimes_us[i];
            if (dt > 0) d.push_back(dt);
        }
        if (!d.empty()) {
            double sum = 0.0;
            for (double x : d) sum += x;
            s.meanDelta_us = sum / (double)d.size();

            double ss = 0.0;
            for (double x : d) ss += (x - s.meanDelta_us)*(x - s.meanDelta_us);
            s.rmsDelta_us = std::sqrt(ss / (double)d.size());
        }
    }

    // Poisson test su conteggi in finestre
    if (liveTime_us > 0.0 && window_us > 0.0) {
        const int nWin = (int)std::floor(liveTime_us / window_us);
        s.nWindows = nWin;

        if (nWin >= 5) {
            std::vector<int> counts(nWin, 0);

            const double t0 = startTimes_us.empty() ? 0.0 : startTimes_us.front();
            // Nota: trasliamo per stabilità numerica; conta in [t0, t0+liveTime]
            for (double ts : startTimes_us) {
                double dt = ts - t0;
                if (dt < 0) continue;
                int ib = (int)std::floor(dt / window_us);
                if (ib >= 0 && ib < nWin) counts[ib]++;
            }

            double mean = 0.0;
            for (int c : counts) mean += (double)c;
            mean /= (double)nWin;

            double var = 0.0;
            for (int c : counts) var += ( (double)c - mean )*( (double)c - mean );
            var /= (double)(nWin - 1);

            s.meanCount = mean;
            s.varCount  = var;

            if (mean > 0.0) {
                s.dispersionIndex = var / mean;
                s.ndf_disp = nWin - 1;
                s.chi2_disp = s.ndf_disp * s.dispersionIndex;

                // p-value: Prob(chi2, ndf)
                s.pval_disp = TMath::Prob(s.chi2_disp, s.ndf_disp);
            }
        }
    }

    return s;
}

void PrintStartBiasReport(const StartProcessStats& s, const PairingStats& pst,
                          double tau_fit_us = std::numeric_limits<double>::quiet_NaN(),
                          double etau_fit_us = std::numeric_limits<double>::quiet_NaN())
{
    std::cout << "\n==================== START PROCESS / BIAS TEST ====================\n";
    std::cout << "N_START (dopo primo reset): " << s.nStart << "\n";
    std::cout << "Livetime stimato:           " << s.liveTime_us << " us  ("
              << s.liveTime_us*1e-6 << " s)\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Rate START:                 " << s.rStart_per_us
              << " 1/us  = " << (s.rStart_per_us*1e3) << " kHz"
              << "   ± " << s.erStart_per_us << " 1/us\n";

    std::cout << "\n[Inter-arrival START]\n";
    std::cout << "  <Δt_START>:               " << s.meanDelta_us << " us\n";
    std::cout << "  RMS(Δt_START):            " << s.rmsDelta_us  << " us\n";
    if (s.meanDelta_us > 0) {
        std::cout << "  (Per Poisson ideale: Δt ~ exp, quindi RMS ~ mean)\n";
    }

    std::cout << "\n[Poisson window test: index of dispersion]\n";
    std::cout << "  Window size:              " << s.window_us << " us\n";
    std::cout << "  N windows:                " << s.nWindows << "\n";
    std::cout << "  Mean count/window:        " << s.meanCount << "\n";
    std::cout << "  Var  count/window:        " << s.varCount  << "\n";
    std::cout << "  Dispersion (var/mean):    " << s.dispersionIndex << "  (Poisson ideal ~1)\n";
    if (s.ndf_disp > 0) {
        std::cout << "  chi2 = ndf*(var/mean):     " << s.chi2_disp << "  ndf=" << s.ndf_disp
                  << "  p-value=" << s.pval_disp << "\n";
        std::cout << "  Interpretazione pratica:\n";
        std::cout << "   - p piccolo (~<1%): rate non stazionario / clustering / dead-time significativo\n";
        std::cout << "   - p ragionevole: compatibile con Poisson stazionario (al livello del test)\n";
    } else {
        std::cout << "  (Test non disponibile: troppo poche finestre o mean~0)\n";
    }

    std::cout << "\n[Impatto START sulla selezione (dal pairing)]\n";
    std::cout << "  Start visti:              " << pst.nStartSeen << "\n";
    std::cout << "  Scartati (new START early): " << pst.nStartDiscardNewStartEarly << "\n";
    std::cout << "  Scartati (new START late):  " << pst.nStartDiscardNewStartLate  << "\n";
    std::cout << "  No early stop:            " << pst.nStartNoEarlyStop << "\n";
    std::cout << "  No final stop:            " << pst.nStartNoFinalStop << "\n";
    std::cout << "  Coppie dt costruite:      " << pst.nPairsBuilt << "\n";

    // Modello di competizione: tau_true = 1 / (1/tau_fit - R)
    if (std::isfinite(tau_fit_us) && tau_fit_us > 0 && s.rStart_per_us > 0) {
        const double lambda_fit  = 1.0 / tau_fit_us;
        const double elambda_fit = (std::isfinite(etau_fit_us) && etau_fit_us > 0)
                                   ? (etau_fit_us / (tau_fit_us*tau_fit_us))
                                   : 0.0;

        const double lambda_true = lambda_fit - s.rStart_per_us;

        std::cout << "\n[Modello competizione: stima correzione]\n";
        std::cout << "  lambda_fit = 1/tau_fit:   " << lambda_fit << " 1/us\n";
        std::cout << "  R_START:                  " << s.rStart_per_us << " 1/us\n";

        if (lambda_true > 0) {
            const double tau_true = 1.0 / lambda_true;

            // Propagazione errori (assumendo indipendenza): sigma_lambda_true^2 = sigma_lambda_fit^2 + sigma_R^2
            const double elambda_true = std::sqrt(elambda_fit*elambda_fit + s.erStart_per_us*s.erStart_per_us);
            const double etau_true = elambda_true / (lambda_true*lambda_true);

            std::cout << "  tau_fit:                  " << tau_fit_us << " ± " << etau_fit_us << " us\n";
            std::cout << "  tau_true(stimato):        " << tau_true   << " ± " << etau_true   << " us\n";
            std::cout << "  Delta(tau):               " << (tau_true - tau_fit_us) << " us\n";
            std::cout << "  Nota critica:\n";
            std::cout << "   - Questa correzione è valida se gli START sono Poisson stazionari e\n";
            std::cout << "     la tua logica scarta eventi quando un nuovo START arriva prima dello STOP.\n";
        } else {
            std::cout << "  [WARN] 1/tau_fit - R <= 0: correzione non fisica (R troppo grande o tau_fit anomalo).\n";
        }
    }

    std::cout << "=====================================================================\n";
}

// =====================================================================
//                         FIT MODELS + DIAGNOSTICA
// =====================================================================

enum class FitModel {
    ExpBkg,           // [0]*exp(-x/[1]) + [2]
    ExpOnly,          // [0]*exp(-x/[1])
    ExpBkgCompete     // [0]*exp(-x*(1/[1] + R)) + [2]  (R fissato) -> [1] = tau_true
};

const char* FitModelName(FitModel m)
{
    switch (m) {
        case FitModel::ExpBkg:       return "ExpBkg";
        case FitModel::ExpOnly:      return "ExpOnly";
        case FitModel::ExpBkgCompete:return "ExpBkgCompete";
        default:                     return "Unknown";
    }
}

struct FitDiagnostics {
    FitModel model = FitModel::ExpBkg;

    double tau  = std::numeric_limits<double>::quiet_NaN();
    double etau = std::numeric_limits<double>::quiet_NaN();

    double devOverNdf = std::numeric_limits<double>::quiet_NaN();
    double rmsRD      = std::numeric_limits<double>::quiet_NaN();

    double nll = std::numeric_limits<double>::quiet_NaN(); // MinFcnValue()

    int nPar = 0;
    int status = -999;

    ResidQuality q;
    TF1* f = nullptr; // ownership esterno (per single-mode)
};

TF1* BuildFitFunction(FitModel model, const char* name, double fitMin, double fitMax, double Rstart_per_us)
{
    TF1* f = nullptr;

    if (model == FitModel::ExpBkg) {
        f = new TF1(name, "[0]*exp(-x/[1]) + [2]", fitMin, fitMax);
        f->SetParNames("N0", "tau", "B");
    }
    else if (model == FitModel::ExpOnly) {
        f = new TF1(name, "[0]*exp(-x/[1])", fitMin, fitMax);
        f->SetParNames("N0", "tau");
    }
    else if (model == FitModel::ExpBkgCompete) {
        // tau param = tau_true, R fissato
        // f(x) = N0 * exp(-x*(1/tau + R)) + B
        f = new TF1(name, "[0]*exp(-x*(1.0/[1] + [3])) + [2]", fitMin, fitMax);
        f->SetParNames("N0", "tau_true", "B", "Rstart");
        f->SetParameter(3, Rstart_per_us);
        f->FixParameter(3, Rstart_per_us);
    }

    return f;
}

void SeedInitialParams(TF1* f, TH1* h, FitModel model)
{
    // Stime iniziali robuste
    const double maxc = h->GetMaximum();
    if (model == FitModel::ExpOnly) {
        f->SetParameter(0, maxc);
        f->SetParameter(1, 2.2);
        f->SetParLimits(1, 0.2, 20.0);
    } else {
        f->SetParameter(0, maxc);
        f->SetParameter(1, 2.2);
        f->SetParLimits(1, 0.2, 20.0);

        // Fondo: media ultimi 10 bin
        double bkgGuess = 0.0;
        const int nb = h->GetNbinsX();
        const int nTail = std::min(10, nb);
        for (int ib = nb - nTail + 1; ib <= nb; ++ib) bkgGuess += h->GetBinContent(ib);
        bkgGuess /= (double)nTail;

        // Parametro B è [2] per ExpBkg e ExpBkgCompete
        f->SetParameter(2, bkgGuess);
        f->SetParLimits(2, 0.0, std::max(10.0, 10.0*bkgGuess + 100.0));
    }
}

FitDiagnostics FitAndDiagnose(TH1* h,
                              double fitMin, double fitMax,
                              FitModel model,
                              double Rstart_per_us = 0.0,
                              bool silent = true,
                              bool keepFunction = false) // keep TF1* for drawing in single mode
{
    FitDiagnostics out;
    out.model = model;
    out.nPar = (model == FitModel::ExpOnly) ? 2 : ((model == FitModel::ExpBkgCompete) ? 3 : 3);

    // Nome univoco funzione
    TString fname = Form("f_%s_%s", FitModelName(model), h->GetName());
    TF1* f = BuildFitFunction(model, fname.Data(), fitMin, fitMax, Rstart_per_us);
    if (!f) return out;

    SeedInitialParams(f, h, model);

    // Opzioni fit:
    //  L: likelihood Poisson
    //  I: confronto con integrale per bin (coerente con MakeDevianceResidualHist)
    //  R: range del TF1
    //  S: ritorna FitResult
    //  0: non disegnare
    TString opt = silent ? "0LIRS" : "LIRS";
    TFitResultPtr r = h->Fit(f, opt.Data());

    // Status e NLL
    out.status = (int)r; // in ROOT, r può convertire a int = status; ma usiamo anche r->Status se disponibile
    if (r.Get()) {
        out.status = r->Status();
        out.nll    = r->MinFcnValue();
    }

    // Tau
    if (model == FitModel::ExpOnly || model == FitModel::ExpBkg) {
        out.tau  = f->GetParameter(1);
        out.etau = f->GetParError(1);
    } else if (model == FitModel::ExpBkgCompete) {
        // tau_true è parametro [1]
        out.tau  = f->GetParameter(1);
        out.etau = f->GetParError(1);
    }

    // Deviance residuals + metriche
    TH1D* hDevRes = MakeDevianceResidualHist(h, f, fitMin, fitMax,
                                            Form("hDevRes_%s", h->GetName()),
                                            "Deviance residuals; t_{decay} [#mu s]; r_{D}");
    out.q = ComputeResidualQuality(hDevRes, fitMin, fitMax, out.nPar);
    out.devOverNdf = out.q.r2OverNdf;
    out.rmsRD      = out.q.rms;

    delete hDevRes;

    if (keepFunction) {
        out.f = f; // ownership passa al chiamante
    } else {
        delete f;
        out.f = nullptr;
    }

    return out;
}

// =====================================================================
//                              SCAN MODE
// =====================================================================

struct ScanRow {
    double tmin  = 0.0;
    double tmax  = 0.0;
    int    nbins = 0;
    FitModel model = FitModel::ExpBkg;

    double tau  = std::numeric_limits<double>::quiet_NaN();
    double etau = std::numeric_limits<double>::quiet_NaN();
    double devOverNdf = std::numeric_limits<double>::quiet_NaN();
    double rmsRD      = std::numeric_limits<double>::quiet_NaN();
    double nll        = std::numeric_limits<double>::quiet_NaN();
};

bool WriteScanCSV(const std::string& csvName, const std::vector<ScanRow>& rows)
{
    std::ofstream fout(csvName.c_str());
    if (!fout.is_open()) {
        std::cerr << "[ERRORE] Impossibile scrivere CSV: " << csvName << "\n";
        return false;
    }

    // Header ESATTAMENTE con le colonne richieste
    fout << "tmin,tmax,nbins,model,tau,etau,dev_over_ndf,rms_rD,nll\n";
    fout << std::fixed << std::setprecision(8);

    for (const auto& r : rows) {
        fout << r.tmin << ","
             << r.tmax << ","
             << r.nbins << ","
             << FitModelName(r.model) << ",";

        // Gestione NaN: scriviamo vuoto (più semplice in analisi successiva)
        auto writeNum = [&](double x){
            if (std::isfinite(x)) fout << x;
        };

        writeNum(r.tau);       fout << ",";
        writeNum(r.etau);      fout << ",";
        writeNum(r.devOverNdf);fout << ",";
        writeNum(r.rmsRD);     fout << ",";
        writeNum(r.nll);
        fout << "\n";
    }

    fout.close();
    return true;
}

void MakeAndSaveGraphs(const std::string& prefix,
                       const std::vector<ScanRow>& rows,
                       const StartProcessStats& sstats,
                       bool alsoMakeParamGraphs = true)
{
    // ROOT file per oggetti (TGraphErrors)
    std::string rootName = prefix + "_scan_graphs.root";
    TFile* froot = new TFile(rootName.c_str(), "RECREATE");

    // Graph tau vs scan index
    TGraphErrors* gTauVsIdx = new TGraphErrors();
    gTauVsIdx->SetName("gTauVsIndex");
    gTauVsIdx->SetTitle("Tau vs scan index;scan index;#tau [#mus]");

    // Graph tau vs tmin / tmax / nbins (scatter; possibili sovrapposizioni)
    TGraphErrors* gTauVsTmin = new TGraphErrors();
    TGraphErrors* gTauVsTmax = new TGraphErrors();
    TGraphErrors* gTauVsNbins= new TGraphErrors();

    gTauVsTmin->SetName("gTauVsTmin");
    gTauVsTmax->SetName("gTauVsTmax");
    gTauVsNbins->SetName("gTauVsNbins");

    gTauVsTmin->SetTitle("#tau vs t_{min};t_{min} [#mus];#tau [#mus]");
    gTauVsTmax->SetTitle("#tau vs t_{max};t_{max} [#mus];#tau [#mus]");
    gTauVsNbins->SetTitle("#tau vs nbins;nbins;#tau [#mus]");

    int ip = 0;
    int itmin = 0, itmax = 0, inb = 0;

    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& r = rows[i];
        if (!std::isfinite(r.tau) || !std::isfinite(r.etau)) continue;

        gTauVsIdx->SetPoint(ip, (double)ip, r.tau);
        gTauVsIdx->SetPointError(ip, 0.0, r.etau);
        ip++;

        if (alsoMakeParamGraphs) {
            gTauVsTmin->SetPoint(itmin, r.tmin, r.tau);
            gTauVsTmin->SetPointError(itmin, 0.0, r.etau);
            itmin++;

            gTauVsTmax->SetPoint(itmax, r.tmax, r.tau);
            gTauVsTmax->SetPointError(itmax, 0.0, r.etau);
            itmax++;

            gTauVsNbins->SetPoint(inb, (double)r.nbins, r.tau);
            gTauVsNbins->SetPointError(inb, 0.0, r.etau);
            inb++;
        }
    }

    // Scrivi su ROOT file
    froot->cd();
    gTauVsIdx->Write();
    gTauVsTmin->Write();
    gTauVsTmax->Write();
    gTauVsNbins->Write();
    froot->Close();
    delete froot;

    // Salva PNG (canvases)
    auto SaveGraphPNG = [&](TGraphErrors* g, const std::string& outPng){
        TCanvas* c = new TCanvas(Form("c_%s", g->GetName()), g->GetTitle(), 900, 700);
        g->SetMarkerStyle(20);
        g->Draw("AP");
        c->SaveAs(outPng.c_str());
        delete c;
    };

    SaveGraphPNG(gTauVsIdx,  prefix + "_tau_vs_index.png");
    SaveGraphPNG(gTauVsTmin, prefix + "_tau_vs_tmin.png");
    SaveGraphPNG(gTauVsTmax, prefix + "_tau_vs_tmax.png");
    SaveGraphPNG(gTauVsNbins,prefix + "_tau_vs_nbins.png");

    // Clean
    delete gTauVsIdx;
    delete gTauVsTmin;
    delete gTauVsTmax;
    delete gTauVsNbins;

    // (Opzionale) Salva un breve riassunto START in txt per tracciabilità
    std::ofstream ftxt((prefix + "_start_bias_summary.txt").c_str());
    if (ftxt.is_open()) {
        ftxt << std::fixed << std::setprecision(8);
        ftxt << "N_START=" << sstats.nStart << "\n";
        ftxt << "LIVE_us=" << sstats.liveTime_us << "\n";
        ftxt << "RSTART_per_us=" << sstats.rStart_per_us << "\n";
        ftxt << "eRSTART_per_us=" << sstats.erStart_per_us << "\n";
        ftxt << "window_us=" << sstats.window_us << "\n";
        ftxt << "dispersion=" << sstats.dispersionIndex << "\n";
        ftxt << "pval_disp=" << sstats.pval_disp << "\n";
        ftxt.close();
    }
}

// Entry point scan
void Mu_life_scan(const char* filename = "FIFOread_Take5.txt")
{
    std::cout << "\n============================================\n";
    std::cout << "[Mu_life_scan] File: " << filename << "\n";
    std::cout << "============================================\n";

    DataBundle db;
    if (!LoadEventsFromFile(filename, db)) return;

    std::cout << "[INFO] Righe lette:               " << db.nRowsRead << "\n";
    std::cout << "[INFO] Reset words viste:         " << db.nResetsSeen << "\n";
    std::cout << "[INFO] Eventi utili post-reset:   " << db.events.size() << "\n";
    std::cout << "[INFO] START estratti:            " << db.startTimes_us.size() << "\n";
    std::cout << "[INFO] Livetime stimato:          " << db.liveTime_us << " us\n";

    // Pairing una volta sola (dt fino a FINAL_STOP_MAX_US)
    BuildDecayPairs(db.events, db.dt_all_us, db.stopBlocks_all, db.pstats, FINAL_STOP_MAX_US);
    std::cout << "[INFO] dt costruiti (0.."<< FINAL_STOP_MAX_US <<" us): " << db.dt_all_us.size() << "\n";

    if (db.dt_all_us.empty()) {
        std::cerr << "[ATTENZIONE] Nessun dt ricostruito: controllare logica o acquisizione.\n";
        return;
    }

    // Analisi START (rate + test Poisson)
    StartProcessStats sstats = AnalyzeStartProcess(db.startTimes_us, db.liveTime_us, /*window_us=*/1.3e7);
    PrintStartBiasReport(sstats, db.pstats);

    // ---------------------------
    // Definizione griglia scan
    // ---------------------------
    // Nota: qui scegli tu la granularità. Evita griglie enormi se non necessario.
    std::vector<double> tminList  = {0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
    std::vector<double> tmaxList  = {12.0, 15.0, 18.0, 20.0};
    std::vector<int>    nbinsList = {150, 200, 250, 300, 400};

    // Modelli:
    //  - ExpBkg: standard (exp + costante)
    //  - ExpOnly: se vuoi testare sensibilità al fondo (spesso peggiora in coda)
    //  - ExpBkgCompete: tau interpretato come tau_true con R_START fissato (modello competizione)
    std::vector<FitModel> modelList = { FitModel::ExpBkg };

    // Se vuoi attivare anche la variante di competizione:
     //modelList = { FitModel::ExpBkg, FitModel::ExpBkgCompete };

    // Output names
    const std::string prefix = BasenameNoExt(filename);
    const std::string csvName = prefix + "_scan.csv";

    std::vector<ScanRow> rows;
    rows.reserve(tminList.size()*tmaxList.size()*nbinsList.size()*modelList.size());

    // Loop scan
    int cfg = 0;
    for (double tmin : tminList) {
        for (double tmax : tmaxList) {
            if (tmax <= tmin) continue;
            if (tmax > FINAL_STOP_MAX_US + 1e-9) continue;

            for (int nbins : nbinsList) {
                if (nbins < 20) continue;

                for (FitModel mdl : modelList) {

                    // Crea hist di questa configurazione e riempilo
                    TString hname = Form("hDecay_scan_%04d", cfg);
                    TH1F* h = new TH1F(hname.Data(),
                                       Form("Muon decay time (scan cfg %d); t_{decay} [#mu s]; Counts", cfg),
                                       nbins, tmin, tmax);
                    h->Sumw2();

                    for (double dt : db.dt_all_us) {
                        h->Fill(dt); // fuori range va in under/over, non pesa nel fit
                    }

                    // Fit + diagnostica (silent)
                    // Per ExpBkgCompete usa R_START in 1/us
                    const double R = sstats.rStart_per_us;
                    FitDiagnostics fd = FitAndDiagnose(h, tmin, tmax, mdl, R, /*silent=*/true, /*keepFunction=*/false);

                    ScanRow r;
                    r.tmin  = tmin;
                    r.tmax  = tmax;
                    r.nbins = nbins;
                    r.model = mdl;
                    r.tau   = fd.tau;
                    r.etau  = fd.etau;
                    r.devOverNdf = fd.devOverNdf;
                    r.rmsRD      = fd.rmsRD;
                    r.nll        = fd.nll;

                    rows.push_back(r);

                    // (Opzionale) output progress minimale
                    if ((cfg % 10) == 0) {
                        std::cout << "[SCAN] cfg=" << cfg
                                  << " t=["<<tmin<<","<<tmax<<"] nb="<<nbins
                                  << " model="<<FitModelName(mdl)
                                  << " tau="<<fd.tau<<" ± "<<fd.etau
                                  << " dev/ndf="<<fd.devOverNdf
                                  << " rmsRD="<<fd.rmsRD
                                  << "\n";
                    }

                    delete h;
                    cfg++;
                }
            }
        }
    }

    // Salva CSV
    if (WriteScanCSV(csvName, rows)) {
        std::cout << "[INFO] CSV scan salvato: " << csvName << "\n";
    }

    // Grafici + ROOT file
    MakeAndSaveGraphs(prefix, rows, sstats, /*alsoMakeParamGraphs=*/true);
    std::cout << "[INFO] Graphs salvati: " << prefix << "_tau_vs_*.png e " << prefix << "_scan_graphs.root\n";

    std::cout << "[INFO] Fine scan.\n";
}

// =====================================================================
//                          MU_LIFE_NEW (single)
// =====================================================================

void Mu_life_new(const char* filename = "FIFOread_Take5.txt",
                 int nbins = 300,
                 double tmin = 0.8,
                 double tmax = 20.0)
{
    std::cout << "\n============================================\n";
    std::cout << "[Mu_life_new] File: " << filename << "\n";
    std::cout << "[Mu_life_new] Range istogramma dt: [" << tmin << ", " << tmax << "] us\n";
    std::cout << "[Mu_life_new] Bins: " << nbins << "\n";
    std::cout << "============================================\n";

    DataBundle db;
    if (!LoadEventsFromFile(filename, db)) return;

    std::cout << "[INFO] Righe lette:               " << db.nRowsRead << "\n";
    std::cout << "[INFO] Reset words viste:         " << db.nResetsSeen << "\n";
    std::cout << "[INFO] Eventi utili post-reset:   " << db.events.size() << "\n";
    std::cout << "[INFO] START estratti:            " << db.startTimes_us.size() << "\n";
    std::cout << "[INFO] Livetime stimato:          " << db.liveTime_us << " us\n";

    // Pairing una sola volta
    BuildDecayPairs(db.events, db.dt_all_us, db.stopBlocks_all, db.pstats, FINAL_STOP_MAX_US);
    std::cout << "[INFO] dt costruiti (0.."<< FINAL_STOP_MAX_US <<" us): " << db.dt_all_us.size() << "\n";
    if (db.dt_all_us.empty()) {
        std::cerr << "[ATTENZIONE] Nessun dt ricostruito: controllare logica o parametri.\n";
        return;
    }

    // Analisi START (stampa + competizione verrà valutata dopo il fit)
    StartProcessStats sstats = AnalyzeStartProcess(db.startTimes_us, db.liveTime_us, /*window_us=*/1.3e7);

    // ------------------------------------------------------------
    // Statistiche combinazioni PMT di stop (blocchi 8–11)
    // ------------------------------------------------------------
    std::map<unsigned int, long long> comboCounts;
    for (std::size_t k = 0; k < db.stopBlocks_all.size(); ++k) {
        const unsigned int sb = db.stopBlocks_all[k];
        if (sb == 0u) continue;
        comboCounts[sb]++;
    }

    std::cout << "\n[INFO] Statistiche combinazioni PMT di stop (blocchi 8–11):\n";
    if (comboCounts.empty()) {
        std::cout << "  Nessun evento di stop con PMT del blocco attivo.\n";
    } else {
        for (const auto &kv : comboCounts) {
            const unsigned int mask = kv.first;
            const long long    n    = kv.second;

            std::cout << "  PMT blocco = { ";
            bool first = true;
            if (mask & BIT_B8)  { std::cout << (first ? "" : ", ") << "8";  first = false; }
            if (mask & BIT_B9)  { std::cout << (first ? "" : ", ") << "9";  first = false; }
            if (mask & BIT_B10) { std::cout << (first ? "" : ", ") << "10"; first = false; }
            if (mask & BIT_B11) { std::cout << (first ? "" : ", ") << "11"; first = false; }
            std::cout << " }  ->  " << n << " eventi"
                      << "   [mask = 0x" << std::hex << mask << std::dec << "]\n";
        }
    }
    std::cout << std::endl;

    // ------------------------------------------------------------
    // Istogrammi: totale + per PMT del blocco
    // ------------------------------------------------------------
    TH1F* hDecay = new TH1F("hDecay",
                           "Muon decay time; t_{decay} [#mu s]; Counts",
                           nbins, tmin, tmax);

    TH1F* hDecay_B8  = new TH1F("hDecay_B8",
                               "Muon decay time (stop PMT 8); t_{decay} [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B9  = new TH1F("hDecay_B9",
                               "Muon decay time (stop PMT 9); t_{decay} [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B10 = new TH1F("hDecay_B10",
                               "Muon decay time (stop PMT 10); t_{decay} [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B11 = new TH1F("hDecay_B11",
                               "Muon decay time (stop PMT 11); t_{decay} [#mu s]; Counts",
                               nbins, tmin, tmax);

    hDecay->Sumw2();
    hDecay_B8->Sumw2();
    hDecay_B9->Sumw2();
    hDecay_B10->Sumw2();
    hDecay_B11->Sumw2();

    for (std::size_t k = 0; k < db.dt_all_us.size(); ++k) {
        const double dt = db.dt_all_us[k];
        const unsigned int sb = db.stopBlocks_all[k];

        hDecay->Fill(dt);
        if (sb & BIT_B8)  hDecay_B8->Fill(dt);
        if (sb & BIT_B9)  hDecay_B9->Fill(dt);
        if (sb & BIT_B10) hDecay_B10->Fill(dt);
        if (sb & BIT_B11) hDecay_B11->Fill(dt);
    }

    std::cout << "[INFO] Entries istogramma totale: " << hDecay->GetEntries() << "\n";
    std::cout << "[INFO] Entries istogramma PMT8:  " << hDecay_B8->GetEntries() << "\n";
    std::cout << "[INFO] Entries istogramma PMT9:  " << hDecay_B9->GetEntries() << "\n";
    std::cout << "[INFO] Entries istogramma PMT10: " << hDecay_B10->GetEntries() << "\n";
    std::cout << "[INFO] Entries istogramma PMT11: " << hDecay_B11->GetEntries() << "\n";

    gStyle->SetOptFit(1);

    // ------------------------------------------------------------
    // Fit istogramma totale + deviance residuals + metriche
    // ------------------------------------------------------------
    FitDiagnostics fdTot = FitAndDiagnose(hDecay, tmin, tmax,
                                         FitModel::ExpBkg,
                                         /*Rstart*/ sstats.rStart_per_us,
                                         /*silent*/ true,
                                         /*keepFunction*/ true);

    TF1* fExpBkg = fdTot.f;

    std::cout << "\n================ RISULTATI FIT (TOTALE) ================\n";
    std::cout << "Model   = " << FitModelName(fdTot.model) << "\n";
    std::cout << "Tau      = " << fdTot.tau  << " ± " << fdTot.etau << " us\n";
    if (fExpBkg) {
        const double B  = fExpBkg->GetParameter(2);
        const double eB = fExpBkg->GetParError(2);
        std::cout << "B (fondo)= " << B    << " ± " << eB   << " counts/bin\n";
    }
    std::cout << "Deviance/ndf = " << fdTot.devOverNdf << "\n";
    std::cout << "RMS(r_D)     = " << fdTot.rmsRD << "\n";
    std::cout << "NLL (FCNmin) = " << fdTot.nll << "\n";
    std::cout << "Fit status   = " << fdTot.status << "\n";
    std::cout << "========================================================\n";

    // Stampa report START completo, includendo la correzione competizione (usa tau fit)
    PrintStartBiasReport(sstats, db.pstats, fdTot.tau, fdTot.etau);

    // Residui per plot
    TH1D* hDevRes = MakeDevianceResidualHist(hDecay, fExpBkg, tmin, tmax,
                              "hDevRes",
                              "Deviance residuals (total); t_{decay} [#mu s]; r_{D}");

    PrintResidualQuality(hDevRes, tmin, tmax, /*nPar=*/3, "(TOTAL)");

    // Distribuzione deviance residuals + fit gaussiano (opzionale ma utile)
    TH1D* hDevResDist = MakeResidualDistribution(hDevRes, tmin, tmax,
                                          60, -6, 6,
                                          "hDevResDist",
                                          "Deviance residual distribution (total);r_{D};counts");
    TF1* fGaus = new TF1("fDevResGaus", "gaus", -4, 4);
    hDevResDist->Fit(fGaus, "0RQ");

    // ------------------------------------------------------------
    // Canvas: sopra fit, sotto deviance residuals (totale)
    // ------------------------------------------------------------
    TCanvas* c1 = new TCanvas("c1", "Muon lifetime + deviance residuals (total)", 1000, 900);

    TPad* pad1 = new TPad("pad1","pad1", 0.0, 0.30, 1.0, 1.0);
    TPad* pad2 = new TPad("pad2","pad2", 0.0, 0.00, 1.0, 0.30);

    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.30);

    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    hDecay->SetMarkerStyle(20);
    hDecay->SetMarkerSize(1.0);
    hDecay->SetMarkerColor(kBlack);

    hDecay->GetXaxis()->SetLabelSize(0.0);
    hDecay->GetXaxis()->SetTitleSize(0.0);

    hDecay->Draw("E1");
    fExpBkg->Draw("same");

    // Box diagnostico coerente con "L": deviance/ndf e RMS residui
    const double height2 = 0.22;
    TPaveText* pGoF = new TPaveText(0.40, 0.80  , 0.62, 0.9 - height2 , "NDC");
    pGoF->SetFillStyle(0);
    pGoF->SetBorderSize(1);
    pGoF->SetTextAlign(12);
    pGoF->SetTextSize(0.03);
    pGoF->AddText(Form("Deviance/ndf = %.2f", fdTot.devOverNdf));
    pGoF->AddText(Form("<r_{D}> = %.3f  RMS = %.3f", fdTot.q.mean, fdTot.q.rms));
    pGoF->Draw();

    // Legenda informativa (range + binning)
    TLegend *leg = new TLegend(0.40, 0.80, 0.62, 0.9);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry((TObject*)0, Form("[%.1f, %.1f] #mus", tmin, tmax), "");
    leg->AddEntry((TObject*)0, Form("Bins: %d", nbins), "");
    leg->Draw();

    pad2->cd();
    hDevRes->SetStats(false);
    hDevRes->SetMinimum(-5);
    hDevRes->SetMaximum(5);

    hDevRes->GetYaxis()->SetNdivisions(505);
    hDevRes->GetYaxis()->SetTitleSize(0.10);
    hDevRes->GetYaxis()->SetLabelSize(0.08);
    hDevRes->GetYaxis()->SetTitleOffset(0.45);

    hDevRes->GetXaxis()->SetTitleSize(0.12);
    hDevRes->GetXaxis()->SetLabelSize(0.10);

    hDevRes->Draw("E1");

    TLine* l0  = new TLine(tmin,  0.0, tmax,  0.0);
    TLine* lp3 = new TLine(tmin,  3.0, tmax,  3.0);
    TLine* lm3 = new TLine(tmin, -3.0, tmax, -3.0);
    l0->SetLineStyle(2);
    lp3->SetLineStyle(3);
    lm3->SetLineStyle(3);
    l0->Draw("same");
    lp3->Draw("same");
    lm3->Draw("same");

    TString outName;
    outName.Form("mu_decay_fit_%.3f_%.3f_%d_bins_withdevres.png", tmin, tmax, nbins);
    c1->SaveAs(outName.Data());

    // Distribuzione residui
    TCanvas* cResDist = new TCanvas("cResDist", "Deviance residual distribution (total)", 700, 600);
    hDevResDist->Draw("E1");
    fGaus->Draw("same");
    TString outResName;
    outResName.Form("devres_distribution_%.3f_%.3f_%d_bins.png", tmin, tmax, nbins);
    cResDist->SaveAs(outResName.Data());

    // ------------------------------------------------------------
    // Subplot PMT: opzionali da terminale
    // ------------------------------------------------------------
    char choice_subplot = 'N';
    std::cout << "Vuoi vedere i subplot PMT8,9,10,11 (con fit e deviance residuals)? [S/N] ";
    std::cin  >> choice_subplot;

    if (choice_subplot == 'S' || choice_subplot == 's') {

        TH1F* hists[4] = {hDecay_B8, hDecay_B9, hDecay_B10, hDecay_B11};
        const char* labels[4] = {"PMT8", "PMT9", "PMT10", "PMT11"};

        for (int ih = 0; ih < 4; ++ih) {
            TH1F* hh = hists[ih];

            if (hh->GetEntries() < 50) {
                std::cout << "[WARN] " << labels[ih] << ": troppe poche entries (" << hh->GetEntries()
                          << "), salto fit e plot.\n";
                continue;
            }

            FitDiagnostics fd = FitAndDiagnose(hh, tmin, tmax,
                                               FitModel::ExpBkg,
                                               sstats.rStart_per_us,
                                               /*silent*/ true,
                                               /*keepFunction*/ true);

            TF1* ff = fd.f;

            TH1D* hp = MakeDevianceResidualHist(hh, ff, tmin, tmax,
                                    Form("hDevRes_%s", labels[ih]),
                                    Form("Deviance residuals (%s); t_{decay} [#mu s]; r_{D}", labels[ih]));

            PrintResidualQuality(hp, tmin, tmax, /*nPar=*/3, Form("(%s)", labels[ih]));

            TCanvas* c = new TCanvas(Form("c_%s", labels[ih]),
                                     Form("Muon lifetime + deviance residuals (%s)", labels[ih]),
                                     1000, 900);

            TPad* p1 = new TPad(Form("pad1_%s", labels[ih]), "pad1", 0.0, 0.30, 1.0, 1.0);
            TPad* p2 = new TPad(Form("pad2_%s", labels[ih]), "pad2", 0.0, 0.00, 1.0, 0.30);
            p1->SetBottomMargin(0.02);
            p2->SetTopMargin(0.02);
            p2->SetBottomMargin(0.30);
            p1->Draw();
            p2->Draw();

            p1->cd();
            hh->SetMarkerStyle(20);
            hh->SetMarkerSize(1.0);
            hh->SetMarkerColor(kBlack);
            hh->GetXaxis()->SetLabelSize(0.0);
            hh->GetXaxis()->SetTitleSize(0.0);
            hh->Draw("E1");
            ff->Draw("same");

            TLegend* l = new TLegend(0.40, 0.80, 0.62, 0.9);
            l->SetBorderSize(1);
            l->SetFillStyle(0);
            l->SetTextSize(0.03);
            l->AddEntry((TObject*)0, Form("%s  [%.1f, %.1f] #mus", labels[ih], tmin, tmax), "");
            l->AddEntry((TObject*)0, Form("Bins: %d", nbins), "");
            l->Draw();

            TPaveText* pDiagH = new TPaveText(0.40, 0.80  , 0.62, 0.9 - height2 , "NDC");
            pDiagH->SetFillStyle(0);
            pDiagH->SetBorderSize(1);
            pDiagH->SetTextAlign(12);
            pDiagH->SetTextSize(0.028);
            pDiagH->AddText(Form("Deviance/ndf = %.3f", fd.devOverNdf));
            pDiagH->AddText(Form("<r_D> = %.3f,  RMS = %.3f", fd.q.mean, fd.q.rms));
            pDiagH->Draw();

            p2->cd();
            hp->SetStats(false);
            hp->SetMinimum(-5);
            hp->SetMaximum(5);

            hp->GetYaxis()->SetNdivisions(505);
            hp->GetYaxis()->SetTitleSize(0.10);
            hp->GetYaxis()->SetLabelSize(0.08);
            hp->GetYaxis()->SetTitleOffset(0.45);

            hp->GetXaxis()->SetTitleSize(0.12);
            hp->GetXaxis()->SetLabelSize(0.10);

            hp->Draw("E1");

            TLine* z0  = new TLine(tmin,  0.0, tmax,  0.0);
            TLine* zp3 = new TLine(tmin,  3.0, tmax,  3.0);
            TLine* zm3 = new TLine(tmin, -3.0, tmax, -3.0);
            z0->SetLineStyle(2);
            zp3->SetLineStyle(3);
            zm3->SetLineStyle(3);
            z0->Draw("same");
            zp3->Draw("same");
            zm3->Draw("same");

            TString outPmt;
            outPmt.Form("mu_decay_fit_%s_%.3f_%.3f_%d_bins_withdevres.png", labels[ih], tmin, tmax, nbins);
            c->SaveAs(outPmt.Data());

            // Cleanup specifico
            delete hp;
            delete ff;
        }
    }

    std::cout << "[INFO] Fine.\n";

    // Cleanup totale (oggetti principali)
    delete fExpBkg;
    delete hDevRes;
    delete hDevResDist;
    delete fGaus;

    // Gli istogrammi li lasciamo vivi se vuoi ispezionarli in sessione;
    // se preferisci cleanup completo, decommenta:
    // delete hDecay; delete hDecay_B8; delete hDecay_B9; delete hDecay_B10; delete hDecay_B11;
}
