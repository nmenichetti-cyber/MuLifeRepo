#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <map>

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
#include "TAxis.h"

// =====================================================================
//                    COSTANTI HARDWARE / DECODIFICA
// =====================================================================
//
// Colonna 1 (CH)  : channel word
// Colonna 2 (CT)  : counter word
//
// Codifica canali (CH):
//   bit0 (1)  : START
//   bit1 (2)  : STOP generale (uscita dual timer stop)
//   bit2 (4)  : PMT8  & gate
//   bit3 (8)  : PMT9  & gate
//   bit4 (16) : PMT10 & gate
//   bit5 (32) : PMT11 & gate
//   bit31     : parola di reset del contatore (2^31)
//
// Il counter è un contatore a 30 bit che conta tick di 5 ns.
// Ogni volta che compare una parola di reset (bit31 = 1) il contatore
// si azzera. Per ottenere il tempo assoluto bisogna sommare,
// per ogni evento, un offset pari a (#reset visti)*2^30*tick.
//
// In questo macro useremo come unità di tempo i microsecondi.
//   1 tick        = 5 ns  = 0.005 µs
//   reset_t_us    = 2^30 * 0.005 µs ≈ 5.37·10^6 µs
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
const double tick_us    = 0.0048892;                             // 5 ns
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
    std::size_t index;      // indice della riga nel file originale
    double      t_us;       // tempo assoluto [µs]
    unsigned int ch;        // channel word "piena"
    bool isStart;
    bool isStop;
    unsigned int stopMask;  // ch & STOP_GENERIC_MASK

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
// Configura l'asse x di un istogramma tempo: range [tmin, tmax] e tick
// principali distanziati di tickStep (default 0.2 µs). Usa Ndivisions
// negativa per forzare esattamente quel passo senza suddivisioni secondarie.
// ---------------------------------------------------------------------
void ConfigureTimeAxis(TH1* h, double tmin, double tmax, double tickStep = 0.8 , double labelSize = 0.00006){
    if (!h) return;

    TAxis* ax = h->GetXaxis();
    ax->SetRangeUser(tmin, tmax);
    
    h->GetXaxis()->SetLabelSize(0.00008);  // testo più piccolo per evitare sovrapposizioni
    int ndiv = std::lround((tmax - tmin) / tickStep);
    if (ndiv < 1) ndiv = 1;

    //ndiv negativo -> nessuna ottimizzazione, nessuna sotto-divisione
    ax->SetNdivisions(-ndiv, /*optim=*/true);

}

// ---------------------------------------------------------------------
// Devince residuals (bin-by-bin) CONSISTENTI con un fit Poisson ("L").
// ---------------------------------------------------------------------
//
// Nel fit con opzione "L" ROOT minimizza la -log L assumendo che ogni
// bin segua una Poisson con media μ_i (attesa dal modello nel bin i).
// In questo contesto, i "pull" gaussiani (N-μ)/σ sono una diagnostica
// utile ma non perfettamente consistente quando N è piccolo.
//
// Una diagnostica più coerente con la Poisson likelihood è il
// *deviance residual* per bin:
//
//    d_i = 2 [ N_i ln(N_i/μ_i) - (N_i - μ_i) ]     (per N_i > 0)
//    d_i = 2 μ_i                                   (per N_i = 0, limite)
//
//    r_i = sign(N_i - μ_i) * sqrt(d_i)
//
// r_i è la radice della contribuzione del bin alla *Poisson deviance*
// (likelihood-ratio) tra il modello "saturated" (μ_i = N_i) e il
// modello fittato (μ_i = μ_i(θ)). In regime asintotico, la somma
// D = Σ d_i si distribuisce ~ χ² e i residui r_i tendono a N(0,1)
// se il modello è corretto.
//
// Coerenza con opzione "I":
// - Con "I" ROOT confronta i bin con la media del modello nel bin,
//   cioè (1/Δx) ∫_bin f(x) dx.
// - Se f(x) è parametrizzata in "counts/bin", allora
//     μ_i = (1/Δx) ∫_bin f(x) dx
//   è l'atteso in counts/bin, direttamente confrontabile con N_i.
//
// Nota tecnica ROOT: in molte versioni TF1::Integral NON è const,
// quindi qui usiamo TF1* (non const TF1*).
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

        // Considera solo i bin nel range di fit
        if (x1 <= fitMin || x0 >= fitMax) continue;

        const double Ni = h->GetBinContent(i);
        const double binw = (x1 - x0);

        // μ_i atteso nel bin (counts/bin) coerente con opzione "I"
        double mu = f->Integral(x0, x1) / binw;

        // Guard-rails numerici: μ deve essere > 0 per log(N/μ)
        if (mu <= 0.0) mu = 1e-12;

        // Contributo di devianza del bin: d_i >= 0
        double di = 0.0;
        if (Ni > 0.0) {
            di = 2.0 * (Ni * std::log(Ni / mu) - (Ni - mu));
        } else {
            // Limite Ni -> 0: Ni*log(Ni/mu) -> 0 e di -> 2*mu
            di = 2.0 * mu;
        }
        if (di < 0.0) di = 0.0; // robustezza numerica

        // Residuo di devianza con segno
        const double sgn = (Ni >= mu) ? 1.0 : -1.0;
        const double ri  = sgn * std::sqrt(di);

        hr->SetBinContent(i, ri);

        // Per un grafico "residui vs x" si usa spesso errore unitario.
        // (Non è un errore di misura: è un residuo normalizzato.)
        hr->SetBinError(i, 1.0);
    }

    return hr;
}

// ---------------------------------------------------------------------
// Metriche di qualità dei deviance residuals: mean, RMS, outlier, somma r^2.
// Per un modello corretto ci si aspetta indicativamente:
//   <r_D> ~ 0, RMS ~ 1, pochi outlier (|r_D|>3).
// La somma r^2 nel range di fit è legata alla devianza Poisson.
// ---------------------------------------------------------------------
struct ResidQuality {
    int    nBinsUsed   = 0;
    double mean        = 0.0;
    double rms         = 0.0;
    int    nAbsGt3     = 0;
    int    nAbsGt5     = 0;
    double fracAbsGt3  = 0.0;
    double fracAbsGt5  = 0.0;
    double sumR2       = 0.0; // somma r^2 (≈ devianza nel range)
    int    ndf         = 0;
    double r2OverNdf   = 0.0;
};

ResidQuality ComputeResidualQuality(const TH1* hres, double fitMin, double fitMax, int nPar=0)
{
    ResidQuality q;

    const int nb = hres->GetNbinsX();

    // Pass 1: media
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

    // Pass 2: RMS, outlier e somma r^2
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
    std::cout << "Sum(r_D^2): " << q.sumR2 << "   (≈ Poisson deviance nel range)\n";
    std::cout << "ndf:        " << q.ndf << "   (bins used - nPar)\n";
    std::cout << "r^2/ndf:    " << q.r2OverNdf << "   (diagnostico; ~1 se OK)\n";
    std::cout << "===================================================================\n";
}

// ---------------------------------------------------------------------
// Istogramma 1D della distribuzione dei deviance residuals (valori dei bin).
// Se il modello è adeguato ci si aspetta una distribuzione circa gaussiana
// centrata in 0 con σ~1 (soprattutto a statistiche sufficienti).
// ---------------------------------------------------------------------
TH1D* MakeResidualDistribution(const TH1* hres, double fitMin, double fitMax,
                               int nBins=60, double xmin=-6.0, double xmax=6.0,
                               const char* name="hDevResDist",
                               const char* title="Distribuzione residui (r_{D});r_{D};counts")
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
// ---------------------------------------------------------------------
// Helper: costruisce e fitta (LIRS) un modello exp+costante su un istogramma.
// Restituisce una TF1* (con nome univoco) già fitatta.
// ---------------------------------------------------------------------
TF1* FitExpPlusConst(TH1* h, double fitMin, double fitMax,
                     const char* baseName="fExpBkg",
                     bool silent=true)
{
    TString fname = Form("%s_%s", baseName, h->GetName());

    TF1* f = new TF1(fname.Data(), "[0]*exp(-x/[1]) + [2]", fitMin, fitMax);
    f->SetParNames("N0", "tau", "B");

    // Stime iniziali robuste
    f->SetParameter(0, h->GetMaximum());
    f->SetParameter(1, 2.2);

    // Stima del fondo dai bin di coda (ultimi 10 bin)
    double bkgGuess = 0.0;
    const int nb = h->GetNbinsX();
    const int nTail = std::min(10, nb);
    for (int ib = nb - nTail + 1; ib <= nb; ++ib) {
        bkgGuess += h->GetBinContent(ib);
    }
    bkgGuess /= (double)nTail;
    f->SetParameter(2, bkgGuess);

    // Fit: L (Poisson), I (integrale per bin), R (usa range TF1), S (FitResult), 0 (no draw)
    TString opt = silent ? "0LIRS" : "LIRS";
    h->Fit(f, opt.Data());

    return f;
}

// =====================================================================
//                          MU_LIFE_NEW
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

    // ------------------------------------------------------------
    // 1) Lettura file grezzo
    // ------------------------------------------------------------
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire il file " << filename << "\n";
        return;
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

    if (CH.empty() || CH.size() != CT.size()) {
        std::cerr << "[ERRORE] File vuoto o colonne di lunghezza diversa.\n";
        return;
    }

    std::cout << "[INFO] Righe lette: " << CH.size() << "\n";

    // ------------------------------------------------------------
    // 2) Costruzione vettore di Event con tempo assoluto
    // ------------------------------------------------------------
    std::vector<Event> events;
    events.reserve(CH.size());

    long long n_reset = -1;     // parte da -1, così il primo reset → offset 0
    bool seenFirstReset = false;

    for (std::size_t i = 0; i < CH.size(); ++i) {
        unsigned int ch = CH[i];
        unsigned int ct = CT[i];

        if (IsResetWord(ch)) {
            seenFirstReset = true;
            n_reset += 1;
            continue;
        }

        if (!seenFirstReset) {
            // Eventi di buffer prima del primo reset: li ignoriamo
            continue;
        }

        unsigned int ctr = (ct & COUNTER_MASK);
        double t_us_abs = (double)ctr * tick_us + (double)n_reset * reset_t_us;

        // Consideriamo solo eventi con almeno un bit significativo
        unsigned int mask = ch & (BIT_START | STOP_GENERIC_MASK);
        if (mask == 0u) continue;

        events.emplace_back(i, t_us_abs, ch);
    }

    std::cout << "[INFO] Eventi dopo il primo reset: " << events.size() << "\n";
    if (events.empty()) {
        std::cerr << "[ERRORE] Nessun evento utile dopo il primo reset.\n";
        return;
    }

    // ------------------------------------------------------------
    // 3) Loop principale: pairing START → STOP (senza goto)
    // ------------------------------------------------------------
    std::vector<double> dt_values;
    std::vector<unsigned int> stopBlocks;

    dt_values.reserve(20000);
    stopBlocks.reserve(20000);

    const std::size_t N = events.size();
    std::size_t i = 0;

    while (i < N) {
        const Event& evStart = events[i];

        if (!evStart.isStart) {
            ++i;
            continue;
        }

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

            // Se nel mezzo appare un nuovo START → scartiamo quello vecchio
            if (ev2.isStart) {
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
            ++i;
            continue;
        }

        // 2) Cerca lo STOP FINALE (BIT_STOP) entro FINAL_STOP_MAX_US
        bool foundFinalStop = false;
        std::size_t idxFinalStop = 0;

        for (std::size_t j = idxEarlyStop + 1; j < N; ++j) {
            const Event& ev2 = events[j];

            const double dt_tmp = ev2.t_us - tStart;
            if (dt_tmp > FINAL_STOP_MAX_US) break;

            if (ev2.isStart) {
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
            ++i;
            continue;
        }

        // 3) Coppia valida: calcola dt e salva maschera blocchi intorno allo stop
        const Event& evStop = events[idxFinalStop];
        const double dt = evStop.t_us - tStart;

        if (dt >= tmin && dt <= tmax) {
            const unsigned int finalBlockMask =
                CollectBlockMask(events, (int)idxFinalStop, FINAL_BLOCK_WINDOW);

            dt_values.push_back(dt);
            stopBlocks.push_back(finalBlockMask);
        }

        i = idxFinalStop + 1;
    }

    std::cout << "[INFO] Coppie START–STOP accettate: " << dt_values.size() << "\n";
    if (dt_values.empty()) {
        std::cerr << "[ATTENZIONE] Nessun dt ricostruito: controllare logica o parametri.\n";
        return;
    }

    // ------------------------------------------------------------
    // 4) Statistiche combinazioni PMT di stop (blocchi 8–11)
    // ------------------------------------------------------------
    std::map<unsigned int, long long> comboCounts;
    for (std::size_t k = 0; k < stopBlocks.size(); ++k) {
        const unsigned int sb = stopBlocks[k];
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
    // 5) Istogrammi: totale + per PMT del blocco
    // ------------------------------------------------------------
    TH1F* hDecay = new TH1F("hDecay",
                           "Istogramma intervalli dt START#minusSTOP; dt [#mu s]; Counts",
                           nbins, tmin, tmax);

    TH1F* hDecay_B8  = new TH1F("hDecay_B8",
                               "Istogramma intervalli dt START#minusSTOP PMT 8; dt [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B9  = new TH1F("hDecay_B9",
                               "Istogramma intervalli dt START#minusSTOP PMT 9; dt [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B10 = new TH1F("hDecay_B10",
                               "Istogramma intervalli dt START#minusSTOP PMT 10; dt [#mu s]; Counts",
                               nbins, tmin, tmax);
    TH1F* hDecay_B11 = new TH1F("hDecay_B11",
                               "Istogramma intervalli dt START#minusSTOP PMT 11; dt [#mu s]; Counts",
                               nbins, tmin, tmax);

    // // Asse x: parte da tmin e tick ogni 0.2 µs (richiesta)
    // ConfigureTimeAxis(hDecay,     tmin, tmax);
    // ConfigureTimeAxis(hDecay_B8,  tmin, tmax);
    // ConfigureTimeAxis(hDecay_B9,  tmin, tmax);
    // ConfigureTimeAxis(hDecay_B10, tmin, tmax);
    // ConfigureTimeAxis(hDecay_B11, tmin, tmax);
    // hDecay->GetYaxis()->SetLabelSize(0.00008);  // testo meno piccolo per evitare sovrapposizioni
    // Abilita la gestione errori per bin (utile per residui)
    hDecay->Sumw2();
    hDecay_B8->Sumw2();
    hDecay_B9->Sumw2();
    hDecay_B10->Sumw2();
    hDecay_B11->Sumw2();

    for (std::size_t k = 0; k < dt_values.size(); ++k) {
        const double dt = dt_values[k];
        const unsigned int sb = stopBlocks[k];

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
    gStyle->SetFitFormat("6.5g");
    gStyle->SetStripDecimals(kTRUE);


    // ------------------------------------------------------------
    // 6) Fit istogramma totale (LIRS) + deviance residuals + metriche
    // ------------------------------------------------------------
    TF1* fExpBkg = FitExpPlusConst(hDecay, tmin, tmax, "fExpBkg", true);

    const double tau  = fExpBkg->GetParameter(1);
    const double etau = fExpBkg->GetParError(1);
    const double B    = fExpBkg->GetParameter(2);
    const double eB   = fExpBkg->GetParError(2);

    std::cout << "\n================ RISULTATI FIT (TOTALE) ================\n";
    std::cout << "Tau      = " << tau  << " ± " << etau << " us\n";
    std::cout << "B (fondo)= " << B    << " ± " << eB   << " counts/bin\n";
    std::cout << "========================================================\n";

    TH1D* hDevRes = MakeDevianceResidualHist(hDecay, fExpBkg, tmin, tmax,
                              "hDevRes",
                              "Deviance residuals; t_{decay} [#mu s]; r_{D}");

    // Metriche quantitative dei residui (più appropriate di un chi2/ndf quando si usa la likelihood Poisson)
    const ResidQuality qTot = ComputeResidualQuality(hDevRes, tmin, tmax, /*nPar=*/3);
    PrintResidualQuality(hDevRes, tmin, tmax, /*nPar=*/3, "(TOTAL)");

    // Istogramma 1D della distribuzione dei deviance residuals + fit gaussiano
    TH1D* hDevResDist = MakeResidualDistribution(hDevRes, tmin, tmax,
                                          60, -6, 6,
                                          "hDevResDist",
                                          "Distribuzione residui (r_{D});r_{D};counts");
    TF1* fGaus = new TF1("fDevResGaus", "gaus", -4, 4);
    hDevResDist->Fit(fGaus, "0RQ"); // 0=no draw, R=range, Q=quiet

    // ------------------------------------------------------------
    // 7) Canvas: sopra fit, sotto deviance residuals (totale)
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

    // estetica: togli label x sopra
    hDecay->GetXaxis()->SetLabelSize(0.0);
    hDecay->GetXaxis()->SetTitleSize(0.0);

    hDecay->Draw("E1");
    fExpBkg->Draw("same");
    
    // Scelgo altezza e gap del box sotto (in NDC)
    const double height2 = 0.22;  // altezza seconda legenda (2 righe circa)

    // Box (testo) con una misura di goodness-of-fit coerente con il fit "L":
    // - Sum(r_D^2) ~ deviance Poisson nel range
    // - (Sum(r_D^2))/ndf dovrebbe essere ~1 se il modello descrive i dati.
    TPaveText* pGoF = new TPaveText(0.40, 0.80  , 0.62, 0.9 - height2 , "NDC");
    pGoF->SetFillStyle(0);
    pGoF->SetBorderSize(1);
    pGoF->SetTextAlign(12);
    pGoF->SetTextSize(0.03);
    pGoF->AddText(Form("Deviance/ndf = %.2f", qTot.r2OverNdf));
    pGoF->AddText(Form("<r_{D}> = %.3f  RMS = %.3f", qTot.mean, qTot.rms));
    pGoF->Draw();

    // -----------------------------------------------------------------
    // Nota grafica richiesta:
    // - questa legenda (range + binning) deve stare *dentro* il grafico
    //   e *sotto* il box dei parametri del fit che ROOT disegna in alto
    //   a destra (gStyle->SetOptFit(1)).
    // - NON riportiamo qui #tau, perché è già nel box del fit.
    //
    // Coordinate in NDC (0..1). Se il box di ROOT dovesse cambiare
    // posizione/dimensione, puoi ritoccare questi numeri.
    // -----------------------------------------------------------------
    TLegend *leg = new TLegend(0.40, 0.80, 0.62, 0.9);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    leg->AddEntry((TObject*)0, Form("[%.1f, %.1f] #mus", tmin, tmax), "");
    leg->AddEntry((TObject*)0, Form("Bins: %d", nbins), "");
    leg->Draw();

    /*
    // Box diagnostico: usa le metriche derivate dai deviance residuals.
    // Questo è coerente con il fit in likelihood (opzione "L").
    TPaveText* pDiag = new TPaveText(0.70, 0.61, 0.93, 0.69, "NDC");
    pDiag->SetFillStyle(0);
    pDiag->SetBorderSize(1);
    pDiag->SetTextAlign(12);
    pDiag->SetTextSize(0.028);
    pDiag->AddText(Form("Deviance/ndf = %.3f", qTot.r2OverNdf));
    pDiag->AddText(Form("<r_D> = %.3f,  RMS = %.3f", qTot.mean, qTot.rms));
    pDiag->Draw();
    */

    pad2->cd();
    hDevRes->SetStats(false);
    hDevRes->SetMinimum(-5);
    hDevRes->SetMaximum(5);

    hDevRes->GetYaxis()->SetNdivisions(505);
    hDevRes->GetYaxis()->SetTitleSize(0.10);
    hDevRes->GetYaxis()->SetLabelSize(0.08);
    hDevRes->GetYaxis()->SetTitleOffset(0.45);

    hDevRes->GetXaxis()->SetTitleSize(0.12);
    hDevRes->GetXaxis()->SetLabelSize(0.08);

    hDevRes->Draw("E1");

    // linee guida
    TLine* l0  = new TLine(tmin,  0.0, tmax,  0.0);
    TLine* lp3 = new TLine(tmin,  3.0, tmax,  3.0);
    TLine* lm3 = new TLine(tmin, -3.0, tmax, -3.0);
    l0->SetLineStyle(2);
    lp3->SetLineStyle(3);
    lm3->SetLineStyle(3);
    l0->Draw("same");
    lp3->Draw("same");
    lm3->Draw("same");

    // Salvataggio canvas fit + deviance residuals
    TString outName;
    outName.Form("mu_decay_fit_%.3f_%.3f_%d_bins_withdevres.png", tmin, tmax, nbins);
    c1->SaveAs(outName.Data());

    // Salvataggio distribuzione dei deviance residuals
    TCanvas* cResDist = new TCanvas("cResDist", "Distribuzione residui (r_{D})", 700, 600);
    hDevResDist->Draw("E1");
    fGaus->Draw("same");
    TString outResName;
    outResName.Form("devres_distribution_%.3f_%.3f_%d_bins.png", tmin, tmax, nbins);
    cResDist->SaveAs(outResName.Data());

    // ------------------------------------------------------------
    // 8) Subplot PMT: opzionali da terminale, CON fit + deviance residuals per ognuno
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

            TF1* ff = FitExpPlusConst(hh, tmin, tmax, "fExpBkg", true);

            TH1D* hp = MakeDevianceResidualHist(hh, ff, tmin, tmax,
                                    Form("hDevRes_%s", labels[ih]),
                                    Form("Deviance residuals (%s); t_{decay} [#mu s]; r_{D}", labels[ih]));

            // Metriche quantitative dei residui per questo istogramma
            const ResidQuality qH = ComputeResidualQuality(hp, tmin, tmax, /*nPar=*/3);
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

            // Pad alto
            p1->cd();
            hh->SetMarkerStyle(20);
            hh->SetMarkerSize(1.0);
            hh->SetMarkerColor(kBlack);
            hh->GetXaxis()->SetLabelSize(0.0);
            hh->GetXaxis()->SetTitleSize(0.0);
            hh->Draw("E1");
            ff->Draw("same");

            // Come richiesto: legenda informativa interna al grafico,
            // sotto il box del fit ROOT. Non ripetiamo #tau qui.
            TLegend* l = new TLegend(0.40, 0.80, 0.62, 0.9);
            l->SetBorderSize(1);
            l->SetFillStyle(0);
            l->SetTextSize(0.03);
            l->AddEntry((TObject*)0, Form("%s  [%.1f, %.1f] #mus", labels[ih], tmin, tmax), "");
            l->AddEntry((TObject*)0, Form("Bins: %d", nbins), "");
            l->Draw();

            // Box diagnostico coerente con likelihood (Poisson deviance / ndf)
            TPaveText* pDiagH = new TPaveText(0.40, 0.80  , 0.62, 0.9 - height2 , "NDC");
            pDiagH->SetFillStyle(0);
            pDiagH->SetBorderSize(1);
            pDiagH->SetTextAlign(12);
            pDiagH->SetTextSize(0.028);
            pDiagH->AddText(Form("Deviance/ndf = %.3f", qH.r2OverNdf));
            pDiagH->AddText(Form("<r_D> = %.3f,  RMS = %.3f", qH.mean, qH.rms));
            pDiagH->Draw();

            // Pad basso
            p2->cd();
            hp->SetStats(false);
            hp->SetMinimum(-5);
            hp->SetMaximum(5);

            hp->GetYaxis()->SetNdivisions(505);
            hp->GetYaxis()->SetTitleSize(0.10);
            hp->GetYaxis()->SetLabelSize(0.08);
            hp->GetYaxis()->SetTitleOffset(0.45);

            hp->GetXaxis()->SetTitleSize(0.12);
            hp->GetXaxis()->SetLabelSize(0.08);

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
        }
    }

    std::cout << "[INFO] Fine.\n";
}
