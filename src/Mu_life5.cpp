#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// ROOT
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TFile.h"
#include <TLegend.h>

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
const unsigned int BIT_STOP  = 1u << 1;    // 2
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
const int   EARLY_STOP_MAX_TICKS  = 10;     // stop "immediato" entro 10 eventi dopo lo start
const double FINAL_STOP_MAX_US    = 20.0;   // stop fisico entro 20 µs dallo start
const int   EARLY_BLOCK_WINDOW    = 2;      // ±2 eventi per stimare i blocchi dello stop "immediato"
const int   FINAL_BLOCK_WINDOW    = 3;      // ±3 eventi per stimare i blocchi dello stop finale

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

// =====================================================================
//                          MU_LIFE_NEW
// =====================================================================

void Mu_life_new(const char* filename = "FIFOread_Take5.txt",
                 int nbins = 80,
                 double tmin = 0.0,
                 double tmax = 20.0)
{
    std::cout << "\n============================================\n";
    std::cout << "[Mu_life_new] File: " << filename << "\n";
    std::cout << "[Mu_life_new] Finestra istogramma dt: ["
              << tmin << ", " << tmax << "] µs\n";
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

        bool isReset = IsResetWord(ch);
        if (isReset) {
            seenFirstReset = true;
            n_reset += 1;
            continue;
        }

        if (!seenFirstReset) {
            // Eventi di buffer prima del primo reset: li ignoriamo
            continue;
        }

        unsigned int ctr = (ct & COUNTER_MASK);
        double t_us = (double)ctr * tick_us + (double)n_reset * reset_t_us;

        // Consideriamo solo eventi con almeno un bit significativo
        unsigned int mask = ch & (BIT_START | STOP_GENERIC_MASK);
        if (mask == 0u) continue;

        events.emplace_back(i, t_us, ch);
    }

    std::cout << "[INFO] Eventi dopo il primo reset: " << events.size() << "\n";
    if (events.empty()) {
        std::cerr << "[ERRORE] Nessun evento utile dopo il primo reset.\n";
        return;
    }

    // ------------------------------------------------------------
    // 3) Loop principale: pairing START → STOP (senza goto)
    // ------------------------------------------------------------
    std::vector<double> dt_values;          // tempi di decadimento
    std::vector<unsigned int> startBlocks;  // blocchi allo stop "immediato"
    std::vector<unsigned int> stopBlocks;   // blocchi allo stop finale

    dt_values.reserve(10000);
    startBlocks.reserve(10000);
    stopBlocks.reserve(10000);

    std::size_t N = events.size();
    std::size_t i = 0;      // indice principale sugli eventi

    while (i < N) {
        const Event& evStart = events[i];

        if (!evStart.isStart) {
            ++i;
            continue;
        }

        std::size_t idxStart = i;
        double tStart = evStart.t_us;

        // Flag per sapere se scartiamo lo start a metà procedura
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
                // ripartiamo a elaborare da questo nuovo start
                i = j;
                discardThisStart = true;
                break;
            }

            // stop generico: almeno un bit tra 2,4,8,16,32
            if ((ev2.stopMask & STOP_GENERIC_MASK) != 0u) {
                foundEarlyStop = true;
                idxEarlyStop = j;
                break;
            }
        }

        if (discardThisStart) {
            // passiamo direttamente al prossimo ciclo del while (nuovo start)
            continue;
        }

        if (!foundEarlyStop) {
            // nessun stop immediato → scarta questo start e vai avanti
            ++i;
            continue;
        }

        // blocchi attivi intorno allo stop immediato
        unsigned int earlyBlockMask =
            CollectBlockMask(events, (int)idxEarlyStop, EARLY_BLOCK_WINDOW);

        // 2) Cerca lo STOP FINALE (bit2 acceso) entro 20 µs dallo start
        bool foundFinalStop = false;
        std::size_t idxFinalStop = 0;

        for (std::size_t j = idxEarlyStop + 1; j < N; ++j) {
            const Event& ev2 = events[j];

            double dt_tmp = ev2.t_us - tStart;
            if (dt_tmp > FINAL_STOP_MAX_US) {
                // oltre la finestra di 20 µs → stop non trovato
                break;
            }

            // se appare un nuovo START prima dello stop finale → scartiamo
            if (ev2.isStart) {
                i = j;  // ripartiamo da questo nuovo start
                discardThisStart = true;
                break;
            }

            // STOP finale: richiediamo il bit2 (STOP generale)
            if ((ev2.ch & BIT_STOP) != 0u) {
                foundFinalStop = true;
                idxFinalStop = j;
                break;
            }
        }

        if (discardThisStart) {
            // abbiamo trovato un nuovo start nel mezzo: questo vecchio è scartato
            continue;
        }

        if (!foundFinalStop) {
            // nessun stop finale → scarta questo start e vai avanti
            ++i;
            continue;
        }

        // 3) Abbiamo una coppia START–STOP: calcola dt e blocchi dello stop finale
        const Event& evStop = events[idxFinalStop];
        double dt = evStop.t_us - tStart;

        if (dt >= tmin && dt <= tmax) {
            unsigned int finalBlockMask =
                CollectBlockMask(events, (int)idxFinalStop, FINAL_BLOCK_WINDOW);

            dt_values.push_back(dt);
            startBlocks.push_back(earlyBlockMask);
            stopBlocks.push_back(finalBlockMask);
        }

        // Dopo aver gestito questa coppia, ripartiamo dall'evento successivo allo stop
        i = idxFinalStop + 1;
    }

    std::cout << "[INFO] Coppie START–STOP accettate: "
              << dt_values.size() << "\n";

    // ------------------------------------------------------------
    // Statistiche sulle combinazioni di PMT del blocco per gli stop
    // ------------------------------------------------------------
    std::map<unsigned int, long long> comboCounts;

    for (std::size_t k = 0; k < stopBlocks.size(); ++k) {
        unsigned int sb = stopBlocks[k];
        if (sb == 0u) continue;  // se per qualche motivo non c'è nessun PMT del blocco, ignora

        comboCounts[sb]++;
    }

    std::cout << "\n[INFO] Statistiche combinazioni PMT di stop (blocchi 8–11):\n";
    if (comboCounts.empty()) {
        std::cout << "  Nessun evento di stop con PMT del blocco attivo.\n";
    } else {
        for (const auto &kv : comboCounts) {
            unsigned int mask = kv.first;
            long long    n    = kv.second;

            std::cout << "  PMT blocco = { ";

            bool first = true;
            if (mask & BIT_B8)  { std::cout << (first ? "" : ", ") << "8";  first = false; }
            if (mask & BIT_B9)  { std::cout << (first ? "" : ", ") << "9";  first = false; }
            if (mask & BIT_B10) { std::cout << (first ? "" : ", ") << "10"; first = false; }
            if (mask & BIT_B11) { std::cout << (first ? "" : ", ") << "11"; first = false; }

            std::cout << " }  ->  " << n << " eventi";

            // (opzionale) stampo anche la maschera in esadecimale
            std::cout << "   [mask = 0x" << std::hex << mask << std::dec << "]\n";
        }
    }
    std::cout << std::endl;

    if (dt_values.empty()) {
        std::cerr << "[ATTENZIONE] Nessun dt ricostruito: controllare logica o parametri.\n";
        return;
    }

    // 4) Istogramma e fit esponenziale + fondo
    // ------------------------------------------------------------
    TH1F* hDecay = new TH1F("hDecay",
                            "Muon decay time; t_{decay} [#mu s]; Counts",
                            nbins, tmin, tmax);

    // Istogrammi separati per i diversi PMT del blocco
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

    // Riempiamo l'istogramma totale e quelli per PMT in base alla maschera dei blocchi
    for (std::size_t k = 0; k < dt_values.size(); ++k) {
        double dt = dt_values[k];
        unsigned int sb = (k < stopBlocks.size()) ? stopBlocks[k] : 0u;

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

    // Modello: N(t) = N0 * exp(-t/tau) + B
    TF1* fExpBkg = new TF1("fExpBkg",
                           "[0]*exp(-x/[1]) +[2]",
                           tmin, tmax);
    fExpBkg->SetParNames("N0", "tau", "B");

    // Stime iniziali
    fExpBkg->SetParameter(0, hDecay->GetMaximum());
    fExpBkg->SetParameter(1, 2.2);   // µs, tempo di vita atteso

    // (segue il resto del codice del fit come prima)


    //Stima grezza del fondo costante dai bin di coda
    double bkgGuess = 0.0;
    int nb = hDecay->GetNbinsX();
    int nTail = std::min(10, nb);
    for (int ib = nb - nTail + 1; ib <= nb; ++ib) {
      bkgGuess += hDecay->GetBinContent(ib);
    }
    bkgGuess /= (double)nTail;
    fExpBkg->SetParameter(2, bkgGuess);

    hDecay->Fit(fExpBkg, "LIR+");

    double tau  = fExpBkg->GetParameter(1);
    double etau = fExpBkg->GetParError(1);
    double B    = fExpBkg->GetParameter(2);
    double eB   = fExpBkg->GetParError(2);

    std::cout << "\n================ RISULTATI FIT ================\n";
    std::cout << "Tau (µ)  = " << tau  << " ± " << etau << " µs\n";
    std::cout << "B (fondo)= " << B    << " ± " << eB   << " counts/bin\n";
    std::cout << "==============================================\n";

    TCanvas* c1 = new TCanvas("c1", "Muon lifetime", 1000, 700);
    hDecay->Draw();
    fExpBkg->Draw("same");
    //c1->SaveAs("Mu_life_new_1.png");
    TLegend *leg = new TLegend(0.3, 0.80, 0.5, 0.90);
    leg->SetBorderSize(1);
    leg->SetFillStyle(0);
    leg->AddEntry((TObject*)0,
                  Form("[%.3f, %.3f] #mus", tmin, tmax),
                  "");
    leg->AddEntry((TObject*)0,
                  Form("Bins: %d", nbins),
                  "");
    leg->Draw();

    TString outName;
    outName.Form("mu_decay_fit_%.3f_%.3f_%d_bins.png", tmin, tmax, nbins);
    //outName.ReplaceAll(".", "p"); // evita i punti nel nome file
    c1->SaveAs(outName.Data());

    // fit con doppio esponenziale per trovare rapporto popolazioni
    
    float tau_meno = 2.043; //vita media mu-
    float tau_plus;
    tau_plus = 2*tau - tau_meno;

    TF1* fExpSum = new TF1("fExpSum",
                           Form("[0]*exp(-x/%f) + [1]*exp(-x/%f) + [2]", tau_meno, tau_plus),
                           tmin, tmax);
    fExpSum->SetParNames("N-", "N+", "C");

    // Stime iniziali
    fExpSum->SetParameter(0, hDecay->GetMaximum()*0.5);
    fExpSum->SetParameter(1, hDecay->GetMaximum()*0.5);
    fExpSum->SetParameter(2, bkgGuess);

    hDecay->Fit(fExpSum, "LIR+");
    double N_minus    = fExpSum->GetParameter(0);
    double eN_minus   = fExpSum->GetParError(0);
    double N_plus  = fExpSum->GetParameter(1);
    double eN_plus = fExpSum->GetParError(1);
    double B_1=      fExpSum->GetParameter(2);
    double errorB_1 = fExpSum->GetParError(2); 

    double abb;
    double error_abb;
    abb = N_plus / N_minus;
    error_abb = (N_plus / N_minus) * sqrt( pow((eN_plus / N_plus),2) + pow((eN_minus / N_minus ),2) );

    std::cout << "\n================ RISULTATI FIT ================\n";
    std::cout << "N_minus  = " << N_minus  << " ± " << eN_minus << " \n";
    std::cout << "N_minus  = " << N_plus  << " ± " << eN_plus << " \n";
    std::cout << "B_1 (fondo)= " << B_1    << " ± " << errorB_1   << " counts/bin\n";
    std::cout << "abbondanze  = " << abb  << " ± " << error_abb << " \n";
    std::cout << "==============================================\n";




    //abilito o meno i subplot dei vari PMT8,9,10,11
    char choice_subplot;
    std::cout << "Vuoi vedere i subplot PMT8,9,10,11?S/N   ";
    std:: cin >> choice_subplot;

    if(choice_subplot == 'S'){

    // Canvas per i singoli PMT (solo istogrammi, senza fit)
    TCanvas* c2 = new TCanvas("c2", "Muon lifetime - stop PMT 8", 800, 600);
    hDecay_B8->Draw();

    TCanvas* c3 = new TCanvas("c3", "Muon lifetime - stop PMT 9", 800, 600);
    hDecay_B9->Draw();

    TCanvas* c4 = new TCanvas("c4", "Muon lifetime - stop PMT 10", 800, 600);
    hDecay_B10->Draw();

    TCanvas* c5 = new TCanvas("c5", "Muon lifetime - stop PMT 11", 800, 600);
    hDecay_B11->Draw();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
}

    TFile* fout = new TFile("Mu_life_new.root", "RECREATE");
    hDecay->Write();
    hDecay_B8->Write();
    hDecay_B9->Write();
    hDecay_B10->Write();
    hDecay_B11->Write();
    fExpBkg->Write();
    c1->Write();
    fout->Close();

    std::cout << "[INFO] Risultati salvati in Mu_life_new.root\n";
}
