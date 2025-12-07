//codice realizzato per l'analisi di waveforms provenienti da un sistema di acquisizione a 4 canali dell'adc converter caen N6725
//il codice è da intendersi come una macro di ROOT e va eseguito all'interno di una shell ROOT
//per eseguire la macro da shell ROOT digitare:
//.L Calib_Cosmic_ray.cpp
// in seguito per eseguire la macro digitare:
//WaveAnalysis("wave0.txt", "wave1.txt", "wave2.txt", "wave3.txt", minStepAmp, binWidth);
//dove wave0.txt, wave1.txt, wave2.txt, wave3.txt sono i file di waveforms acquisiti
//minStepAmp è la soglia minima di ampiezza del gradino per essere considerato valido (default 200 conteggi adc della somma dei 4 canali)
//binWidth è la larghezza dei bin dell'istogramma delle ampiezze dei gradini (default 100 conteggi adc)
//la macro legge i 4 file di waveforms, combina gli eventi con lo stesso event number
//e per ogni evento fa il fit della somma dei 4 canali con una sigmoide per estrarre l'ampiezza del gradino
//infine crea un istogramma con la distribuzione delle ampiezze dei gradini







#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <TSystem.h>  // per gSystem->ProcessEvents()
#include <chrono>

struct EventWave {
    int eventNumber;                 // "Event Number: X"
    std::vector<unsigned short> ch;  // ch[i] = valore ADC al campione i
};

// Evento con 4 canali e waveform sommata
struct MultiChannelEvent {
    int eventNumber = -1;
    std::vector<unsigned short> ch[4];      // ch[0], ch[1], ch[2], ch[3]
    std::vector<unsigned int>   sum;        // somma campione per campione
};

// ======================
//   LETTURA FILE WAVEFORMS
// ======================
// Legge un file di waveforms e li ritorna in una mappa
//  key =   numero evento
//  value = EventWave (struttura con eventNumber e vettore di campioni  ADC)



std::map<int, EventWave> ReadWaveFile(const std::string &fname)     // legge un file di waveforms e li ritorna in una mappa 
{                                                                  // key = numero evento, value = EventWave




    std::ifstream in(fname);            // apre il file in lettura
    if (!in.is_open()) {            // controlla che il file sia stato aperto correttamente
        std::cerr << "Impossibile aprire " << fname << std::endl;
        return {};
    }




std::map<int, EventWave> events;         // mappa di eventi

    int currentEvent = -1;           // numero dell'evento corrente
    EventWave *currentWave = nullptr;  // puntatore all'EventWave corrente

    std::string line;                 // variabile per leggere le righe del file

    while (std::getline(in, line)) {
        if (line.empty())     {   // salta righe vuote
            continue;
        }

    // 🔹 1) Caso: riga con "Event Number"
    if (line.find("Event Number") != std::string::npos) { // trova la riga con "Event Number"
        size_t pos = line.find(':');       // trova i due punti dopo "Event Number"
        if (pos != std::string::npos) {     // se trova i due punti trova la posizione
            int ev;
            std::stringstream ss(line.substr(pos+1));   // estrae la sottostringa dopo i due punti e la mette in uno stringstream
            ss >> ev;    // legge il numero dell'evento
        currentEvent = ev;   // aggiorna il numero dell'evento corrente
        EventWave &ew = events[currentEvent];    // crea un nuovo EventWave nella mappa per questo evento
        ew.eventNumber = currentEvent;    // imposta il numero dell'evento
        ew.ch.clear();   // pulisce il vettore dei campioni ADC
        currentWave = &ew;  // aggiorna il puntatore all'EventWave corrente
    }
    continue;  // vai alla prossima riga del file
}

    // 🔹 2) Caso: riga con campione ADC (numero intero)
    if (std::isdigit(line[0]) || line[0] == '-') {
        if (!currentWave) {
        // se NON abbiamo ancora visto un Event Number, questi dati non sappiamo dove metterli
            continue;
    }


        unsigned int adc; // valore ADC
        std::stringstream ss(line); // string stream leggo la linea del file
        ss >> adc;  // estraggo il valore ADC
        currentWave->ch.push_back(static_cast<unsigned short>(adc)); // aggiungo il valore ADC al vettore dei campioni
        continue; // vai alla prossima riga del file
        }
        
        // qualsiasi altra cosa (spazi strani ecc.) viene ignorata
    }

    return events;      // ritorna la mappa di eventi
}


// ======================
//   ANALISI DEI GRADINI
// ======================
//
// Analizza tutti gli eventi della mappa "merged":
//  - per ogni evento usa la waveform "sum" (somma dei 4 canali)
//  - applica un cut su |ampiezza stimata| >= minStepAmp
//  - fa il fit con una sigmoide
//       f(x) = y0 + amp / (1 + exp(-(x - x0)/width))
//    permettendo anche amp < 0 (gradini verso il basso)
//  - riempie un istogramma con |amp| (altezza del gradino)
//
//  N.B. minStepAmp è in conteggi ADC della SOMMA (4 canali).
//
void AnalyzeSteps(const std::map<int, MultiChannelEvent> &merged,
                  double minStepAmp = 200.0,// soglia su |ampiezza|
                  double binWidth = 100.0)   // larghezza bin istogramma
{
    if (merged.empty()) {
        std::cout << "[AnalyzeSteps] La mappa 'merged' è vuota.\n";
        return;
    }

    // 1) vettore per salvare |amp| di tutti gli eventi buoni
    std::vector<double> ampValues;
    ampValues.reserve(merged.size()); // per evitare riallocazioni continue

    int nFittedEvents = 0;

    // 2) Vettori x,y riutilizzati ad ogni evento
    std::vector<double> x, y;

    // 3) TGraph riutilizzato ad ogni evento (senza new/delete nel loop)
    TGraph g;  // grafico "stack" che vive per tutta la funzione

    // 4) TF1 sigmoide riutilizzata
    //    Range iniziale  [0, 1] verrà adattato ad ogni evento
    TF1 fsig("fsig", "[0] + [1]/(1 + exp(-(x-[2])/[3]))", 0.0, 1.0);
    fsig.SetParName(0, "y0");
    fsig.SetParName(1, "amp");
    fsig.SetParName(2, "x0");
    fsig.SetParName(3, "width");

    // Loop sugli eventi
    for (const auto &kv : merged) {
        int ev = kv.first;
        const MultiChannelEvent &m = kv.second;

        const auto &wave = m.sum;
        size_t N = wave.size();
        if (N < 10) continue;

        // ---------------- Baseline e ampiezza stimata ----------------
        double minVal = wave[0];
        double maxVal = wave[0];
        for (size_t i = 1; i < N; ++i) {
            if (wave[i] < minVal) minVal = wave[i];
            if (wave[i] > maxVal) maxVal = wave[i];
        }

        size_t nBaseline = std::min<size_t>(50, N/4);
        double baseline = 0.0;
        for (size_t i = 0; i < nBaseline; ++i) baseline += wave[i];
        baseline /= nBaseline;
        double y0_init = baseline;

        double upAmp   = maxVal - y0_init;
        double downAmp = minVal - y0_init;
        double amp_init = (std::fabs(upAmp) >= std::fabs(downAmp)) ? upAmp : downAmp;
        double ampAbs   = std::fabs(amp_init);
        if (ampAbs < minStepAmp) continue;  // gradino troppo piccolo

        // ---------------- Stima x0, width (come prima) ---------------
        double midLevel = y0_init + 0.5 * amp_init;
        double y10      = y0_init + 0.1 * amp_init;
        double y90      = y0_init + 0.9 * amp_init;

        size_t iMid = 0;
        if (amp_init > 0) {
            for (size_t i = 0; i < N; ++i)
                if (wave[i] >= midLevel) { iMid = i; break; }
        } else {
            for (size_t i = 0; i < N; ++i)
                if (wave[i] <= midLevel) { iMid = i; break; }
        }
        double x0_init = static_cast<double>(iMid);

        size_t i10 = 0, i90 = N - 1;
        if (amp_init > 0) {
            for (size_t i = 0; i < N; ++i)
                if (wave[i] >= y10) { i10 = i; break; }
            for (size_t i = i10; i < N; ++i)
                if (wave[i] >= y90) { i90 = i; break; }
        } else {
            for (size_t i = 0; i < N; ++i)
                if (wave[i] <= y10) { i10 = i; break; }
            for (size_t i = i10; i < N; ++i)
                if (wave[i] <= y90) { i90 = i; break; }
        }

        double width_init = 5.0;
        if (i90 > i10) {
            width_init = (static_cast<double>(i90) - static_cast<double>(i10)) / 4.4;
            if (width_init <= 0.0) width_init = 5.0;
        }

        // ---------------- Riusa vettori x,y e TGraph -----------------
        x.resize(N);
        y.resize(N);
        for (size_t i = 0; i < N; ++i) {
            x[i] = static_cast<double>(i);
            y[i] = static_cast<double>(wave[i]);
        }

        g.Set(N); // imposta il numero di punti
        for (int i = 0; i < (int)N; ++i) {
            g.SetPoint(i, x[i], y[i]);
        }
        g.SetTitle(Form("Evento %d;indice campione;somma ADC", ev));

        // ---------------- Settaggio TF1 e fit ------------------------
        fsig.SetRange(0.0, static_cast<double>(N - 1));
        fsig.SetParameter(0, y0_init);
        fsig.SetParameter(1, amp_init);
        fsig.SetParameter(2, x0_init);
        fsig.SetParameter(3, width_init);
        fsig.SetParLimits(1, -2.0 * ampAbs, 2.0 * ampAbs);
        fsig.SetParLimits(3, 0.1, static_cast<double>(N));

        int fitStatus = g.Fit(&fsig, "QN");  // niente draw, niente stampe
        if (fitStatus != 0) continue;

        double ampFit    = fsig.GetParameter(1);
        double ampFitAbs = std::fabs(ampFit);

        ampValues.push_back(ampFitAbs);
        ++nFittedEvents;
    }

    std::cout << "[AnalyzeSteps] Eventi fittati con successo: "
              << nFittedEvents << std::endl;

    if (ampValues.empty()) {
        std::cout << "[AnalyzeSteps] Nessuna ampiezza valida, non creo istogramma.\n";
        return;
    }

    // ---------- Costruzione istogramma con range automatico ----------
    auto [minIt, maxIt] = std::minmax_element(ampValues.begin(), ampValues.end());
    double minAmp = *minIt;
    double maxAmp = *maxIt;
    double margin = 1000.0;

    double histMin = minAmp;          // se vuoi puoi mettere 0.0
    double histMax = maxAmp + margin;

    
    
    int nbins = std::max(1, int((histMax - histMin) / binWidth));

    TH1D *hAmp = new TH1D("hAmp",
                          "Distribuzione |amp| (altezza gradino);|amp| [ADC];entries",
                          nbins, histMin, histMax);

    for (double a : ampValues) {
        hAmp->Fill(a);
    }

    TCanvas *cAmp = new TCanvas("cAmp", "Distribuzione |amp|", 800, 600);
    hAmp->Draw();

    double mean  = hAmp->GetMean();
    double sigma = hAmp->GetRMS();

    std::cout << "Distribuzione |amp|:\n"
              << "  Min  = " << minAmp  << "\n"
              << "  Max  = " << maxAmp  << "\n"
              << "  Mean = " << mean    << "\n"
              << "  RMS  = " << sigma   << std::endl;
}



// ======================================================================
//  Browser interattivo degli eventi con UN SOLO CANVAS
//  - Mostra la somma dei 4 canali per un evento alla volta
//  - Usa sempre lo stesso TCanvas, ridisegnando la waveform
//  - Dopo ogni evento chiede se continuare o fermarsi
//
//  Argomenti:
//    merged : mappa <eventNumber, MultiChannelEvent> già riempita
//    dt     : scala sull'asse x (1.0 = indice; 4.0 = 4 ns, ecc.)
// ======================================================================
void BrowseSummedEventsSingleCanvas(const std::map<int, MultiChannelEvent> &merged,
                                    double dt = 1.0)
{
    if (merged.empty()) {
        std::cout << "[BrowseSummedEvents] La mappa 'merged' è vuota.\n";
        return;
    }

    std::cout << "[BrowseSummedEvents] Numero di eventi disponibili: "
              << merged.size() << std::endl;

    // Creo UN SOLO canvas, che userò per tutti gli eventi
    TCanvas *c = new TCanvas("c_sum_browser",
                             "Somma 4 canali - browser eventi",
                             800, 600);

    // Puntatore al grafico corrente (così posso cancellarlo prima di disegnare il successivo)
    TGraph *g_sum = nullptr;

    char choice = 'y';

    // Ciclo sugli eventi in ordine di eventNumber
    for (auto it = merged.begin();
         it != merged.end() && (choice == 'y' || choice == 'Y');
         ++it)
    {
        int evToPlot = it->first;
        const MultiChannelEvent &m = it->second;

        size_t Nsum = m.sum.size();
        if (Nsum == 0) {
            std::cerr << "[BrowseSummedEvents] Evento " << evToPlot
                      << " ha sum vuota, salto.\n";
            continue;
        }

        // Cancello il grafico precedente, se esiste
        if (g_sum) {
            delete g_sum;
            g_sum = nullptr;
        }

        // Costruisco i vettori x,y per questo evento
        std::vector<double> x_sum(Nsum), y_sum(Nsum);
        for (size_t i = 0; i < Nsum; ++i) {
            x_sum[i] = static_cast<double>(i) * dt;      // indice (o tempo)
            y_sum[i] = static_cast<double>(m.sum[i]);    // somma ADC
        }

        // Creo un nuovo TGraph per questo evento (solo in memoria, non nuovo canvas)
        g_sum = new TGraph(Nsum, x_sum.data(), y_sum.data());
        g_sum->SetTitle(Form("Evento %d, somma canali;indice campione;ADC sommati",
                             evToPlot));
        g_sum->SetLineWidth(2);

        // Disegno sul canvas esistente
        c->cd();
        g_sum->Draw("AL");
        c->SetTitle(Form("Somma 4 canali - evento %d", evToPlot));
        c->Modified();
        c->Update();
        gSystem->ProcessEvents();

        // Interazione con l'utente
        std::cout << "\n[BrowseSummedEvents] Mostrato evento " << evToPlot << ".\n"
                  << "Vuoi vedere l'evento successivo? (y/n): ";
        std::cin >> choice;
    }

    std::cout << "[BrowseSummedEvents] Fine esplorazione eventi.\n";
}





// ======================
//   FUNZIONE "MAIN" ROOT
// ======================
//
// Questa è la funzione che chiami da ROOT.
//
// Esempio di uso da shell ROOT:
//   .L WaveAnalysis.C+
//   WaveAnalysis("wave0.txt", "wave1.txt", "wave2.txt", "wave3.txt",
//                200.0, 5000.0);
//
// dove:
//   minStepAmp    = 200.0 (cut su max-min)
//   maxAmpForHist = 5000.0 (range asse x istogramma di amp)

void WaveAnalysis(const char* file0,
                  const char* file1,
                  const char* file2,
                  const char* file3,
                  double minStepAmp    = 200.0,
                  double binWidth       = 100.0)
                  {
                    auto start = std::chrono::high_resolution_clock::now();
    // 1) Leggi i 4 file e costruisci la map<eventNumber, EventWave>
    //auto capisce in automatico il tipo di variabile in questo caso viene compilato come std::map<int, EventWave>
    auto events0 = ReadWaveFile(file0); // legge un file di waveforms e li ritorna in una mappa
    auto events1 = ReadWaveFile(file1);
    auto events2 = ReadWaveFile(file2);
    auto events3 = ReadWaveFile(file3);

 // 2) Crea una struttura con i 4 canali + somma per ogni eventNumber comune
    std::map<int, MultiChannelEvent> merged; //merge degli eventi

    for (const auto &kv : events0) {
        int ev = kv.first; // estrae il numero dell'evento  

        // Verifica che lo stesso eventNumber esista in tutti e 4 i file
        if (events1.count(ev) == 0 ||
            events2.count(ev) == 0 ||
            events3.count(ev) == 0)
            continue;
        // Prendi l’evento numero ev dalla mappa events0 e creane un riferimento costante con nome w0
        const auto &w0 = events0.at(ev); // vettore dei campioni ADC di wave0 → canale 0
        const auto &w1 = events1.at(ev); 
        const auto &w2 = events2.at(ev);
        const auto &w3 = events3.at(ev);

        // per sicurezza prendo la lunghezza minima tra i canali
        size_t N = std::min({w0.ch.size(), w1.ch.size(), w2.ch.size(), w3.ch.size()});      // numero di campioni da considerare

        //m è un multichannel event la struttura dichiarata ad inzio codice che ha come argomenti i 4 canali e la somma e eventNumber
        MultiChannelEvent m; // crea un nuovo MultiChannelEvent
        m.eventNumber = ev; // imposta il numero dell'evento
        m.ch[0].assign(w0.ch.begin(), w0.ch.begin() + N); // copia i campioni del canale 0
        m.ch[1].assign(w1.ch.begin(), w1.ch.begin() + N); // copia i campioni del canale 1
        m.ch[2].assign(w2.ch.begin(), w2.ch.begin() + N); // copia i campioni del canale 2
        m.ch[3].assign(w3.ch.begin(), w3.ch.begin() + N); // copia i campioni del canale 3
        m.sum.resize(N);

        // Calcola la somma dei 4 canali
        for (size_t i = 0; i < N; ++i) {
            m.sum[i] = static_cast<unsigned int>(m.ch[0][i]) + //faccio il casting a int per sicurezza di non avere overflow
                       static_cast<unsigned int>(m.ch[1][i]) + // gli unsigned int short sono al massimo 65535 quindi la somma di 4 valori potrebbe superare questo limite
                       static_cast<unsigned int>(m.ch[2][i]) +
                       static_cast<unsigned int>(m.ch[3][i]);
        }

        merged[ev] = std::move(m); // aggiungi l'evento combinato alla mappa "merged"
    }

    // 3) A questo punto "merged" contiene, per ogni eventNumber comune ai 4 file:
    //    - m.ch[c][i] con c=0..3, i=0..N-1
    //    - m.sum[i] somma dei 4 canali
    //
    // Qui puoi:
    //  - fare grafici
    //  - salvare su TTree
    //  - calcolare l'altezza del gradino, ecc.

    std::cout << "Eventi combinati (presenti in tutti e 4 i file): "
              << merged.size() << std::endl;

    // Esempio: stampa per il primo evento trovato
    if (!merged.empty()) {
        const auto &first = merged.begin()->second;
        std::cout << "Primo eventNumber combinato: " << first.eventNumber << "\n";
        std::cout << "Numero di campioni: " << first.sum.size() << "\n";
        std::cout << "Campioni [i, ch0, ch1, ch2, ch3, sum] per i=0..4:\n";
        for (size_t i = 0; i < std::min<size_t>(5, first.sum.size()); ++i) {
            std::cout << i << "  "
                      << first.ch[0][i] << "  "
                      << first.ch[1][i] << "  "
                      << first.ch[2][i] << "  "
                      << first.ch[3][i] << "  "
                      << first.sum[i] << "\n";
        }
    }

    // Informo l'utente su quanti eventi ci sono
    std::cout << "[BrowseSummedEvents] Numero di eventi disponibili: "
              << merged.size() << std::endl;



              //se volessi controllare il dataset che sto usando abilitare questa funzione
              //questa funzione apre un browser interattivo per vedere gli eventi uno alla volta
              //purtroppo il plot del primo grafico rimane vuoto
              //con l'inserimento del primo Y/y appare il primo grafico correttamente
              //una volta dato nuovamente y come input il canvas verrà chiuso e ne apparirà uno nuovo con il secondo evento e così via


    //BrowseSummedEventsSingleCanvas(merged, 1.0);  // dt = 1.0 → asse x = indice


    // funzione per l'analisi dei gradini
    // questa funzione crea un istogramma con la distribuzione delle ampiezze dei gradini
    //le ampiezze vengono calcolate facendo il fit con una sigmoide per ogni evento
    //si scartano gli eventi con ampiezza minore di minStepAmp
    //l'istogramma va da minAmp a maxAmp + 1000
    //e  ha come numero di bin n_bin = (maxAmp - minAmp) / binWidth
    //è molto lento abbiate pazienza non è ottimizzato per un cazzo
   
    AnalyzeSteps(merged, minStepAmp, binWidth);
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    std::cout << "Tempo totale di esecuzione: " << elapsed << " secondi.\n";

}
 
