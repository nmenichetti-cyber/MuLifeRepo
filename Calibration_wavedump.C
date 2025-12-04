void analizza_array() {
    // Definizione dei 4 array di esempio
    // Legge i dati dai file esterni
    Int_t size1 = 0, size2 = 0, size3 = 0, size4 = 0;
    const int n = size1; // Dimensione massima degli array
    Double_t array1[n], array2[n], array3[n], array4[n];

    // Legge il primo file
    ifstream file1("array1.txt");
    if(!file1.is_open()) {
        cout << "Errore: impossibile aprire array1.txt" << endl;
        return;
    }
    while(file1 >> array1[size1] && size1 < n) {
        size1++;
    }
    file1.close();
    cout << "Letti " << size1 << " valori da array1.txt" << endl;
    
    // Legge il secondo file
    ifstream file2("array2.txt");
    if(!file2.is_open()) {
        cout << "Errore: impossibile aprire array2.txt" << endl;
        return;
    }
    while(file2 >> array2[size2] && size2 < n) {
        size2++;
    }
    file2.close();
    cout << "Letti " << size2 << " valori da array2.txt" << endl;
    
    // Legge il terzo file
    ifstream file3("array3.txt");
    if(!file3.is_open()) {
        cout << "Errore: impossibile aprire array3.txt" << endl;
        return;
    }
    while(file3 >> array3[size3] && size3 < n) {
        size3++;
    }
    file3.close();
    cout << "Letti " << size3 << " valori da array3.txt" << endl;
    
    // Legge il quarto file
    ifstream file4("array4.txt");
    if(!file4.is_open()) {
        cout << "Errore: impossibile aprire array4.txt" << endl;
        return;
    }
    while(file4 >> array4[size4] && size4 < n) {
        size4++;
    }
    file4.close();
    cout << "Letti " << size4 << " valori da array4.txt" << endl << endl;
    Double_t array = array1[] + array2[] + array3[] + array[]  
    

    // Vettori per salvare le differenze
    vector<Double_t> diff;
    vector<Int_t> index;
    
    // Salva i primi valori di riferimento iniziali
    Double_t ref = array[0];
    
    cout << "Valori di riferimento iniziali:" << endl;
    cout << "Array: " << ref1 << endl;
    
    // Scorre gli array e salva le differenze > 2000, resettando quando torna al valore iniziale
    for(int i = 1; i < n; i++) {
        // Array
        if(array[i] == ref) {
            cout << "Array: valore tornato a " << ref1 << " all'indice " << i << " - RESET" << endl;
        } else {
            Double_t d = array[i] - ref;
            if(d > 2000) {
                diff.push_back(d);
                index.push_back(i);
            }
        }
        }
    }
    
    cout << endl;
    double mean = 0;
    for(int j = 1; j<n; j++){
    mean = (mean + array[i])/n; 
    }

    Scale_factor = 60/mean;

    // Stampa i risultati
    cout << "Il valore è: " << Scale_factor << endl;
    
