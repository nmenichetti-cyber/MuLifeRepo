#!/usr/bin/env python3

def merge_by_reset_period(file1, file2, output_file):
    """
    Unisce due file allineando i reset_count.
    Gestisce correttamente gli eventi prima del primo reset.
    """
    
    # 1. Trova l'ultimo reset_count del file1
    last_reset1 = -1  # Inizializza a -1 per gestire caso senza reset
    with open(file1, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            channel = int(parts[0])
            # Controlla se è una riga di reset (bit31=1)
            if (channel >> 31) & 1:
                last_reset1 = int(parts[1])
    
    # Se non ci sono reset nel file1, assumi 0
    if last_reset1 == -1:
        last_reset1 = 0
    
    # 2. Trova il PRIMO reset_count del file2
    first_reset2 = None
    events_before_first_reset = []  # Eventi prima del primo reset nel file2
    
    with open(file2, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            channel = int(parts[0])
            counter = int(parts[1])
            
            # Controlla se è una riga di reset
            if (channel >> 31) & 1:
                if first_reset2 is None:
                    first_reset2 = int(parts[1])
                # Salva la riga di reset per dopo
                events_before_first_reset.append((channel, counter, True))
            else:
                # Se non abbiamo ancora trovato il primo reset,
                # questi eventi vanno dopo il file1 ma PRIMA del primo reset del file2
                events_before_first_reset.append((channel, counter, False))
    
    # 3. Calcola l'offset per i reset_count del file2
    # Il file2 continua dal reset_count successivo all'ultimo del file1
    reset_offset = (last_reset1 + 1) - first_reset2 if first_reset2 is not None else 0
    
    print(f"Ultimo reset file1: {last_reset1}")
    print(f"Primo reset file2: {first_reset2}")
    print(f"Reset offset: {reset_offset}")
    
    # 4. Scrivi il file unito
    with open(output_file, 'w') as out:
        # 4a. Scrivi tutto il file1 così com'è
        with open(file1, 'r') as f1:
            out.write(f1.read())
        
        # 4b. Scrivi gli eventi del file2 PRIMA del primo reset (se ci sono)
        # Questi vanno DOPO il file1 ma PRIMA del primo reset del file2
        for channel, counter, is_reset in events_before_first_reset:
            if is_reset:
                if first_reset2 is not None and counter == first_reset2:
                    # Questo è il primo reset: applica l'offset
                    new_counter = counter + reset_offset
                    out.write(f"{channel} {new_counter}\n")
                else:
                    # Altri reset prima del primo? Non dovrebbero esistere
                    # ma se ci sono, applicali comunque
                    new_counter = counter + reset_offset
                    out.write(f"{channel} {new_counter}\n")
            else:
                # Eventi normali prima del primo reset
                out.write(f"{channel} {counter}\n")
        
        # 4c. Continua con il resto del file2 (dopo il primo reset)
        # Riapriamo il file2 e saltiamo la parte già scritta
        if first_reset2 is not None:
            already_written = len(events_before_first_reset)
            with open(file2, 'r') as f2:
                # Salta le righe già scritte
                for _ in range(already_written):
                    next(f2)
                
                # Scrivi il resto applicando l'offset a TUTTI i reset
                for line in f2:
                    parts = line.strip().split()
                    if len(parts) < 2:
                        out.write(line)
                        continue
                    
                    channel = int(parts[0])
                    counter = int(parts[1])
                    
                    if (channel >> 31) & 1:
                        new_counter = counter + reset_offset
                        out.write(f"{channel} {new_counter}\n")
                    else:
                        out.write(line)
    
    print(f"Merge completato!")
    print(f"Eventi prima del primo reset in file2: {len([e for e in events_before_first_reset if not e[2]])}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Uso: python3 merge_muon.py <file1> <file2> <output>")
        print("Esempio: python3 merge_muon.py dati1.txt dati2.txt unito.txt")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output = sys.argv[3]
    
    merge_by_reset_period(file1, file2, output)
