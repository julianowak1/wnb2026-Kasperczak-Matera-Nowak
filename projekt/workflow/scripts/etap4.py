import csv

zmapowane = {}
with open(snakemake.input.snpy, "r", encoding="utf-8") as plik:
    czytnik = csv.reader(plik, delimiter="\t")
    naglowek = next(czytnik)  # pomijamy pierwszy wiersz z nazwami kolumn
    for wiersz in czytnik:
        if not wiersz:
            continue
        id_markera = wiersz[0]
        chrom = wiersz[1]
        poz = wiersz[2]
        nic = wiersz[3]
        ref_dart = wiersz[4]
        alt_dart = wiersz[5]
        zmapowane[id_markera] = [chrom, poz, nic, ref_dart, alt_dart]

with open(snakemake.input.dart_csv, "r", encoding="utf-8") as plik:
    wiersze_dart = list(csv.reader(plik))

nr_wiersza_naglowka = int(snakemake.params.header_row) - 1
nr_pierwszego_wiersza = int(snakemake.params.data_start_row) - 1
nr_kolumny_id = int(snakemake.params.id_col) - 1
nr_kolumny_probek = int(snakemake.params.geno_start_col) - 1

nazwy_probek = []
for nazwa in wiersze_dart[nr_wiersza_naglowka][nr_kolumny_probek:]:
    nazwy_probek.append(nazwa.strip())

liczba_probek = len(nazwy_probek)

with open(snakemake.output[0], "w", encoding="utf-8", newline="") as wyjscie:
    zapis = csv.writer(wyjscie, delimiter="\t")
    zapis.writerow(["marker_id", "chrom", "pos", "strand", "ref_dart", "alt_dart"] + nazwy_probek)

    i = nr_pierwszego_wiersza
    while i + 1 < len(wiersze_dart):
        wiersz_ref = wiersze_dart[i]
        wiersz_snp = wiersze_dart[i + 1]

        id_markera = wiersz_ref[nr_kolumny_id].strip()

        if id_markera in zmapowane:
            dane_markera = zmapowane[id_markera]
            wartosci_ref = wiersz_ref[nr_kolumny_probek : nr_kolumny_probek + liczba_probek]
            wartosci_snp = wiersz_snp[nr_kolumny_probek : nr_kolumny_probek + liczba_probek]

            genotypy = []
            for j in range(liczba_probek):
                r = wartosci_ref[j].strip()
                s = wartosci_snp[j].strip()

                if r == "1" and s == "0":
                    gt = "0/0"
                elif r == "1" and s == "1":
                    gt = "0/1"
                elif r == "0" and s == "1":
                    gt = "1/1"
                else:
                    gt = "./."

                genotypy.append(gt)

            zapis.writerow([id_markera] + dane_markera + genotypy)
        i = i + 2
