import csv

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def odwroc(baza):
    return COMP.get(baza.upper(), baza.upper())


def main():
    blast_plik = snakemake.input.blast
    snp_plik = snakemake.input.snp_info
    wyjscie = snakemake.output[0]

    min_pident = float(snakemake.params.min_pident)
    min_aln_len = int(snakemake.params.min_aln_len)

    snpy = {}
    with open(snp_plik, "r", encoding="utf-8") as f:
        for linia in f:
            linia = linia.strip()
            if not linia:
                continue
            czesci = linia.split("\t")
            if len(czesci) < 2 or ":" not in czesci[1] or ">" not in czesci[1]:
                continue
            marker_id = czesci[0].strip()
            poz_str, allele_str = czesci[1].strip().split(":", 1)
            ref_sonda, alt_sonda = allele_str.split(">", 1)
            snpy[marker_id] = (int(poz_str), ref_sonda.strip(), alt_sonda.strip())

    najlepsze = {}
    with open(blast_plik, "r", encoding="utf-8") as f:
        for linia in f:
            linia = linia.strip()
            if not linia:
                continue
            kol = linia.split("\t")
            qseqid = kol[0]
            sseqid = kol[1]
            pident = float(kol[2])
            dlugosc = int(kol[3])
            qstart = int(kol[6])
            qend = int(kol[7])
            sstart = int(kol[8])
            send = int(kol[9])
            evalue = float(kol[10])
            bitscore = float(kol[11])

            if pident < min_pident or dlugosc < min_aln_len:
                continue

            wynik = (bitscore, pident, -evalue)
            if qseqid not in najlepsze or wynik > najlepsze[qseqid]["wynik"]:
                najlepsze[qseqid] = {
                    "wynik": wynik,
                    "chrom": sseqid,
                    "qstart": qstart,
                    "qend": qend,
                    "sstart": sstart,
                    "send": send,
                    "pident": pident,
                    "dlugosc": dlugosc,
                }

    with open(wyjscie, "w", encoding="utf-8", newline="") as out:
        pisz = csv.writer(out, delimiter="\t")
        pisz.writerow(["marker_id", "chrom", "pos", "strand", "ref_dart", "alt_dart", "pident", "aln_len"])

        for marker_id, hit in najlepsze.items():
            if marker_id not in snpy:
                continue
            poz_w_sondzie, ref_s, alt_s = snpy[marker_id]
            if not (hit["qstart"] <= poz_w_sondzie <= hit["qend"]):
                continue

            przesuniecie = poz_w_sondzie - hit["qstart"]
            if hit["sstart"] <= hit["send"]:
                nic = "+"
                poz_genom = hit["sstart"] + przesuniecie
                ref_out, alt_out = ref_s, alt_s
            else:
                nic = "-"
                poz_genom = hit["sstart"] - przesuniecie
                ref_out, alt_out = odwroc(ref_s), odwroc(alt_s)

            pisz.writerow([
                marker_id,
                hit["chrom"],
                poz_genom,
                nic,
                ref_out,
                alt_out,
                hit["pident"],
                hit["dlugosc"],
            ])


if __name__ == "__main__":
    main()
