# Walidacja wyniku — Blok 4

# 15.1 
- Liczba wierszy nagłówkowych: 4
- Liczba wariantów: 1543
- Unikalna liczba kolumn: 103
- Liczba powtórzonych markerów (ID): 0
# 15.2
- bcftools view final.vcf > /dev/null: plik poprawny formalnie
- number of samples: 94
- number of records / number of SNPs: 1543
- no-ALTs / MNPs / indels / others / multiallelic sites: 0
- Ts/Tv: 1.70
# 15.3
- genom = DArT REF: 367
- genom = DArT ALT: 338
- poza {REF,ALT}: 838 (54.3%)
- spójnych (REF lub ALT): 45.7%
# Wniosek: udział poza {REF, ALT} wynosi 54.3% (838 wierszy), co według kryteriów kontroli sugeruje możliwy błąd pozycji lub korekty nici minus / różnice w danych i wymaga weryfikacji.
# Ćwiczenie 4.3
Przy stałym start = 26 kod dalej działa poprawnie składniowo i przypisuje poprawne etykiety genotypów, zmieniając jedynie ich rozkład w wyniku. Nie pojawia się żaden komunikat błędu, bo awk nie pilnuje, czy czytana kolumna odpowiada właściwej próbce.
#Ćwiczenie 4.4
Udział rekordów poza {REF, ALT} wzrósł z 49.4% w step1_core.tsv do 54.3% po korekcie nici w step2_core_rc.tsv. Po ujednoliceniu orientacji z nicią referencyjną udział ten powinien spadać, więc odwrotny kierunek zmiany wskazuje na potrzebę weryfikacji logiki korekty nici lub orientacji w danych.
