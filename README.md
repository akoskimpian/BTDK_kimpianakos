# Kimpián Ákos BTDK dolgozatához felhasznált scriptjei


Előkészületek: 
A következő eszközök szükségesek a bash scriptek lefuttatásához: <br> wiggletools, wigToBigWig, bigWigToBedGraph, fseq2, bedtools

### A "scripts" mappában található kódokat használtam fel elemzésem során az itteni felsorolásuk sorrendjében.

bigwig_median.sh:
  &emsp;A felhasználandó bigwig formátumú fájlokból egy wig fájlt készít, ami a genom minden pozícióban a felhasznált bigwig fájlokból számolt mediánt tartalmazza.

filter_and_convert_median.sh:
  &emsp;Eltávolítja a véletlenszerű (random) scaffoldokat a wig fájlból, majd a fájlt BigWig formátumba konvertálja. Az eredmény a konszenzus fájl.
  
chunk_creation.sh:
  &emsp;Számítógépem kapacitáskorlátai és a fájl nagy mérete miatt a bigwig fájlt 10 millió sort tartalmazó "chunk"-okra osztottam.

chunk_callpeak.sh:
  &emsp;Végigiterál az összes chunk-on, majd elvégzi a "peak calling" folyamatát, az eredmény az egyes chunkokhoz tartozó, a csúcsrégiókat tartalmazó bed fájlok.

fseq_results_merge.sh:
  &emsp;A chunk-ok peakjeit egyetlen bed fájlba egyesíti.

open_and_in_promoter.sh:
  &emsp;Átfedéseket azonosít promóter régiókkal, majd ezeket az átfedő régiókat output-olja.

tf_bs_how_many_unique_new_blacklist.sh:
  &emsp;Eltávolítja a "feketelistás" régiókat, majd a promóterekben található nyitott régiókban (az előző script output-ja) 199 transzkripciós faktor kötőhelyeit keresi. Az output az, hogy az egyes vizsgált transzkripciós faktorok hány különböző gén promóteréhez kötnek, azaz hány TF-célgén interakció valósulhat meg a vizsgált minták alapján.
