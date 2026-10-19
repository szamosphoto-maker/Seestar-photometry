# Seestar Photometry – Windows EXE (saját használatra)

A Python fotometriai motorja az eredeti `seestar_photometry(1).py` fájlból származik, változatlan számításokkal.
A `seestar_photometry_gui.py` külön kezelőfelületet ad hozzá.

## Privát GitHub repository

1. Hozz létre új repositoryt, állítsd **Private** értékre.
2. Töltsd fel a csomag tartalmát. A `build-windows.yml` a `.github/workflows/` könyvtárba kerüljön (a ZIP-ben már ott van).
3. GitHub **Actions → Build Seestar Photometry Windows EXE → Run workflow**.
4. Sikeres futás után töltsd le a **Seestar_Photometry_Windows** artifactot.
5. Csomagold ki, és indítsd a `Seestar_Photometry.exe`-t. Python telepítése nem szükséges.

## Használat

- Válaszd ki a CSV-k és a `Csillagnév.txt` adatlapok közös könyvtárát.
- Válassz célcsillagot, majd indítsd a feldolgozást.
- A program **nem olvas almappákat**, az eredményeket a könyvtár `Eredmenyek` almappájába menti.
- A keresési sugár alapértelmezése 8 ívmásodperc.
- A nyers, tisztított és fázisfénygörbék, illetve a CSV-k és a minőségjelentés az eredeti motor szerint készülnek.

## Megjegyzés

A GUI fordítása GitHubon történik; a build éles Windows-futását az első artifact letöltése után ellenőrizni kell.
