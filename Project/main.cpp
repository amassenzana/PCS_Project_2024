#include "DFNLibrary.hpp"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
    if(argc<2) {
        cerr << "Manca il parametro path" << endl;
        return -1;
    }
    string path = argv[1];

    DFNLibrary::DFN dfn;

    if(!importDFN(path, dfn)) {
        cerr << "Errore nell'importazione dei dati" << endl;
        return -2;
    }
    else cout << "Importazione dati eseguita correttamente" << endl;

    dfn.computeDFN();
    dfn.output();
    dfn.plotFracture();

    cout << "Programma terminato 0" << endl;
    return 0;
}
