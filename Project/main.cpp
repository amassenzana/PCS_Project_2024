#include "DFNLibrary.hpp"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
    if(argc<2) cerr << "Manca il parametro path" << endl;
    string path = argv[1];

    DFNLibrary::DFN dfn;

    if(!importDFN(path, dfn)) cerr << "Errore nell'importazione dei dati" << endl;
    else cout << "Importazione dati eseguita correttamente" << endl;


    dfn.computeDFN();

    cout << "Programma terminato 0" << endl;
    return 0;
}
