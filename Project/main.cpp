#include "DFNLibrary.hpp"
#include <iostream>

using namespace std;

int main(int argc, char** argv){
    string path = argv[1];

    DFNLibrary::DFN dfn;

    if(!importDFN(path, dfn))
        cerr << "Errore" << endl;

    cout << "Importazione dati eseguita correttamente" << endl;


    return 0;
}
