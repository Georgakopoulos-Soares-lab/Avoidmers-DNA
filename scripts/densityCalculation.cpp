#include <iostream>
#include <fstream>
#include <set>

using namespace std;

int isAba(const string& x){
    int N, mid;
    N = x.length();

    if (N <= 2) {
        return 0;
    }

    mid = (N+1)/2;
    for(int i=1; i<mid; i++) {
        if (x.compare(0, i, x, N - i, i) == 0) {
            return 1; 
        }
    }

    return 0;
}

int findZiminDensity(string sequence) {
    int N, totalAba;
    string chunk;
    N = sequence.length();
    totalAba = 0;

    for(int l=1; l<=N; l++) {
        for(int i=0; i<N-l+1; i++) {
                chunk = sequence.substr(i, l);
                if (isAba(chunk)==1) {
                    totalAba++;
                }
        }
    }
    return totalAba;
}

int main() {
    string line;
    double ziminDensity, kmerDiversity;
    int totalSubsequences, totalAba, N;

    ifstream file("zimin_sequences.txt");
    ofstream fout("zimin_enriched.txt");

    if (!file.is_open()) {
        std::cerr << "Could not open the file." << std::endl;
        return 1;
    }

    while(getline(file, line)) {
        N = line.length();
        totalSubsequences = N * (N+1) / 2;
        totalAba = findZiminDensity(line);
        ziminDensity = totalAba / double(totalSubsequences);
        fout << line << "\t" << totalAba << "\t" << ziminDensity << "\n";
    }

    file.close();
    fout.close();
    return 0;
}
