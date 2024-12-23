#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <stdio.h>
#include <ctime>
#include <chrono>
#include <map>

using namespace std;
int NMAX = 100000;

// Random number in ranger
int getRandomInRange(int min, int max) {
    // Initialize random number generator with a random seed
    random_device rd;  
    mt19937 gen(rd()); // Mersenne Twister engine seeded with rd()
    
    // Define the range [min, max]
    uniform_int_distribution<> dist(min, max);

    // Generate and return a random number in the specified range
    return dist(gen);
}


class model {
    public:
    int L;
    int N;
    float h;
    vector<int> psi;
    vector<int> tau;
    vector<int> plaquettes;

    random_device rd;  // Used to obtain a seed for the random number engine
    mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<int> rand_x;
    uniform_int_distribution<int> rand_y;
    uniform_int_distribution<int> rand_z;
    uniform_int_distribution<int> rand_dir;
    uniform_real_distribution<float> rand_prob;

    model(int L, double h) {
        this->L = L;
        this->N = 3*pow(L,3);
        this->psi = vector<int> (3*pow(L,3), 1);
        this->plaquettes = vector<int> (3*pow(L,3), 1);
        this->tau = vector<int>(pow(L,3),1);

        this->h = h;

        this->gen = std::mt19937(rd());
        this->rand_x = std::uniform_int_distribution<int>(0, L-1);
        this->rand_y = std::uniform_int_distribution<int>(0, L-1);
        this->rand_z = std::uniform_int_distribution<int>(0, L-1);
        this->rand_dir = std::uniform_int_distribution<int>(0, 2);
    }

    int getIdx(int x, int y, int z, int d) {
        return 3 * (x + this->L * y + this->L * this->L * z) + d;
    }

    int getTauIdx(int x, int y, int z) {
        return x + this->L * y + this->L * this->L * z;
    }

    void flipVertex(int x, int y, int z) {
        this->tau[getTauIdx(x,y,z)] *= -1;
        int L = this->L;

        this->psi[getIdx(x,y,z,0)]              *= -1;
        this->psi[getIdx((x+L-1)%L, y, z, 0)]   *= -1;
        this->psi[getIdx(x,y,z,1)]              *= -1;
        this->psi[getIdx(x,(y+L-1)%L, z, 1)]    *= -1;
        this->psi[getIdx(x,y,z,2)]              *= -1;
        this->psi[getIdx(x, y, (z+L-1)%L, 2)]   *= -1;
    }

    void sweep (double T) {
        float beta = 1/T;
        
        for (int n = 0; n < this->N; n++) {
            int x = rand_x(gen);
            int y = rand_y(gen);
            int z = rand_z(gen);
            int d = rand_dir(gen);

            int s = getIdx(x,y,z,d);
            int f1 = 0;
            int f2 = 0;
            int f3 = 0;
            int f4 = 0;

            int t1 = getTauIdx(x,y,z);
            int t2 = 0;

            if (d == 0) {
                f1 = getIdx(x,y,z,d);
                f2 = getIdx(x, (y + this->L - 1) % (this->L),z,d);
                f3 = getIdx(x,y,z,d+2);
                f4 = getIdx(x, y, (z + this->L - 1) % (this->L),d+2);
                t2 = getTauIdx((x+1)%this->L,y,z);

            }

            else if (d == 1) {
                f1 = getIdx(x,y,z,d-1);
                f2 = getIdx((x + this->L - 1) % (this->L),y,z,d-1);
                f3 = getIdx(x,y,z,d);
                f4 = getIdx(x, y, (z + this->L - 1) % (this->L),d);
                t2 = getTauIdx(x,(y+1)%this->L,z);
            }

            else {
                f1 = getIdx(x,y,z,d-1);
                f2 = getIdx(x, (y + this->L - 1) % (this->L),z,d-1);
                f3 = getIdx(x,y,z,d);
                f4 = getIdx((x + this->L - 1) % (this->L), y, z,d);
                t2 = getTauIdx(x,y,(z+1)%this->L);
            }

            float dS = 2 * (this->plaquettes[f1] + 
                            this->plaquettes[f2] + 
                            this->plaquettes[f3] + 
                            this->plaquettes[f4] + 
                            this->h * this->psi[s] * this->tau[t1] * this->tau[t2]); 

            if (dS < 0){
                    this->psi[s] *= -1;
                    this->plaquettes[f1] *= -1;
                    this->plaquettes[f2] *= -1;
                    this->plaquettes[f3] *= -1;
                    this->plaquettes[f4] *= -1;
            }
            else {
                float p = exp(-beta * abs(dS));
                if (rand_prob(gen) < p) {
                    this->psi[s] *= -1;
                    this->plaquettes[f1] *= -1;
                    this->plaquettes[f2] *= -1;
                    this->plaquettes[f3] *= -1;
                    this->plaquettes[f4] *= -1;
                }
            }

            if (rand_prob(gen) < 0.5) {
                flipVertex(rand_x(gen), rand_y(gen), rand_z(gen));
            }


        }

        
        for (int n = 0; n < pow(L,3); n++) {
            int x = rand_x(gen);
            int y = rand_y(gen);
            int z = rand_z(gen);

            float energy = 0;

            int t = getTauIdx(x,y,z);

            int s1 = getIdx(x,y,z,0);
            int t1 = getTauIdx((x+1)%this->L,y,z);

            energy += this->h * this->tau[t] * this->tau[t1] * this->psi[s1];

            int s2 = getIdx(x,y,z,1);
            int t2 = getTauIdx(x,(y+1)%this->L, z);
            energy += this->h * this->tau[t] * this->tau[t2] * this->psi[s2];

            int s3 = getIdx(x,y,z,2);
            int t3 = getTauIdx(x, y, (z+1)%this->L);
            energy += this->h * this->tau[t] * this->tau[t3] * this->psi[s3];

            int s4 = getIdx((x+this->L-1)%this->L,y,z,0);
            int t4 = getTauIdx((x+this->L-1)%this->L,y,z);
            energy += this->h * this->tau[t] * this->tau[t4] * this->psi[s4];

            int s5 = getIdx(x,(y+this->L-1)%this->L,z,1);
            int t5 = getTauIdx(x,(y+this->L-1)%this->L,z);
            energy += this->h * this->tau[t] * this->tau[t5] * this->psi[s5];

            int s6 = getIdx(x,y,(z+this->L-1)%this->L,2);
            int t6 = getTauIdx(x,y,(z+this->L-1)%this->L);
            energy += this->h * this->tau[t] * this->tau[t6] * this->psi[s6];

            float dS = 2 * energy;

            if (dS < 0){
                this->tau[t] *= -1;
            }
            else {
                float p = exp(-beta * abs(dS));
                if (rand_prob(gen) < p) {
                    this->tau[t] *= -1;
                }
            }
            
        }
        

    }

    float calculateEnergy() {
        float energy = 0;

        for (int z = 0; z < this-> L; z++) {
            for (int y = 0; y < this->L; y++) {
                for (int x = 0; x < this-> L;x++) {
                    int t0 = getTauIdx(x,y,z);
                    int t1 = getTauIdx((x+1)%this->L,y,z);


                    int s1 = getIdx(x,y,z,0);

                    int t2 = getTauIdx(x,(y+1)%this->L,z);
                    int s2 = getIdx(x,y,z,1);

                    int t3 = getTauIdx(x,y,(z+1)%this->L);
                    int s3 = getIdx(x,y,z,2);

                    energy += this->h * this->psi[s1] * this->tau[t0] * this->tau[t1];
                    energy += this->h * this->psi[s2] * this->tau[t0] * this->tau[t2];
                    energy += this->h * this->psi[s3] * this->tau[t0] * this->tau[t3];


                }
            }
        }

        for (auto i : this->plaquettes)
            energy += i;
        
        return energy;
    }

};



int main(int argc, char* argv[]) {

    time_t now = time(0);
    tm *ltm = localtime(&now);

    // Simulation Params
    int L = 8;
    float T = 1;
    float h = 0.2;
    int maxl = 20;
    int binsteps = 1;
    int steps = 1000 * maxl;
    string outputDir = "./";

    for (int i = 0; i < argc; i++) {
        string arg = "";
        arg += argv[i];
        if (arg.substr(0,2).compare("-T") == 0) {
            T = stof(arg.substr(3));
        }
        else if (arg.substr(0,2).compare("-L") == 0) {
            L = stoi(arg.substr(3));
        }
        else if (arg.substr(0,2).compare("-h") == 0) {
            h = stof(arg.substr(3));
        }
        else if (arg.substr(0,5).compare("-maxl") == 0) {
            maxl = stoi(arg.substr(6));
        }
        else if (arg.substr(0,9).compare("-binsteps") == 0) {
            binsteps = stoi(arg.substr(10));
        }
        else if (arg.substr(0,10).compare("-outputDir") == 0) {
            outputDir = arg.substr(11) + "/";
        }

    }

    vector<float> E;

    model system(L, h);
    
    // equilibrate
    for (int i = 0; i < 1000; i++) {
        system.sweep(T);
        cout << i << '\r' << flush;
    }


    
    cout << "equilibrated" << flush;
    cout << "" << endl;

    for (int i = 0; i < steps; i++) {
        cout << i << '\r' << flush;
        system.sweep(T);
        if (i % maxl == 0) {
        float energy = system.calculateEnergy();
        E.push_back(energy);
        }
    }

    string filename = outputDir + "3DGauge_MC_taus_L=" + to_string(L) + "_T=" +
        to_string(T).substr(0,4) + "_h="+ to_string(h).substr(0,4) + "_steps=" + to_string(steps) + "_maxl=" + to_string(maxl)  + "_binsteps="
        + to_string(binsteps) + "_" + to_string(1 + ltm->tm_mon) + "_" + to_string(ltm->tm_mday) + "_" + to_string(1900 + ltm->tm_year) + ".txt";

    ofstream outFile(filename);

    if (outFile.is_open()) {
        // Write each vector to the file
        for (const auto& vec : {E}) {
            for (const auto& val : vec) {
                outFile << val << " ";
            }
            outFile << "\n";  // New line after each vector
        }

        outFile.close();
        std::cout << "Data has been written to ";
        cout << filename << endl;
    } else {
        std::cerr << "Unable to open file for writing.\n";
    }
    

    return 0;
}