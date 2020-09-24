// this is madness // yes this is indeed madness

#include <iostream>
#include <fstream> // file
#include <vector> // vector
#include <valarray> // numeric vector
#include <sstream> // stringstream
#include <algorithm> // l3iba f j3iba
#include <numeric> // iota
#include <queue> // queue
#include "solver.h" // matrix solver

using namespace std;

// reads file and initializes the conectivity_vector with all information
void init(fstream& file_input, int& gates, int& nets, int& pads, vector<vector<int>>& conectivity_vector, vector<vector<double>>& pads_vector_coords, vector<int>& pad_netlists) {

    file_input >> gates;
    file_input >> nets;

    int var;
    double var_cord;
    conectivity_vector.resize(nets);
    stringstream s;
    string line;

    getline(file_input, line);

    for (int i = 0; i < gates; i++) {
        getline(file_input, line);
        s.str(line);
        s >> var;
        s >> var;

        while (s >> var) {
            conectivity_vector[var - 1].push_back(i);
        }
        // sorted by default!!!!!!!

        s.str(std::string()); //clear stringstream // starburst string streamu
        s.clear();
    }


    file_input >> pads;
    pads_vector_coords.resize(pads);
    pad_netlists.resize(pads);

    getline(file_input, line);
    for (int i = 0; i < pads; i++) {
        getline(file_input, line);
        s.str(line);
        s >> var;
        // net
        s >> var;
        pad_netlists[i] = (var - 1); // net
        // coords
        s >> var_cord;
        pads_vector_coords[i].push_back(var_cord);
        s >> var_cord;
        pads_vector_coords[i].push_back(var_cord);

        s.str(std::string()); //clear stringstream // starburst string streamu
        s.clear();
    }
}

// write result to gate
void write_result(fstream& file_output, valarray<double> x, valarray<double> y) {
    for (size_t i = 0; i < x.size(); i++) {
        file_output << i + 1 << " " << x[i] << " " << y[i] << "\n";
    }
}

// makes the new pad vector with either the gates that become pads on the frontier or the pads that are translated
vector<double> become_pad(double x, double y, int x_div, int y_div, vector<double> up_left, int acc, int vertical) {

    // acc and vertical are for gates that are inside the partition that need to be moved left or right or up or down
    // I just drew the examples up to 64 partitions and found this relations that I use to partition and find coordinates....
    // they work 

    double new_x, new_y;
    bool inside = (x >= up_left[0]) && (x <= (up_left[0] + 100.0 / x_div)) && (y >= up_left[1]) && (y <= (up_left[1] + 100.0 / y_div)); // pads in frontiere will always be inside

    if(!inside){ // outside
        if (x < up_left[0]) { // left
            new_x = up_left[0];
        }
        else if (x > (up_left[0] + 100.0 / x_div)) { // right
            new_x = up_left[0] + 100.0 / x_div;
        }
        else {
            new_x = x; // between
        }

        if (y < up_left[1]) { // up
            new_y = up_left[1];
        }
        else if (y > (up_left[1] + 100.0 / y_div)) { // down
            new_y = up_left[1] + 100.0 / y_div;
        }
        else {
            new_y = y; // between
        }
    }
    else { // inside
        if (!vertical) {
            if (((acc + 1) % 2) == 0) { // move to right
                new_x = up_left[0];
            }
            else { // move to left
                new_x = up_left[0] + 100.0 / x_div;
            }

            new_y = y; // leave it
        }
        else {
            new_x = x; // leave it

            if (((acc + 1) % 2) == 0) {  // move to up
                new_y = up_left[1];
            }
            else {
                new_y = up_left[1] + 100.0 / y_div; // move to down
            }
        }
    }

    return vector<double>{new_x, new_y};
}

// intersection using binary search, netlist is already sorted by how I read the file, get the gates present in a netlist, the rest is added to gates_to_pads in order to become pad
void find_gates_gates(vector<int> gates, vector<int> netlist, vector<int>& gates_indices, vector<int>& gates_to_pads) {

    // clear from previous iteration
    gates_indices.clear();
    gates_to_pads.clear();

    for (size_t i = 0; i < gates.size(); i++) { // for each gate, see if there is a match in netlist using binary search
        if (binary_search(netlist.begin(), netlist.end(), gates[i])) {
            gates_indices.push_back(i);
            netlist.erase(remove(netlist.begin(), netlist.end(), gates[i]), netlist.end()); // erase the value of the gate from the netlist, still sorted
        }
    }

    // tzz a hmed, u never gonna find the gate inside the gate vector Zzzzzz
    // if there were any gates, all of them got removed, the rest will become pads
    for (size_t i = 0; i < netlist.size(); i++) {
        gates_to_pads.push_back(netlist[i]); // these gates dont have an index because they dont appear on gates vector
    }
}


void make_matricesA_B_V3(vector<vector<int>> conectivity_vector, vector<vector<double>> pads_vector_coords, vector<int> pad_netlists,
    coo_matrix& A, valarray<double>& Bx, valarray<double>& By, vector<int> gates, valarray<double> x, valarray<double> y, int x_div, int y_div, vector<double> up_left, int acc, int vertical) {

    // problem with V1 : too much shit going on and bad, didn't really need to make a conectivity vector each time
    // problem with V2: if some gates are connected again in another and doing columns and rows of C wrongly
    // problem with V3: imagine writing a V4

    vector<int> R, C;
    vector<double> Dat;
    vector<int> gates_to_pads;
    vector<int> gates_indices;
    vector<double> netlist_weight(conectivity_vector.size());
    double temp;
    vector<vector<double>> Matrix(gates.size(), vector<double>(gates.size(), 0.0));
    int pad = 0;

    for (size_t i = 0; i < conectivity_vector.size(); i++) {
        // find gates indices and gates to pads for this netlist
        pad = 0;
        find_gates_gates(gates, conectivity_vector[i], gates_indices, gates_to_pads);
        
        for (size_t j = 0; j < pad_netlists.size(); j++) {
            if (pad_netlists[j] == i)
                pad++;
        }
        pad += gates_to_pads.size();
        if ((gates_indices.size() + pad) == 1) { // error in the data file netlist with only 1 gate and no other conection, does nothing we skip
            continue;
        }
        temp = 1.0 / (gates_indices.size() + pad - 1);
        netlist_weight[i] = temp;

        // fill matrix A with weights and store them
        for (size_t j = 0; j < gates_indices.size(); j++) {
            Matrix[gates_indices[j]][gates_indices[j]] += pad * temp; //  we add the number of the pads to the diagonal element
            for (size_t k = j + 1; k < gates_indices.size(); k++) {
                // diagonal
                Matrix[gates_indices[j]][gates_indices[j]] += temp;
                Matrix[gates_indices[k]][gates_indices[k]] += temp;
                // non diagonal
                Matrix[gates_indices[j]][gates_indices[k]] -= temp;
                Matrix[gates_indices[k]][gates_indices[j]] -= temp;
            }
        }

        // every gate to pad becomes a pad
        for (size_t j = 0; j < gates_to_pads.size(); j++) { // never goes inside when first iteration, 100%
            pad_netlists.push_back(i);
            pads_vector_coords.push_back(become_pad(x[gates_to_pads[j]], y[gates_to_pads[j]], x_div, y_div, up_left, acc, vertical));
        }

        // netlist becomes gate_indices
        conectivity_vector[i] = gates_indices; // we don't need conectivity vector[i] anymore, so we use the (netlist,gates) pair to make the Bx and By // pass by copy
    }

    // Pads and Bx, By
    Bx.resize(gates.size(), 0.0);
    By.resize(gates.size(), 0.0);


    for (size_t i = 0; i < pad_netlists.size(); i++) {
        for (size_t j = 0; j < conectivity_vector[pad_netlists[i]].size(); j++) {
            if ((pads_vector_coords[i][0] >= up_left[0]) && (pads_vector_coords[i][0] <= (up_left[0] + 100.0 / x_div)) &&
                (pads_vector_coords[i][1] >= up_left[1]) && (pads_vector_coords[i][1] <= (up_left[1] + 100.0 / y_div))) { // if the pad inside nothing happens

                Bx[conectivity_vector[pad_netlists[i]][j]] += pads_vector_coords[i][0] * netlist_weight[pad_netlists[i]];
                By[conectivity_vector[pad_netlists[i]][j]] += pads_vector_coords[i][1] * netlist_weight[pad_netlists[i]];
            }
            else { // if it's outside, we need to move it
                pads_vector_coords[i] = become_pad(pads_vector_coords[i][0], pads_vector_coords[i][1], x_div, y_div, up_left, acc, vertical);

                Bx[conectivity_vector[pad_netlists[i]][j]] += pads_vector_coords[i][0] * netlist_weight[pad_netlists[i]];
                By[conectivity_vector[pad_netlists[i]][j]] += pads_vector_coords[i][1] * netlist_weight[pad_netlists[i]];
            }
        }
    }

    // we make R, C and Dat;
    for (size_t i = 0; i < gates.size(); i++) {
        for (size_t j = 0; j < gates.size(); j++) {
            if (Matrix[i][j] != 0.0) {
                R.push_back(i);
                C.push_back(j);
                Dat.push_back(Matrix[i][j]);
            }
        }
    }

    A.n = gates.size();
    A.nnz = R.size();

    A.row.resize(A.nnz);
    A.col.resize(A.nnz);
    A.dat.resize(A.nnz);

    copy(begin(R), end(R), begin(A.row));
    copy(begin(C), end(C), begin(A.col));
    copy(begin(Dat), end(Dat), begin(A.dat));

}

// stores the new values of X and Y of the gates in gates_indices in valarrays x and y
void store(valarray<double>& x, valarray<double>& y, valarray<double> x_temp, valarray<double> y_temp, vector<int> gates_indices) {
    for (size_t i = 0; i < x_temp.size(); i++) {
        x[gates_indices[i]] = x_temp[i];
        y[gates_indices[i]] = y_temp[i];
    }
}

// copies the gates between start and end to temp and returns it
vector<int> half_gates_indices(vector<int> keep_track_indices, vector<int> gates_indices, size_t start, size_t end) {
    vector<int> temp(end - start);
    //copy(begin(gates_indices) + start, begin(gates_indices) + end, begin(temp)); this line was so cool sniff
    for (size_t i = start; i < end; i++) {
        temp[i - start] = gates_indices[keep_track_indices[i]];
    }
    return temp;
}

// computes the up left corner of the partition coordinates
void up_left_corner_coord(vector<double>& up_left, int acc, int v_h_control) {
    // unit test good!
    double x_temp = 50.0, y_temp = 50.0;
    double temp;
    size_t v = 1;

    while (v_h_control != 1) {
        v_h_control /= 2;

        if (v) {
            temp = x_temp;
            x_temp /= 2;
        }
        else {
            temp = y_temp;
            y_temp /= 2;
        }
        v = ~v & 1;
        up_left[v] += temp * (acc >= v_h_control); // partition numbering starts at 1 that why 
        if (acc >= v_h_control) {
            acc -= v_h_control;
        }
        if (acc == 0) {
            acc = 0;
        }
    }
}

// places the gates using a queue for iteration
void place(vector<vector<int>> conectivity_vector, vector<vector<double>> pads_vector_coords, vector<int> pad_netlists, valarray<double>& x, valarray<double>& y, int g, int depth) {

    // queue, each time we divide and do the stuff and keepdoing the stuff until depth reached
    queue<vector<int>> queue_BFS; // queue will have the indices of the gates;

    coo_matrix A;
    valarray<double> Bx, By;

    vector<int> gates(g);
    iota(gates.begin(), gates.end(), 0);

    queue_BFS.push(gates);
    valarray<double> x_temp, y_temp;

    int vertical = 1;
    int v_h_control = 1;
    int acc = 0; // not starting at for i = 0;, then starts counting at 0

    int x_div = 1; // grid : x division
    int y_div = 1; // y division

    // part of correction we need a vector to keep track ofindices bcs pfor example if we have gates 0 and 2, size = 2, but we can't do x[2] where x is a 2 size vector
    vector<int> keep_track_indices;
    vector<double> up_left(2, 0); // P(x,y) of up left corner of the partition

    // BFS
    while(v_h_control < (depth * 2)) {

        cout << acc << "\n";

        gates = queue_BFS.front();

        keep_track_indices.resize(gates.size()); 
        iota(keep_track_indices.begin(), keep_track_indices.end(), 0);

        x_temp.resize(gates.size());
        y_temp.resize(gates.size());

        up_left_corner_coord(up_left, acc, v_h_control);
        make_matricesA_B_V3(conectivity_vector, pads_vector_coords, pad_netlists, A, Bx, By, gates, x, y, x_div, y_div, up_left, acc, vertical);

        // solve
        x_temp.resize(gates.size());
        y_temp.resize(gates.size());
        A.solve(Bx, x_temp);
        A.solve(By, y_temp);

        // store the new coordinates in x and y
        store(x, y, x_temp, y_temp, gates);

        // sort
        sort(keep_track_indices.begin(), keep_track_indices.end(),
            [x_temp, y_temp, vertical](int i1, int i2) {
                if (vertical) {
                    return (100000 * x_temp[i1] + y_temp[i1]) < (100000 * x_temp[i2] + y_temp[i2]);
                }
                return (100000 * y_temp[i1] + x_temp[i1]) < (100000 * y_temp[i2] + x_temp[i2]);
            });

        // push in queue
        queue_BFS.push(half_gates_indices(keep_track_indices, gates, 0, gates.size() / 2));
        queue_BFS.push(half_gates_indices(keep_track_indices, gates, gates.size() / 2, gates.size()));

        // Vertical and horizontal alignements follow powers of 2
        if (acc == (v_h_control - 1)) {
            acc = 0;
            v_h_control <<= 1;
            if (vertical) {
                x_div *= 2;
            }
            else {
                x_div *= 2;
                y_div *= 2;
            }
            vertical = (~vertical) & 1;
        }
        else {
            acc++;
        }

        queue_BFS.pop();
    }
    // I think that's it!
}


void unit_test() {
    
    // test up top left function: success
    vector<double> up_left(2, 0);
    cout << " testing up top left function: \n";
    
    for (int i = 0; i < 64; i++) {
        up_left = { 0.0, 0.0 };
        cout << " up_left_corner_coord of acc = " << i << ", v_h_control = 64";
        up_left_corner_coord(up_left, i, 64);
        cout << " (" << up_left[0] << ", " << up_left[1] << ")" << "\n";
    }
    
    // become pad test
    up_left_corner_coord(up_left, 49, 64);
    vector<double> pads_vector_coords_temp;
    vector<int> pad_netlists_temp;
    pads_vector_coords_temp = become_pad(35.0, 70.0, 8, 8, up_left, 49, 1);
    cout << " (35,70) becomes (50,75) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // right

    pads_vector_coords_temp = become_pad(70.0, 70.0, 8, 8, up_left, 49, 1);
    cout << " (70,70) becomes (62.5,70) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // left

    pads_vector_coords_temp = become_pad(0.0, 0.0, 8, 8, up_left, 49, 1);
    cout << " (0,0) becomes (50,62.5) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // up left

    pads_vector_coords_temp = become_pad(100.0, 100.0, 8, 8, up_left, 49, 1);
    cout << " (100,100) becomes (62.5,75) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // down right

    pads_vector_coords_temp = become_pad(60, 80.0, 8, 8, up_left, 49, 1);
    cout << " (60,80) becomes (60,75) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // down

    pads_vector_coords_temp = become_pad(60, 40.0, 8, 8, up_left, 49, 1);
    cout << " (60,40) becomes (60,62.5) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // up

    pads_vector_coords_temp = become_pad(100.0, 0.0, 8, 8, up_left, 49, 1);
    cout << " (100,0) becomes (62.5,62.5) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // up right

    pads_vector_coords_temp = become_pad(0.0, 100.0, 8, 8, up_left, 49, 1);
    cout << " (0,100) becomes (50,75) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // down left

    pads_vector_coords_temp = become_pad(60.0, 70.0, 8, 8, up_left, 49, 1);
    cout << " (60,70) becomes (60,62.5) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // inside

    up_left = { 25.0, 50,0 };
    pads_vector_coords_temp = become_pad(30.0, 70.0, 4, 2, up_left, 3, 0);
    cout << " (30,70) becomes (25,70) :" << " " << pads_vector_coords_temp[0] << " " << pads_vector_coords_temp[1] << "\n"; // inside

}


int main() {
    fstream file_input, file_output;
    int nets, gates, pads;
    vector<vector<int>> conectivity_vector;
    vector<vector<double>> pads_vector_coords;
    vector<int> pad_netlists;

    //unit_test();

    file_input.open("files/primary1");

    init(file_input, gates, nets, pads, conectivity_vector, pads_vector_coords, pad_netlists);

    file_input.close();
    
    valarray<double> x(gates), y(gates);

    place(conectivity_vector, pads_vector_coords, pad_netlists, x, y, gates, 4);

    file_output.open("results/primary1");

    write_result(file_output, x, y);

    file_output.close();

    return 0;
}