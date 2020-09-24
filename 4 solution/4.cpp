// this is madness // yes this is indeed madness

#include <iostream>
#include <fstream> // file
#include <vector> // vector
#include <sstream> // stringstream
#include <queue> // prio q
#include <math.h> // abs
using namespace std;

class wavefrontCell{
public:
    unsigned int x, y, l;
    int pathcost; // cost of path
    vector<unsigned int> path; // N = 1, S = 2, E = 3, W = 4, U = 5, D = 6

    wavefrontCell(){
        this->x = 0;
        this->y = 0;
        this->l = 0;
        this->pathcost = 0;
        this->path = vector<unsigned int>(1,0);
    }

    wavefrontCell(unsigned int x, unsigned int y, unsigned int layer, int pathcost) {
        this->x = x;
        this->y = y;
        this->l = layer;
        this->pathcost = pathcost;
        this->path = vector<unsigned int>(1, 0);
    }

    wavefrontCell(unsigned int x, unsigned int y, unsigned int layer, int pathcost, vector<unsigned int> path) {
        this->x = x;
        this->y = y;
        this->l = layer;
        this->pathcost = pathcost;
        this->path = path;
    }
};

struct CompareWavefrontCellPathcost {
    bool operator()(wavefrontCell& A, wavefrontCell& B) {
        return A.pathcost > B.pathcost;
    }
};

class Net {
public:
    unsigned int netId;
    int pathcost;
    vector<unsigned int> path;
    unsigned int xs, ys, ls, xt, yt, lt;

    Net(stringstream& s) {
        s >> netId;

        s >> ls;
        ls = ls - 1;
        s >> xs;
        s >> ys;

        s >> lt;
        lt = lt - 1;
        s >> xt;
        s >> yt;

        pathcost = 0;
        path = vector<unsigned int>(0);
    }

    Net(){
        netId = 0;
        xs = 0;
        ys = 0;
        ls = 0;
        xt = 0;
        yt = 0;
        lt = 0;
        pathcost = 0;
        path = vector<unsigned int>(0);
    }
};

void init(fstream& file_grid, fstream& file_nets, unsigned int& bend_penalty, unsigned int& via_penalty, vector<vector<vector<int>>>& grid, vector<Net>& netlist) {

    size_t x_size, y_size;

    file_grid >> x_size;
    file_grid >> y_size;
    file_grid >> bend_penalty;
    file_grid >> via_penalty;

    grid.resize(2, vector<vector<int>>(y_size, vector<int>(x_size))); // layer X y_size X x_size

    stringstream s;
    string line;
    int temp;

    getline(file_grid, line);

    for (size_t layer = 0; layer < 2; layer++) {
        for (size_t y = 0; y < y_size; y++) {
            getline(file_grid, line);
            s.str(line);

            for (size_t x = 0; x < x_size; x++) {
                s >> grid[layer][y][x];
            }

            s.str(std::string()); //clear stringstream // starburst string streamu
            s.clear();

        }
    }

    file_nets >> x_size; // just reusing x_size for the netlist size
    netlist.resize(x_size);
    getline(file_nets, line);

    for (size_t i = 0; i < x_size; i++) {
        getline(file_nets, line);
        s.str(line);
        
        netlist[i] = Net(s);

        s.str(std::string()); //clear stringstream // starburst string streamu
        s.clear();
    }
}

void write_solution(fstream& file_output, vector<Net> netlist) {
    file_output << netlist.size() << endl;
    unsigned int x, y, layer;
    for (size_t i = 0; i < netlist.size(); i++) {
        file_output << netlist[i].netId << endl;

        x = netlist[i].xs;
        y = netlist[i].ys;
        layer = netlist[i].ls;

        file_output << layer + 1 << " " << x << " " << y << endl;

        for (size_t j = 1; j < netlist[i].path.size(); j++) {
            switch (netlist[i].path[j]) { // Start = 0, N = 1, S = 2, E = 3, W = 4, U = 5, D = 6
            case 1:
                y--;
                break;
            case 2:
                y++;
                break;
            case 3:
                x--;
                break;
            case 4:
                x++;
                break;
            case 5:
                file_output << 3 << " " << x << " " << y << endl;
                layer++;
                break;
            case 6:
                file_output << 3 << " " << x << " " << y << endl;
                layer--;
                break;
            }

            file_output << layer + 1 << " " << x << " " << y << endl;
        }

        file_output << 0 << endl;
    }
}

void BacktraceCleanup(wavefrontCell C, vector<vector<vector<int>>>& grid, Net& net) {
    net.path = C.path;
    unsigned int x, y, layer;
    x = net.xs;
    y = net.ys;
    layer = net.ls;

    grid[layer][y][x] = -1;


    for (size_t i = 1; i < net.path.size(); i++) { // the first element is 0, bcs when i call path.back() it throws an error, so this is a little quick patch to go around the error

        switch (net.path[i]) { // Start = 0, N = 1, S = 2, E = 3, W = 4, U = 5, D = 6
        case 1:
            y--;
            break;
        case 2: 
            y++;
            break;
        case 3:
            x--;
            break;
        case 4:
            x++;
            break;
        case 5:
            layer++;
            break;
        case 6:
            layer--;
            break;
        }

        grid[layer][y][x] = -1;  // blocked
    }
}

vector<unsigned int> addElementAtEnd(vector<unsigned int> v, unsigned int elem) {
    v.push_back(elem);
    return v;
}

// I don't like this many if elses, it works but there has not be a better way
void expand(priority_queue<wavefrontCell, vector<wavefrontCell>, CompareWavefrontCellPathcost>& wavefront, wavefrontCell C, vector<vector<vector<int>>> grid, Net net, vector<vector<vector<int>>> pathCostToCell, unsigned int bend_penalty, unsigned int via_penalty) {
    // N = 1, S = 2, E = 3, W = 4, U = 5, D = 6
    unsigned int pathcost;
    if (C.y != 0) {
        if (C.path.back() != 2 && grid[C.l][C.y - 1][C.x] != -1) {

            pathcost = C.pathcost + grid[C.l][C.y - 1][C.x] + bend_penalty * (C.path.back() == 3 || C.path.back() == 4) + (net.xt - C.x) * (net.xt - C.x) + (net.yt - C.y - 1) * (net.yt - C.y - 1);

            if (pathcost < pathCostToCell[C.l][C.y - 1][C.x]) {
                wavefront.emplace(C.x, C.y - 1, C.l, pathcost, addElementAtEnd(C.path, 1));
                pathCostToCell[C.l][C.y - 1][C.x] = pathcost;
            }
        }
    }

    if (C.y + 1 < grid[C.l].size()) {
        if (C.path.back() != 1 && grid[C.l][C.y + 1][C.x] != -1) {

            pathcost = C.pathcost + grid[C.l][C.y + 1][C.x] + bend_penalty * (C.path.back() == 3 || C.path.back() == 4) + (net.xt - C.x) * (net.xt - C.x) + (net.yt - C.y + 1) * (net.yt - C.y + 1);

            if (pathcost < pathCostToCell[C.l][C.y + 1][C.x]) {
                wavefront.emplace(C.x, C.y + 1, C.l, pathcost, addElementAtEnd(C.path, 2));
                pathCostToCell[C.l][C.y + 1][C.x] = pathcost;
            }
        }
    }

    if (C.x != 0) {
        if (C.path.back() != 4 && grid[C.l][C.y][C.x - 1] != -1) {

            pathcost = C.pathcost + grid[C.l][C.y][C.x - 1] + bend_penalty * (C.path.back() == 2 || C.path.back() == 1) + (net.xt - C.x - 1) * (net.xt - C.x - 1) + (net.yt - C.y) * (net.yt - C.y);

            if (pathcost < pathCostToCell[C.l][C.y][C.x - 1]) {
                wavefront.emplace(C.x - 1, C.y, C.l, pathcost, addElementAtEnd(C.path, 3));
                pathCostToCell[C.l][C.y][C.x - 1] = pathcost;
            }
        }
    }

    if (C.x + 1 < grid[C.l][C.y].size()) {
        if (C.path.back() != 3 && grid[C.l][C.y][C.x + 1] != -1) {

            pathcost = C.pathcost + grid[C.l][C.y][C.x + 1] + bend_penalty * (C.path.back() == 2 || C.path.back() == 1) + (net.xt - C.x + 1) * (net.xt - C.x + 1) + (net.yt - C.y) * (net.yt - C.y);

            if (pathcost < pathCostToCell[C.l][C.y][C.x + 1]) {
                wavefront.emplace(C.x + 1, C.y, C.l, pathcost, addElementAtEnd(C.path, 4));
                pathCostToCell[C.l][C.y][C.x + 1] = pathcost;
            }
        }
    }

    if (C.path.back() != 6 && C.l != net.lt && grid[1][C.y][C.x] != -1) {
        
        pathcost = C.pathcost + grid[1][C.y][C.x] + via_penalty + (net.xt - C.x) * (net.xt - C.x) + (net.yt - C.y) * (net.yt - C.y) + 1;
        
        if (pathcost < pathCostToCell[1][C.y][C.x]) {
            wavefront.emplace(C.x, C.y, 1, pathcost, addElementAtEnd(C.path, 5));
            pathCostToCell[1][C.y][C.x] = pathcost;
        }
    }

    if (C.path.back() != 5 && C.l != net.lt && grid[0][C.y][C.x] != -1) {
        
        pathcost = C.pathcost + grid[0][C.y][C.x] + via_penalty + (net.xt - C.x) * (net.xt - C.x) + (net.yt - C.y) * (net.yt - C.y) + 1;
       
        if (pathcost < pathCostToCell[0][C.y][C.x]) {
            wavefront.emplace(C.x, C.y, 0, pathcost, addElementAtEnd(C.path, 6));
            pathCostToCell[0][C.y][C.x] = pathcost;
        }
    }
}

void route(vector<vector<vector<int>>> grid, vector<Net>& netlist, unsigned int bend_penalty, unsigned int via_penalty) {

    priority_queue<wavefrontCell, vector<wavefrontCell>, CompareWavefrontCellPathcost> wavefront;
    wavefrontCell C;
    vector<vector<vector<int>>> pathCostToCell;

    for (size_t i = 0; i < netlist.size(); i++) {
        pathCostToCell = vector<vector<vector<int>>>(2, vector<vector<int>>(grid[0].size(), vector<int>(grid[0][0].size(), 10000)));

        wavefront.emplace(netlist[i].xs, netlist[i].ys, netlist[i].ls, grid[netlist[i].ls][netlist[i].ys][netlist[i].xs]);

        cout << netlist[i].netId << endl;
        while (!wavefront.empty()) {

            C = wavefront.top();
            if (C.x == netlist[i].xt && C.y == netlist[i].yt && C.l == netlist[i].lt) { // success
                BacktraceCleanup(C, grid, netlist[i]);
                break;
            }
            wavefront.pop();
            expand(wavefront, C, grid, netlist[i], pathCostToCell, bend_penalty, via_penalty);
        }
        wavefront = priority_queue<wavefrontCell, vector<wavefrontCell>, CompareWavefrontCellPathcost>(); // reset wavefront for next net
    }
}

int main() {

    fstream file_grid, file_nets, file_output;
    size_t x_size, y_size, net_size;
    unsigned int bend_penalty, via_penalty; 
    vector<vector<vector<int>>> grid;
    vector<Net> netlist;

    std::string file = "fract2";

    file_grid.open("files/" + file + ".grid");
    file_nets.open("files/" + file + ".nl");

    init(file_grid, file_nets, bend_penalty, via_penalty, grid, netlist);

    file_grid.close();
    file_nets.close();

    route(grid, netlist, bend_penalty, via_penalty);

    file_output.open("solution/" + file);

    write_solution(file_output, netlist);

    file_output.close();

    return 0;
}