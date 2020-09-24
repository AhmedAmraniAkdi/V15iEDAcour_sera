// BooleanComplement.cpp : This file contains the 'main' function. Program execution begins and ends there.
// TODO : take away total_vars , any list already has size property BUT HEAR ME ON THIS ONE what if the list is emppty, how will u know what the number of vars
// matrix was the better choice, i started taking away variables from the structs, because i didnt need them and all is left is a struct of a vector LUL

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include<algorithm>
#include<math.h>
#include<numeric>
using namespace std;

// maybe continued with the idea of matrix? tried linked list was bad, vector of struct* better, but why no just matrix? like i did in the beggining // lists, deques, other data structures?
struct cubelist {
    cubelist(int n) {
        list.resize(n);
        all_dont_care = 0;
    }
    int all_dont_care;
    vector<int> list{};
};

// works correcyly
// vars go from 0 to n-1
// -1 for 0,
// 1 for 0,
// 0 for don't care
// assignement says won't give all dont care cube.
void init_cubelist(fstream& file_input, vector<struct cubelist*>& R, int &total_vars) {

    int val;
    file_input >> val;
    total_vars = val;
    file_input >> val;
    R.resize(val);

    if (val == 0) {
        return;
    }

    struct cubelist* temp = NULL;

    stringstream s;
    string line;

    // get everyline, then get the value for each array index
    int i = 0;
    while (getline(file_input, line)) {
        if (line.compare("") == 0 || line.compare(" ") == 0) {
            continue;
        }
        temp = new cubelist(total_vars);
        s.str(line);
        s >> val;
        if (val == 0) {
            temp->all_dont_care = 1;
        }
        while (s >> val) {
            if (val < 0) {
                temp->list[val * (-1) - 1] = -1;
            }
            if (val > 0) {
                temp->list[val - 1] = 1;
            }
        }
        R[i] = temp;
        i++;
        s.str(std::string()); //clear stringstream // starburst string streamu
        s.clear();
    }
}

/**************************************************************************************/
/**************************************************************************************/
// writing
void writing_cubelist(fstream& file_output, vector<struct cubelist*>R, int total_vars) {

    file_output << total_vars << "\n";
    file_output << R.size() << "\n";
    
    int non_dont_care; // so this is why you need the number of non dont cares LUL

    for (size_t i = 0; i < R.size(); i++) {
        non_dont_care = 0;
        for (int j = 0; j < total_vars; j++) {
            if (R[i]->list[j] != 0) {
                non_dont_care++;
            }
        }
        file_output << non_dont_care << " ";
        for (int j = 0; j < total_vars; j++) {
            if (R[i]->list[j] == -1) {
                file_output << (-1)*(j+1) << " " ;
            }
            if (R[i]->list[j] == 1) {
                file_output << (j+1) << " ";
            }
        }
        file_output << "\n";
    }
}

// end writing
/**************************************************************************************/
/**************************************************************************************/
/*TERMINATION FUNCTIONS*/ // works correctly


// writes the cubelist to the file with in the assignement format
//void write_output(){}
// deletes all the cubelists from memory then clears the null pointers; R becomes {}
void clear_root(vector<struct cubelist*> &R) {
    for (size_t i = 0; i < R.size(); i++) {
        delete R[i];
    }
    R.clear();
}

// checks whether the cubelist is 1
int check_is_one(vector<struct cubelist*> R) {
    // check = 1 means F is one;
    int check = 0;
    for (size_t i = 0; i < R.size(); i++) {
        if (R[i]->all_dont_care == 1) {
            check = 1;
            return check;
        }
    }
    return check;
}

// negates a cube using morgan laws
void negate_one_cube(vector<struct cubelist*>& R, int total_vars) {
    vector<int> temp_arr = R[0]->list;
    clear_root(R); // will this clear the created structs with objects
    struct cubelist* temp = NULL;
    for (int i = 0; i < total_vars; i++) {
        if (temp_arr[i] != 0) {
            temp = new cubelist(total_vars);
            temp->list[i] = temp_arr[i] * (-1); //negate
            R.push_back(temp);
        }

    }
}

// computes the terminations conditions
int termination(vector<struct cubelist*>& R, int total_vars) {
    // terminarion:
    // empty list means F = 0, return all dont care cube
    if (R.size() == 0) {
        struct cubelist* temp = new cubelist(total_vars);
        temp->all_dont_care = 1;
        R.push_back(temp);
        return 1;
    }
    // is_one measn F = 1, return empty list
    if (check_is_one(R)) {
        clear_root(R);
        return 1;
    }
    // contains just one cube, negate it
    if (R.size() == 1) {
        negate_one_cube(R, total_vars);
        return 1;
    }
    return 0;
}

/*END TERMINATION FUNCTIONS*/
/**************************************************************************************/
/**************************************************************************************/
/*MOST BINATE*/ // works correctly
int splitting_variable(vector<struct cubelist*> R, int total_vars) {
    //debuging
    //cout << "one time I enter\n";
    
    vector<int> true_counts(total_vars, 0);
    vector<int> complement_counts(total_vars, 0);
    vector<int> true_plus_comple_counts(total_vars, 0);
    vector<int> dont_care_counts(total_vars, 0);

    //int binate_count = 0;
    // binate for the indexes of the binate variables
    vector<int> binate;

    for (int i = 0; i < total_vars; i++) {
        for (size_t j = 0; j < R.size(); j++) {
            if (R[j]->list[i] == -1) {
                complement_counts[i]++;
                true_plus_comple_counts[i]++;
                continue;
            }
            if (R[j]->list[i] == 1) {
                true_counts[i]++;
                true_plus_comple_counts[i]++;
                continue;
            }
            if (R[j]->list[i] == 0) {
                dont_care_counts[i]++;
                continue;
            }
        }
        // if either true or complemenet counts are 0 while the other isn't 0 means unate
        if ( !( (true_counts[i] == 0 && complement_counts[i] !=0) || (true_counts[i] != 0 && complement_counts[i] == 0)) 
                && dont_care_counts[i] != R.size()) {
            //binate_count++; // if binate says 0 means all unate
            binate.push_back(i);
        }
    }

    if (binate.size() == 1) {
        return binate[0];
    }

    if (binate.size() != 0) {
        stable_sort(binate.begin(), binate.end(),
            [true_plus_comple_counts](int i1, int i2) {
                return true_plus_comple_counts[i1] > true_plus_comple_counts[i2];
            });
        // if there isnt repetition
        if (true_plus_comple_counts[binate[0]] != true_plus_comple_counts[binate[1]]) {
            return binate[0]; // good
        } // if there is repetition, we see the smallest abs(true - complement)
        else { 
            int i = 1;
            int still_tie = 1;
            int chosen = binate[0];
            while (((size_t)i < binate.size()) && (true_plus_comple_counts[binate[i]] == true_plus_comple_counts[binate[i - 1]])) { //order important
                int a = abs(true_counts[binate[i]] - complement_counts[binate[i]]);
                int b = abs(true_counts[chosen] - complement_counts[chosen]);
                if (a == b && still_tie) {
                    still_tie = 1;
                }
                else {
                    still_tie = 0;
                    if (a < b) {
                        chosen = binate[i];
                    }
                }
                i++;
            }
            if (still_tie) { // if there is still a tie, we return the smallest index 
                return binate[0]; // let's pray it's the first one for NOW ------ WAIT IT'S A STABLE SORTING ALGHORITHM WEEEEEEEEEEEEEEEEEEE // WORKS GOOD
            }
            return chosen; // good
        }
    }
    else { // ALL UNATE
        vector<int> indexes;
        for (size_t i = 0; i < dont_care_counts.size(); i++) {
            if (dont_care_counts[i] != R.size()) {
                indexes.push_back(i);
            }
        }
        stable_sort(indexes.begin(), indexes.end(),
            [true_plus_comple_counts](int i1, int i2) {
                return true_plus_comple_counts[i1] > true_plus_comple_counts[i2];
            });
        return indexes[0]; // stable sorting, first one is the smallest index with the most occurrences on the cubelists // good
    }
}


/*END MOST BINATE*/
/**************************************************************************************/
/**************************************************************************************/
// OR(F1, F2) , concatenates 2 cubelists and returns it

void OR(vector<struct cubelist*> &P, vector<struct cubelist*> N) { // fukcyck struct, having to clear memory manually... MATRIX MY LOVE
    P.reserve(P.size() + N.size());
    P.insert(P.end(), N.begin(), N.end());
    //clear_root(N);
}

// END OR
/**************************************************************************************/
/**************************************************************************************/
// AND(variable, F) , adds the variable with it sign on every cube, goes from don'tcare to the sign
// I didnt expect AND and OR to be this simple compared to the other ones.
void AND(vector<struct cubelist*> &P, int var, int sign) {
    for (size_t i = 0; i < P.size(); i++) {
        P[i]->all_dont_care = 0;
        P[i]->list[var] = sign;
    }
}

// END AND
/**************************************************************************************/
/**************************************************************************************/
// POSITIVE COFACTOR evaluates the variable to 1 or 0, if it goes to 0 then erase the cube, if it one become dont care
// checks for all dont care cube
vector<struct cubelist*> positive_cofactor(vector<struct cubelist*> F, int var) {
    //vector<struct cubelist*> P = F;
    vector<struct cubelist*> P; //pointers bad bad bad bad bad bad bad bad bad -> memory management!!!!!!
    struct cubelist* temp; // see hmed we don't need total_vars matrix cubesxvarsx2 would suffice
    int all_dont_care;
    for (size_t i = 0; i < F.size(); i++) {
        /*if (F[i]->list[var] == -1) { // if the cube becomes 0 simply dont add it
            delete P[i];
            P.erase(P.begin() + i);
            i--;
            continue;
        }*/
        if (F[i]->list[var] == 1) {
            temp = new cubelist(F[0]->list.size());
            temp->list = F[i]->list;
            temp->list[var] = 0;
            all_dont_care = 1;
            // we need to check for all dont care cube
            for (size_t j = 0; j < temp->list.size(); j++) {
                if (temp->list[j] != 0) {
                    all_dont_care = 0;
                    break;
                }
            }
            if (all_dont_care == 1) {
                temp->all_dont_care = 1;
            }
            P.push_back(temp);
        }
        if (F[i]->list[var] == 0) {
            temp = new cubelist(F[0]->list.size());
            temp->list = F[i]->list;
            P.push_back(temp);
        }
    }
    /*
    for (size_t i = 0; i < P.size(); i++) {
        if (P[i]->list[var] == -1) {
            delete P[i];
            P.erase(P.begin() + i);
            i--;
            continue;
        }
        if (P[i]->list[var] == 1) {
            P[i]->list[var] = 0;
            all_dont_care = 1;
            // we need to check for all dont care cube
            for (size_t j = 0; j < P[i]->list.size(); j++) {
                if (P[i]->list[j] != 0) {
                    all_dont_care = 0;
                    break;
                }
            }
            if (all_dont_care == 1) {
                P[i]->all_dont_care = 1;
            }
        }
    }*/
    return P;
}
// end positive cofactor
/**************************************************************************************/
/**************************************************************************************/
// negative cofactor, same as positive cofactor, swapped the ifs
vector<struct cubelist*> negative_cofactor(vector<struct cubelist*> F, int var) {
    vector<struct cubelist*> N;
    struct cubelist* temp; 
    int all_dont_care;
    for (size_t i = 0; i < F.size(); i++) {
        if (F[i]->list[var] == -1) {
            temp = new cubelist(F[0]->list.size());
            temp->list = F[i]->list;
            temp->list[var] = 0;
            all_dont_care = 1;
            // we need to check for all dont care cube
            for (size_t j = 0; j < temp->list.size(); j++) {
                if (temp->list[j] != 0) {
                    all_dont_care = 0;
                    break;
                }
            }
            if (all_dont_care == 1) {
                temp->all_dont_care = 1;
            }
            N.push_back(temp);
        }
        if (F[i]->list[var] == 0) {
            temp = new cubelist(F[0]->list.size());
            temp->list = F[i]->list;
            N.push_back(temp);
        }
    }
    /*vector<struct cubelist*> N = F;
    int all_dont_care;

    for (size_t i = 0; i < N.size(); i++) {
        if (N[i]->list[var] == 1) {
            delete N[i];
            N.erase(N.begin() + i);
            i--;
            continue;
        }
        if (N[i]->list[var] == -1) {
            N[i]->list[var] = 0;
            all_dont_care = 1;
            // we need to check for all dont care cube
            for (size_t j = 0; j < N[i]->list.size(); j++) {
                if (N[i]->list[j] != 0) {
                    all_dont_care = 0;
                    break;
                }
            }
            if (all_dont_care == 1) {
                N[i]->all_dont_care = 1;
            }
        }
    }*/
    return N;
}
// end negative cofactor
/**************************************************************************************/
/**************************************************************************************/
// stuff to make better, pass by reference, pass by value, variables to pass to function
// computes the complement of the function I hope it works LUL
vector<struct cubelist*> complement(vector<struct cubelist*> R, int total_vars) {
    if (termination(R, total_vars)) {
        return R;
    }
    else {
        int x = splitting_variable(R, total_vars);
        //cout << x << "\n";
        vector<struct cubelist*> P = complement(positive_cofactor(R, x), total_vars);
        vector<struct cubelist*> N = complement(negative_cofactor(R, x), total_vars);
        AND(P, x, 1);
        AND(N, x, -1);
        OR(P, N);
        return P;
    }
}


int main()
{
    cout << "testing reading!\n";

    int total_vars;
    vector<struct cubelist*> R;
    fstream file_input, file_output;
    //file_input.open("part1.pcn"); //normal input test
    //file_input.open("part1 - Copie.pcn");// 1 cube test termination :good
    //file_input.open("part1 - Copie - Copie.pcn"); //zero cubes termination :good
    //file_input.open("part1 - Copie - Copie (2).pcn"); //all dont care cube termination :good
    file_input.open("part2.pcn");
    if (!file_input.is_open()) {
        cout << "error";
        return 0;
    }

    init_cubelist(file_input, R, total_vars);
    file_input.close();
    /*
    cout << "printing";
    cout << "\n" << total_vars << " " << R.size();
    for (size_t i = 0; i < R.size(); i++) {
        cout << "\n";
        for (int j = 0; j < total_vars; j++) {
            cout << R[i]->list[j] << " ";
        }
    }*/
    cout << "\n";
    vector<struct cubelist*> result = complement(R, total_vars);

    /*
    cout << "\nafter complementing";
    cout << "\n" << total_vars << " " << result.size();
    for (size_t i = 0; i < result.size(); i++) {
        cout << "\n";
        for (int j = 0; j < total_vars; j++) {
            cout << result[i]->list[j] << " ";
        }
    }*/
    //file_output.open("result_please_work.pcn"); // ma boi works
    file_output.open("part2r.pcn"); // ma boi works
    writing_cubelist(file_output, result, total_vars);
    file_output.close();

    return 0;
}

// das wright, happiness flows through my veins it worked, this abomination worked

// need to work on knowing better the data structures
// nothing gives you more power than debuggin a 15 variable 154 termes boolean function