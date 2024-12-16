#include <iostream>
#include <vector>
#include <ctime>

#include "graph.h"

using namespace std;

void gen_graph(const int n, const int last, const int h) {
    double time_gen = 0.0;
    int amount_eul = 0;
    int amount_ham = 0;
    int amount_union = 0;
    int m = n;

    printf("_________________________________________________________________________\n");
    printf("| Amount vertis | Amount edge | Amount eler | Amount gamilt | All graph |\n");
    while (m <= last) {
        while (time_gen < 0.1) {
            const clock_t start = clock();
            const vector<vector<int> > graph = generate_rand_graph(n, m);
            const clock_t end = clock();
            time_gen += static_cast<double>(end - start) / CLOCKS_PER_SEC;

            if (is_hamiltonian_graph(graph))
                amount_ham++;

            if (is_eulerian_graph(graph))
                amount_eul++;

            amount_union++;
        }
        printf("| %13d | %11d | %11d | %13d | %9d |\n", n, m, amount_eul, amount_ham, amount_union);
        m += h;
        time_gen = 0;
        amount_eul = 0;
        amount_ham = 0;
        amount_union = 0;
    }
}


int main() {
    // cout << "Table for 8 vertics" << endl;
    // gen_graph(8, 28, 1);
    // cout << endl;
    //
    // cout << "Table for 9 vertics" << endl;
    // gen_graph(9, 36, 1);
    // cout << endl;
    //
    // cout << "Table for 10 vertics" << endl;
    // gen_graph(10, 45, 2);
    // cout << endl;


    const vector<vector<int>> a = {
        {0, 1, 1, 0, 0},
        {1, 0, 1, 1, 1},
        {1, 1, 0, 1, 1},
        {0, 1, 1, 0, 0},
        {0, 1, 1, 0, 0},
    };


    const vector<vector<int>> b = {
        {0, 1, 0, 1, 1},
        {1, 0, 1, 0, 0},
        {0, 1, 0, 1, 0},
        {1, 0, 1, 0, 1},
        {1, 0, 0, 1, 0},
    };


    const vector<vector<int>> c = {
        {0, 1, 0, 1},
        {1, 0, 1, 0},
        {0, 1, 0, 1},
        {1, 0, 1, 0},
    };

    cout << "Uele" << endl;
    vector<vector<int>> res = get_all_eulerian_cycle(a);
    print_matrix(res);
    cout << endl;

    cout << "Ham" << endl;
    res = get_all_hamiltonian_cycle(b);
    print_matrix(res);
    cout << endl;

    cout << "Uele and ham" << endl;
    res = get_all_eulerian_cycle(c);
    print_matrix(res);
    cout << endl;
    res = get_all_hamiltonian_cycle(c);
    print_matrix(res);
    cout << endl;

    return 0;
}
