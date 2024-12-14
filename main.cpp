#include <iostream>
#include <vector>
#include <ctime>

#include "graph.h"

using namespace std;


bool is_eulerian_graph(const vector<vector<int>> graph) {
    bool flag = true;

    if (is_contact(graph)) {
        int j = 1;
        int i = 1;
        int count = 0;

        while (i < graph.size() && flag) {
            while (j < graph[0].size() && flag) {
                if (graph[i][j])
                    count++;
                j++;
            }
            if ((count + 1) % 2 == 0)
                flag = false;
            i++;
            j = 0;
        }
    } else {
        flag = false;
    }

    return flag;
}


// bool is_hamiltonian_graph(const vector<vector<int>> graph) {
//     if (is_contact(graph)) {
//         const vector<int> a(graph.size(), 0);
//
//         const vector<int> vertics(graph.size()+1, 0);
//
//         if (find_hamiltonian_cycle(a, vertics, 1, graph.size(), graph[0].size(), graph))
//             return true;
//         return false;
//     }
//
//     return false;
// }


// bool is_hamiltonian_graph(const vector<vector<int>> graph) {
//     vector<vector<int>> result;
//     vector<int>
//
// }


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
            const vector<vector<int>> graph = generate_rand_graph(n, m);
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
    // printf("Table for 8 vertics: \n");
    // gen_graph(8, 28, 1);
    // printf("\n");

    // printf("Table for 9 vertics: \n");
    // gen_graph(9, 36, 1);
    // printf("\n");


    // printf("Table for 10 vertics: \n");
    // gen_graph(10, 45, 2);
    // printf("\n");

    return 0;
}
