//
// Created by bnkr on 14.12.2024.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include <set>
#include <utility>
#include <algorithm>
#include <iostream>
#include <random>
#include <unordered_set>
#include <utility>

using namespace std;


// Печатает элементы вектора vec через пробел
void print_vec(const vector<int> &vec) {
    for (const int x: vec)
        std::cout << x << " ";
    std::cout << std::endl;
}

// Печатает элементы матрицы а
void print_matrix(const vector<vector<int> > &a) {
    for (const auto &row: a) {
        for (const auto &x: row)
            cout << x << " ";
        cout << endl;
    }
    cout << endl;
}


// Возвращает список вершин, смежных с вершиной w[i-1] в графе graph
vector<int> get_g(const vector<vector<int>> &graph, const int i, const vector<int> &w) {
    const int index = w[i - 1] - 1;
    vector<int> result;
    for (int j = 0; j < graph[index].size(); j++)
        if (graph[index][j])
            result.push_back(j + 1);
    return result;
}


// // Возвращает разность двух векторов vec1 и vec2 (элементы vec1, которых нет в vec2)
// vector<int> diff_vector(const vector<int>& vec1, const vector<int>& vec2) {
//     vector<int> result;
//     for (const auto &x: vec1) {
//         bool flag = true;
//         for (const auto &y: vec2)
//             if (x == y) {
//                 flag = false;
//                 break;
//             }
//         if (flag)
//             result.push_back(x);
//     }
//
//     return result;
// }


vector<int> diff_vector(const vector<int>& vec1, const vector<int>& vec2) {
    unordered_set<int> set2(vec2.begin(), vec2.end());
    vector<int> result;
    for (const auto &x: vec1)
        if (set2.find(x) == set2.end()) // Проверка на отсутствие в set2
            result.push_back(x);
    return result;
}


// Возвращает список вершин, соединенных с вершиной w_i через ребра e
vector<int> get_y_from_e(const int w_i, const vector<vector<int>> &e) {
    vector<int> y;
    for (const auto& elem: e)
        if (elem[0] == w_i)
            y.push_back(elem[1]);
    return y;
}


// Умножает две квадратные матрицы mat1 и mat2
vector<vector<int>> multiplies_square_matrix(const vector<vector<int>> &mat1,
                                             const vector<vector<int>> &mat2) {
    vector<vector<int>> res(mat1.size(), vector<int>(mat1.size(), 0));

    for (int k = 0; k < mat1.size(); k++)
        for (int i = 0; i < mat1.size(); i++)
            for (int j = 0; j < mat1.size(); j++)
                res[i][j] += mat1[i][k] * mat2[k][j];
    return res;
}


// Проверяет, включены ли все элементы vec1 в vec2 (нестрогое включение)
bool not_strikly_include(const vector<int> &vec1, const vector<int> &vec2) {
    for (const auto x: vec1) {
        bool x_not_in_vec2 = true;
        for (const auto y: vec2)
            if (x == y) {
                x_not_in_vec2 = false;
                break;
            }
        if (x_not_in_vec2)
            return false;
    }
    return true;
}


// Проверяет, равны ли два вектора v1 и v2
bool vector_equal(const vector<int> &v1, const vector<int> &v2) {
    for (int i = 0; i < v1.size(); i++)
        if (v1[i] != v2[i])
            return false;

    return true;
}


int count_edges(const vector<vector<int>> &graph) {
    int total_edges = 0;
    for (const auto &row : graph) {
        for (const auto &val : row) {
            total_edges += val;
        }
    }
    return total_edges / 2; // Для неориентированного графа делим на 2
}


// cycle
// Находит все циклы в графе graph
void get_all_cycle(const int i,
                   const vector<vector<int>> &graph,
                   vector<int> w,
                   vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    for (const auto &x: g) {
        w[i] = x;
        if (x == w[0] && i > 2) {
            result.push_back(w);
        } else {
            get_all_cycle(i + 1, graph, w, result);
        }
    }
}


// Находит все простые циклы в графе graph
void get_all_simple_cycle(const int i,
                          vector<int> v,
                          const vector<vector<int>> &graph,
                          vector<int> w,
                          vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    g = diff_vector(g, v);

    for (const auto &x: g) {
        w[i] = x;
        if (x == w[0] && i > 2) {
            result.push_back(w);
        } else {
            vector<int> sub_v = v;
            sub_v.push_back(x);
            get_all_simple_cycle(i+1, sub_v, graph, w, result);
        }
    }
}


vector<vector<int> > generate_zero_graph(const int n) {
    vector<vector<int> > a(n, vector<int>(n, 0));
    return a;
}


vector<vector<int>> generate_rand_graph(const int n, const int m) {
    srand(static_cast<unsigned int>(clock())); // Инициализация генератора случайных чисел
    vector<vector<int>> a = generate_zero_graph(n);

    set<pair<int, int>> edges; // Хранение добавленных рёбер

    while (edges.size() < m) {
        int i = rand() % n;
        int j = rand() % n;

        if (i != j) { // Избегаем петель (i == j)
            int u = min(i, j); // Минимальная вершина
            int v = max(i, j); // Максимальная вершина
            pair<int, int> edge = {u, v}; // Упорядочиваем вершины

            if (edges.find(edge) == edges.end()) { // Проверяем, что такого ребра ещё нет
                edges.insert(edge);
                a[u][v] = 1;
                a[v][u] = 1; // Граф неориентированный
            }
        }
    }

    return a;
}


// Хэш-функция для пары вершин (i, j)
struct PairHash {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
    }
};

// Функция для генерации случайной матрицы смежности
vector<vector<int>> generate_adj_matrix(int n, int m) {
    vector<vector<int>> adj_matrix(n, vector<int>(n, 0));
    unordered_set<pair<int, int>, PairHash> edges;

    random_device rd;
    mt19937 rng(rd());
    uniform_int_distribution<int> dist(0, n - 1);

    while (edges.size() < m) {
        int u = dist(rng);
        int v = dist(rng);
        if (u != v) {
            pair<int, int> edge = minmax(u, v);
            if (edges.insert(edge).second) {
                adj_matrix[edge.first][edge.second] = 1;
                adj_matrix[edge.second][edge.first] = 1;
            }
        }
    }

    return adj_matrix;
}


// Функция для выполнения DFS
void dfs(int node, const vector<vector<int>> &graph, vector<bool> &visited) {
    visited[node] = true;
    for (int neighbor = 0; neighbor < graph[node].size(); ++neighbor) {
        if (graph[node][neighbor] == 1 && !visited[neighbor]) {
            dfs(neighbor, graph, visited);
        }
    }
}

// Функция для проверки связности графа
bool is_contact(const vector<vector<int>> &graph) {
    int n = graph.size(); // Количество вершин в графе
    vector<bool> visited(n, false);

    // Запускаем DFS с первой вершины
    dfs(0, graph, visited);

    // Проверяем, были ли посещены все вершины
    for (bool wasVisited : visited) {
        if (!wasVisited) {
            return false; // Если хоть одна вершина не посещена, граф несвязный
        }
    }
    return true;
}


int degree_vertics(const int v, const vector<vector<int>> &graph) {
    int res = 0;
    for (int j = 0; j < graph.size(); j++)
        res += graph[v-1][j];
    return res;
}



// hamilton
void get_hamiltonian_cycle(const int i,
                           vector<int> v,
                           const vector<vector<int>> &graph,
                           vector<int> w,
                           vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    g = diff_vector(g, v);

    for (const auto &x: g) {
        w.push_back(x);
        if (x == w[0] && i == graph.size()) { //
            result.push_back(w);
        } else {
            v.push_back(x);
            get_hamiltonian_cycle(i+1, v, graph, w, result);
        }
        v.pop_back();
        w.pop_back();
    }
}


vector<vector<int>> get_all_hamiltonian_cycle(const vector<vector<int>> &graph) {
    vector<vector<int>> result;
    for (int i = 1; i <= graph.size(); i++) {
        vector<vector<int>> sub_res;
        vector<int> w = {i};
        vector<int> v;

        get_hamiltonian_cycle(1, v, graph, w, sub_res);

        for (const auto &row: sub_res)
            result.push_back(row);
    }

    return result;
}


bool find_hamiltonian_graph(const int i,
                           vector<int> v,
                           const vector<vector<int>> &graph,
                           vector<int> w) {
    vector<int> g = get_g(graph, i, w);
    g = diff_vector(g, v);

    for (const auto &x: g) {
        w.push_back(x);
        if (x == w[0] && i == graph.size())
            return true;

        v.push_back(x);
        if (find_hamiltonian_graph(i+1, v, graph, w))
            return true;

        v.pop_back();
        w.pop_back();
    }
    return false;
}


bool is_hamiltonian_graph(const vector<vector<int>> &graph) {
    if (!is_contact(graph))
        return false;


    const int n = static_cast<int>(graph.size());
    int amount_true = 0;
    for (int i = 1; i <= graph.size(); i++)
        if (degree_vertics(i, graph) >= n / 2)
            amount_true++;

    if (amount_true == n)
        return true;


    vector<int> w = {1};
    vector<int> v;

    if (find_hamiltonian_graph(1, v, graph, w))
        return true;

    return false;
}




// eulerian
void get_eulerian_cycle(const int i,
                        vector<vector<int>> e,
                        const vector<vector<int>> &graph,
                        vector<int> w,
                        const int amount_edges,
                        vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    vector<int> y = get_y_from_e(w[i - 1], e);
    g = diff_vector(g, y);

    for (const auto &x : g) {
        w.push_back(x);
        if (x == w[0] && i == amount_edges) {
            result.push_back(w);
        } else {
            e.push_back({w[i - 1], x});
            e.push_back({x, w[i - 1]});
            get_eulerian_cycle(i + 1, e, graph, w, amount_edges, result);
        }
        e.pop_back();
        e.pop_back();
        w.pop_back();
    }
}



vector<vector<int>> get_all_eulerian_cycle(const vector<vector<int>> &graph) {
    vector<vector<int>> result;
    for (int i = 1; i <= graph.size(); i++) {
        vector<vector<int>> sub_res;
        vector<int> w = {i};
        vector<vector<int>> e;
        const int amount_edges = count_edges(graph);

        get_eulerian_cycle(1, e, graph, w, amount_edges, sub_res);

        for (const auto &row: sub_res)
            result.push_back(row);
    }

    return result;
}



bool is_eulerian_graph(const vector<vector<int>> &graph) {
    if (!is_contact(graph))
        return false;

    for (int i = 0; i < graph.size(); i++)
        if (degree_vertics(i+1, graph) % 2)
            return false;

    return true;
}


#endif //GRAPH_H