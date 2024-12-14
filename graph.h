//
// Created by bnkr on 14.12.2024.
//

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>

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


// Возвращает разность двух векторов vec1 и vec2 (элементы vec1, которых нет в vec2)
vector<int> diff_vector(const vector<int>& vec1, const vector<int> &vec2) {
    vector<int> result;
    for (auto x: vec1) {
        bool push_x = true;
        for (auto y: vec2)
            if (x == y) {
                push_x = false;
                break;
            }

        if (push_x)
            result.push_back(x);
    }

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


vector<vector<int> > generate_zero_graph(const int n, const int m) {
    vector<vector<int> > a;
    for (int i = 0; i < n; i++) {
        vector<int> row(m, 0);
        a.push_back(row);
    }
    return a;
}


vector<vector<int> > generate_rand_graph(const int n, const int m) {
    vector<vector<int> > a = generate_zero_graph(n, m);
    srand(clock());
    for (int j = 0; j < m; j++) {
        const int i = rand() % n;
        a[i][j] = 1;
        int k = i;
        while (i == k)
            k = rand() % n;
        a[k][j] = 1;
    }
    return a;
}


void free_matrix(vector<vector<int>> &a) {
    const int size = a.size();
    for (int i = 0; i < size; i++)
        a[i].clear();
}


bool is_contact(const vector<vector<int> > &graph) {
    vector<int> a(graph.size(), 0);
    for (int i = 0; i < graph.size(); i++)
        if (a[i] == i)
            for (int j = 0; j < graph[0].size(); j++)
                if (graph[i][j])
                    for (int k = 0; k < graph.size(); k++)
                        if (k != i && graph[k][j] == 1 && a[k] == 0)
                            a[k] = i + 1;

    for (int i = 0; i < graph.size(); i++) {
        if (a[i] == 0)
            return false;
        i++;
    }

    return true;
}


bool is_connected_vertices(const int i, const int j, const vector<vector<int>> &graph) {
    for (int k = 0; k < graph[0].size(); k++)
        if (graph[i][k] && graph[j][k])
            return true;
    return false;
}


bool find_hamiltonian_cycle(vector<int> a,
                            vector<int> vertics,
                            const int i,
                            const int n,
                            const int m,
                            const vector<vector<int>> &graph) {
    for (int x = 0; x < n; x++) {
        if (is_connected_vertices(a[i-1], x, graph) && x != a[i-1] && vertics[x] == 0) {
            a[i] = x;
            if (a[i-1] == a[0] && i == n)
                return true;

            vertics[x] = 1;
            if (find_hamiltonian_cycle(a, vertics, i+1, n, m, graph))
                return true;
            vertics[x] = 0;
        }
    }

    return false;
}


void get_all_eulerian_graph(const int i,
                            vector<vector<int>> e,
                            const vector<vector<int>> &graph,
                            vector<int> w,
                            vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    const vector<int> y = get_y_from_e(w[i-1], e);
    g = diff_vector(g, y);

    for (const auto x: g) {
        w[i] = x;
        if (x == w[0] && i == e.size()) {
            result.push_back(w);
        } else {
            vector<vector<int>> sub_e = e;
            vector<int> elem = {w[i-1], x};
            sub_e.push_back(elem);
            get_all_eulerian_graph(i+1, sub_e, graph, w, result);
        }
    }
}


void get_all_hahamiltonian_cycle(const int i,
                                 vector<int> v,
                                 const vector<vector<int>> &graph,
                                 vector<int> w,
                                 vector<vector<int>> &result) {
    vector<int> g = get_g(graph, i, w);
    g = diff_vector(g, v);

    for (const auto &x: g) {
        w[i] = x;
        if (x == w[0] && i == v.size()) {
            result.push_back(w);
        } else {
            vector<int> sub_v = v;
            sub_v.push_back(x);
            get_all_hahamiltonian_cycle(i+1, sub_v, graph, w, result);
        }
    }
}

#endif //GRAPH_H