#include<bits/stdc++.h> // Includes all necessary headers

using namespace std;

// --- Data Structures ---
struct Cell {
    int r, c;
    int h_level_min;
    int h_level_max;
};

struct DroneType {
    int d_level_min;
    int d_level_max;
};

struct Placement {
    int r, c;
    int level;
};

struct PlacementOption {
    int r, c;
    int level;
    int radius;
    set<pair<int, int>> covered_cells;
};

struct Edge {
    int to_node;
    int cost;
};

// --- Global Variables ---
int n, m;
int min_drones_found = INT_MAX;
vector<vector<Placement>> optimal_solutions;
vector<PlacementOption> all_options;

// --- NEW: Globals for cycle finding ---
set<vector<int>> unique_cycles;
map<vector<int>, int> cycle_costs;


// --- NEW: HUNGARIAN ALGORITHM (MUNKRES) IMPLEMENTATION ---
// This is a standard, minimal C++ implementation to solve
// the Assignment Problem in O(k^3) time.
// It replaces your O(k!) find_min_cost_valid_match function.
//---------------------------------------------------------

class Munkres {
private:
    vector<vector<int>> matrix;
    vector<int> row_covered, col_covered, path_row, path_col;
    int path_start_row, path_start_col, n_dim;

    void step1(int &step) {
        for (int i = 0; i < n_dim; ++i) {
            int min_val = matrix[i][0];
            for (int j = 1; j < n_dim; ++j) {
                min_val = min(min_val, matrix[i][j]);
            }
            for (int j = 0; j < n_dim; ++j) {
                matrix[i][j] -= min_val;
            }
        }
        step = 2;
    }

    void step2(int &step) {
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[i][j] == 0 && row_covered[j] == 0 && col_covered[i] == 0) {
                    row_covered[j] = 1;
                    col_covered[i] = 1;
                }
            }
        }
        for (int i = 0; i < n_dim; ++i) {
            col_covered[i] = 0;
            row_covered[i] = 0;
        }
        step = 3;
    }

    void step3(int &step) {
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[i][j] == 0) {
                    bool found = false;
                    for (int k = 0; k < n_dim; ++k) {
                        if (matrix[k][j] == 0 && row_covered[k] == 1) {
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        for (int k = 0; k < n_dim; ++k) {
                            if (matrix[i][k] == 0 && col_covered[k] == 1) {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found) {
                        row_covered[i] = 1;
                        col_covered[j] = 1;
                    }
                }
            }
        }
        int col_count = 0;
        for (int j = 0; j < n_dim; ++j) {
            if (col_covered[j] == 1) col_count++;
        }
        if (col_count >= n_dim) step = 7;
        else step = 4;
    }

    void find_zero(int &row, int &col) {
        row = -1; col = -1;
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[i][j] == 0 && row_covered[i] == 0 && col_covered[j] == 0) {
                    row = i; col = j;
                    return;
                }
            }
        }
    }

    int find_star_in_row(int row) {
        for (int j = 0; j < n_dim; ++j) {
            if (row_covered[j] == 1 && matrix[row][j] == 0) return j;
        }
        return -1;
    }

    int find_star_in_col(int col) {
        for (int i = 0; i < n_dim; ++i) {
            if (col_covered[i] == 1 && matrix[i][col] == 0) return i;
        }
        return -1;
    }

    void step4(int &step) {
        int row = -1, col = -1;
        while (true) {
            find_zero(row, col);
            if (row == -1) {
                step = 6;
                return;
            }
            matrix[row][col] = 2; // Prime the zero
            int star_col = find_star_in_row(row);
            if (star_col != -1) {
                row_covered[row] = 1;
                col_covered[star_col] = 0;
            } else {
                path_start_row = row;
                path_start_col = col;
                step = 5;
                return;
            }
        }
    }

    void step5(int &step) {
        vector<pair<int, int>> path;
        path.push_back({path_start_row, path_start_col});
        while (true) {
            int row = find_star_in_col(path.back().second);
            if (row == -1) break;
            path.push_back({row, path.back().second});
            int col = -1;
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[row][j] == 2) {
                    col = j;
                    break;
                }
            }
            path.push_back({row, col});
        }
        for (auto p : path) {
            if (matrix[p.first][p.second] == 1) matrix[p.first][p.second] = 0;
            else if (matrix[p.first][p.second] == 2) matrix[p.first][p.second] = 1;
        }
        for (int i = 0; i < n_dim; ++i) {
            row_covered[i] = 0;
            col_covered[i] = 0;
        }
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[i][j] == 1) {
                    row_covered[j] = 1;
                    col_covered[i] = 1;
                }
            }
        }
        step = 3;
    }

    void step6(int &step) {
        int min_val = INT_MAX;
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (row_covered[i] == 0 && col_covered[j] == 0) {
                    min_val = min(min_val, matrix[i][j]);
                }
            }
        }
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (row_covered[i] == 1) matrix[i][j] += min_val;
                if (col_covered[j] == 0) matrix[i][j] -= min_val;
            }
        }
        step = 4;
    }

public:
    int solve(vector<vector<int>> &cost_matrix, vector<int> &assignment) {
        n_dim = cost_matrix.size();
        matrix = cost_matrix;
        row_covered.assign(n_dim, 0);
        col_covered.assign(n_dim, 0);
        path_row.assign(n_dim * 2, 0);
        path_col.assign(n_dim * 2, 0);

        int step = 1;
        while (step != 7) {
            switch (step) {
                case 1: step1(step); break;
                case 2: step2(step); break;
                case 3: step3(step); break;
                case 4: step4(step); break;
                case 5: step5(step); break;
                case 6: step6(step); break;
            }
        }

        assignment.assign(n_dim, -1);
        int total_cost = 0;
        for (int i = 0; i < n_dim; ++i) {
            for (int j = 0; j < n_dim; ++j) {
                if (matrix[i][j] == 1) {
                    assignment[i] = j;
                    total_cost += cost_matrix[i][j];
                }
            }
        }
        return total_cost;
    }
};
//---------------------------------------------------------
// --- END OF HUNGARIAN ALGORITHM IMPLEMENTATION ---
//---------------------------------------------------------


/**
 * @brief Step 1: Pre-processes the grid.
 */
vector<PlacementOption> generatePlacementOptions(const vector<vector<Cell>>& grid, const DroneType& drone) {
    vector<PlacementOption> options;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            const Cell& cell = grid[r][c];
            
            int valid_min_level = max(cell.h_level_min, drone.d_level_min);
            int valid_max_level = min(cell.h_level_max, drone.d_level_max);

            for (int current_level = valid_min_level; current_level <= valid_max_level; ++current_level) {
                PlacementOption option;
                option.r = r;
                option.c = c;
                option.level = current_level;
                option.radius = current_level;

                for (int rr = 0; rr < n; ++rr) {
                    for (int cc = 0; cc < m; ++cc) {
                        if (abs(r - rr) + abs(c - cc) <= option.radius) {
                            option.covered_cells.insert({rr, cc});
                        }
                    }
                }
                options.push_back(option);
            }
        }
    }
    return options;
}

/**
 * @brief Step 3: The core backtracking recursive function. (Unchanged)
 */
void solve(set<pair<int, int>> cells_to_cover, 
           set<pair<int, int>> available_placement_cells, 
           vector<Placement>& current_solution) {

    // 1. Base Case
    if (cells_to_cover.empty()) {
        int count = current_solution.size();
        if (count < min_drones_found) {
            min_drones_found = count;
            optimal_solutions.clear();
            optimal_solutions.push_back(current_solution);
        } else if (count == min_drones_found) {
            optimal_solutions.push_back(current_solution);
        }
        return;
    }

    // 2. Pruning
    if (current_solution.size() >= min_drones_found) {
        return;
    }

    // 3. Recursive Step
    pair<int, int> target_cell = *cells_to_cover.begin();
    vector<PlacementOption> candidate_options;
    for (const auto& option : all_options) {
        if (option.covered_cells.count(target_cell) && 
            available_placement_cells.count({option.r, option.c})) {
            candidate_options.push_back(option);
        }
    }
    for (const auto& option : candidate_options) {
        set<pair<int, int>> next_cells_to_cover = cells_to_cover;
        set<pair<int, int>> next_available_cells = available_placement_cells;
        for (const auto& covered_coord : option.covered_cells) {
            next_cells_to_cover.erase(covered_coord);
        }
        next_available_cells.erase({option.r, option.c});
        current_solution.push_back({option.r, option.c, option.level});
        solve(next_cells_to_cover, next_available_cells, current_solution);
        current_solution.pop_back();
    }
}

// --- Manhattan distance helper (Unchanged) ---
int manhattan_dist(const Placement& a, const Placement& b) {
    return abs(a.r - b.r) + abs(a.c - b.c);
}

// --- REMOVED: find_min_cost_valid_match ---
// Your old O(k!) recursive function was here and is now deleted.


// --- CYCLE FINDING FUNCTIONS (Unchanged) ---

int get_edge_cost(int u, int v, const vector<vector<Edge>>& graph) {
    for (const auto& edge : graph[u]) {
        if (edge.to_node == v) {
            return edge.cost;
        }
    }
    return 0; 
}

void dfs_find_all_cycles(int u, int p, 
                         const vector<vector<Edge>>& graph,
                         vector<int>& path, 
                         vector<bool>& visited_in_path) {
    
    path.push_back(u);
    visited_in_path[u] = true;

    for (const auto& edge : graph[u]) {
        int v = edge.to_node;
        if (v == p) continue;

        if (visited_in_path[v]) {
            vector<int> cycle_nodes;
            int total_cost = edge.cost; 
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                cycle_nodes.push_back(*it);
                if (*it == v) break; 
            }
            reverse(cycle_nodes.begin(), cycle_nodes.end()); 

            for (size_t i = 0; i < cycle_nodes.size() - 1; ++i) {
                total_cost += get_edge_cost(cycle_nodes[i], cycle_nodes[i+1], graph);
            }

            vector<int> normalized_cycle = cycle_nodes;
            sort(normalized_cycle.begin(), normalized_cycle.end());

            if (unique_cycles.find(normalized_cycle) == unique_cycles.end()) {
                unique_cycles.insert(normalized_cycle);
                cycle_costs[normalized_cycle] = total_cost;
            }
        } else {
            dfs_find_all_cycles(v, u, graph, path, visited_in_path);
        }
    }
    path.pop_back();
    visited_in_path[u] = false;
}

void findAllAndPrintCycles(const vector<vector<Edge>>& graph) {
    int K = graph.size();
    if (K == 0) return;

    unique_cycles.clear();
    cycle_costs.clear();
    vector<bool> visited_in_path(K, false);
    vector<int> path;

    for (int i = 0; i < K; ++i) {
        fill(visited_in_path.begin(), visited_in_path.end(), false);
        path.clear();
        dfs_find_all_cycles(i, -1, graph, path, visited_in_path);
    }

    cout << "\n=================================================" << endl;
    cout << "Cycle Detection Results:" << endl;

    if (unique_cycles.empty()) {
        cout << "Result: NO cycle was found in the solution graph." << endl;
    } else {
        cout << "Result: Found " << unique_cycles.size() << " unique cycle(s)." << endl;
        
        int cycle_num = 1;
        for (const auto& cycle : unique_cycles) {
            cout << "\n  --- Cycle " << cycle_num++ << " ---" << endl;
            cout << "    Nodes (" << cycle.size() << "): ";
            for (size_t i = 0; i < cycle.size(); ++i) {
                cout << (cycle[i] + 1) << (i == cycle.size() - 1 ? "" : " - ");
            }
            cout << endl;
            cout << "    Cost: " << cycle_costs[cycle] << endl;
        }
    }
    cout << "=================================================" << endl;
}


int main() {
    
    // --- Configuration (User Input) ---
    cout << "Enter grid dimensions (n m): ";
    cin >> n >> m;
    DroneType drone;
    cout << "Enter drone's min and max allowed levels (d_min d_max): ";
    cin >> drone.d_level_min >> drone.d_level_max;
    cout << "---" << endl;
    vector<vector<Cell>> grid(n, vector<Cell>(m));
    cout << "Setting all cell min/max levels to match drone levels (" 
         << drone.d_level_min << ", " << drone.d_level_max << ")..." << endl;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            grid[r][c].r = r; grid[r][c].c = c;
            grid[r][c].h_level_min = drone.d_level_min;
            grid[r][c].h_level_max = drone.d_level_max;
        }
    }
    cout << "---" << endl;

    // --- 1, 2, 3: Pre-processing, Init, Solve (Unchanged) ---
    all_options = generatePlacementOptions(grid, drone);
    cout << "Generated " << all_options.size() << " total possible placement options." << endl;
    cout << "---" << endl;
    set<pair<int, int>> cells_to_cover;
    set<pair<int, int>> available_placement_cells;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            cells_to_cover.insert({r, c});
            available_placement_cells.insert({r, c});
        }
    }
    vector<Placement> current_solution;
    cout << "Running solver..." << endl;
    solve(cells_to_cover, available_placement_cells, current_solution);
    cout << "---" << endl;


    // --- 4. Print Results (Unchanged) ---
    if (min_drones_found == INT_MAX) {
        cout << "No solution found." << endl;
    } else {
        cout << "Minimum Drones Required: " << min_drones_found << endl;
        cout << endl;
        cout << "Found " << optimal_solutions.size() << " Optimal Solution(s):" << endl;

        for (int i = 0; i < optimal_solutions.size(); ++i) {
            cout << "  Solution " << (i + 1) << ":" << endl;
            for (const auto& placement : optimal_solutions[i]) {
                cout << "    - Drone at (" << placement.r << ", " << placement.c 
                     << ") at level " << placement.level << endl;
            }
        }

        // --- 5. MODIFIED: Build Solution Graph (Now uses Hungarian) ---
        cout << "=================================================" << endl;
        cout << "Building Solution Graph..." << endl;
        cout << "=================================================" << endl;

        int K = optimal_solutions.size(); 
        if (K == 0) {
            cout << "No solutions to build graph from." << endl;
            return 0;
        }

        vector<vector<Edge>> solution_graph(K);
        int num_drones = optimal_solutions[0].size(); 
        
        // Define a large cost for "invalid" moves (dist > 1)
        const int INF_COST = 999999;

        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                
                // 1. Create the k x k cost matrix
                vector<vector<int>> cost_matrix(num_drones, vector<int>(num_drones));
                for (int drone_i = 0; drone_i < num_drones; ++drone_i) {
                    for (int drone_j = 0; drone_j < num_drones; ++drone_j) {
                        
                        int dist = manhattan_dist(optimal_solutions[i][drone_i], 
                                                  optimal_solutions[j][drone_j]);
                        
                        if (dist <= 1) {
                            cost_matrix[drone_i][drone_j] = dist; // Valid move (cost 0 or 1)
                        } else {
                            cost_matrix[drone_i][drone_j] = INF_COST; // Invalid move
                        }
                    }
                }

                // 2. Solve with Hungarian Algorithm
                Munkres hungarian_solver;
                vector<int> assignment; // This will store the assignment, e.g., assignment[0] = 2
                int total_cost = hungarian_solver.solve(cost_matrix, assignment);

                // 3. Add edge if a valid assignment was found
                if (total_cost < INF_COST) {
                    // A valid mapping (where all moves <= 1) was found
                    solution_graph[i].push_back({j, total_cost});
                    solution_graph[j].push_back({i, total_cost});
                }
            }
        }

        // --- Print the final graph (Unchanged) ---
        for (int i = 0; i < K; ++i) {
            cout << "\n## Node " << (i + 1) << " (Solution " << (i + 1) << ")" << endl;
            cout << "   Placements:" << endl;
            for (const auto& p : optimal_solutions[i]) {
                cout << "     - (" << p.r << "," << p.c << ") L" << p.level << endl;
            }

            cout << "   Edges:" << endl;
            if (solution_graph[i].empty()) {
                cout << "     - None" << endl;
            } else {
                for (const auto& edge : solution_graph[i]) {
                    cout << "     - to Node " << (edge.to_node + 1) 
                         << " (Cost: " << edge.cost << ")" << endl;
                }
            }
        }

        // --- 6. Find and Print All Cycles (Unchanged) ---
        findAllAndPrintCycles(solution_graph);
    }

    return 0;
}