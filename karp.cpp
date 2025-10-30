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
vector<vector<Placement>> optimal_solutions; // Raw output from solver
vector<PlacementOption> all_options;

// A very large number to represent infinity
const double INF = 1e18;
// A large number for the Hungarian Algorithm's invalid moves
const int INF_COST = 999999;


// --- CORRECTED HUNGARIAN ALGORITHM (MUNKRES) IMPLEMENTATION ---
class Munkres {
private:
    int n; 
    vector<vector<int>> cost_matrix;
    vector<int> match_x, match_y; 
    vector<int> lx, ly;           
    vector<bool> S, T;           
    vector<int> slack;
    vector<int> slack_x;
    vector<int> prev_path;

    void update_labels() {
        int delta = INT_MAX;
        for (int y = 0; y < n; ++y) {
            if (!T[y]) {
                delta = min(delta, slack[y]);
            }
        }
        for (int x = 0; x < n; ++x) {
            if (S[x]) lx[x] -= delta;
        }
        for (int y = 0; y < n; ++y) {
            if (T[y]) ly[y] += delta;
            else slack[y] -= delta;
        }
    }

    void add_to_tree(int x, int prev_x) {
        S[x] = true;
        prev_path[x] = prev_x;
        for (int y = 0; y < n; ++y) {
            if (lx[x] + ly[y] - cost_matrix[x][y] < slack[y]) {
                slack[y] = lx[x] + ly[y] - cost_matrix[x][y];
                slack_x[y] = x;
            }
        }
    }

public:
    int solve(vector<vector<int>> &input_cost_matrix) {
        n = input_cost_matrix.size();
        if (n == 0) return 0;
        cost_matrix.assign(n, vector<int>(n));
        const int INTERNAL_INF = numeric_limits<int>::max() / 2;
        int max_cost = 0;
        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                if (input_cost_matrix[i][j] != INF_COST && input_cost_matrix[i][j] > max_cost) {
                    max_cost = input_cost_matrix[i][j];
                }
            }
        }
        max_cost += 1; 

        for(int i=0; i<n; ++i) {
            for(int j=0; j<n; ++j) {
                if (input_cost_matrix[i][j] == INF_COST) {
                    cost_matrix[i][j] = -INTERNAL_INF; 
                } else {
                    cost_matrix[i][j] = max_cost - input_cost_matrix[i][j]; 
                }
            }
        }

        lx.assign(n, 0);
        ly.assign(n, 0);
        match_x.assign(n, -1);
        match_y.assign(n, -1);

        for (int x = 0; x < n; ++x) {
            for (int y = 0; y < n; ++y) {
                lx[x] = max(lx[x], cost_matrix[x][y]);
            }
        }

        for(int i=0; i<n; ++i) {
            S.assign(n, false);
            T.assign(n, false);
            prev_path.assign(n, -1);
            vector<int> q(n*2); 
            int wr = 0, rd = 0;

            int root = i;
            if (match_x[root] != -1) continue; 

            q[wr++] = root;
            S[root] = true;

            slack.assign(n, INTERNAL_INF);
            slack_x.assign(n, -1);
            for(int y=0; y<n; ++y) {
                 if (lx[root] + ly[y] - cost_matrix[root][y] < slack[y]) {
                     slack[y] = lx[root] + ly[y] - cost_matrix[root][y];
                     slack_x[y] = root;
                 }
            }

            int x, y;
            bool augmented = false;
            while (!augmented) {
                while (rd < wr && !augmented) {
                    x = q[rd++];
                    for (y = 0; y < n; ++y) {
                        if (cost_matrix[x][y] == lx[x] + ly[y] && !T[y]) {
                            if (match_y[y] == -1) {
                                while (x != -1) {
                                    int prev_y = match_x[x];
                                    match_y[y] = x;
                                    match_x[x] = y;
                                    y = prev_y;
                                    x = prev_path[x];
                                }
                                augmented = true;
                                break; 
                            }
                            T[y] = true;
                            q[wr++] = match_y[y];
                            add_to_tree(match_y[y], x);
                        }
                    }
                }
                if(augmented) break; 

                update_labels(); 
                wr = 0, rd = 0; 
                for (y = 0; y < n; ++y) {
                    if (!T[y] && slack[y] == 0) {
                        if (match_y[y] == -1) {
                            x = slack_x[y];
                            while (x != -1) {
                                int prev_y = match_x[x];
                                match_y[y] = x;
                                match_x[x] = y;
                                y = prev_y;
                                x = prev_path[x];
                            }
                            augmented = true;
                            break; 
                        }
                        T[y] = true;
                        if (!S[match_y[y]]) {
                            q[wr++] = match_y[y];
                            add_to_tree(match_y[y], slack_x[y]);
                        }
                    }
                }
                 if(augmented) break;
            }
        }

        int total_cost = 0;
        bool possible = true;
        for (int x = 0; x < n; ++x) {
            if (match_x[x] == -1 || input_cost_matrix[x][match_x[x]] == INF_COST) {
                possible = false; 
                total_cost = INF_COST;
                break;
            }
            total_cost += input_cost_matrix[x][match_x[x]];
        }

        return total_cost;
    }
};
//---------------------------------------------------------
// --- END OF HUNGARIAN ALGORITHM IMPLEMENTATION ---
//---------------------------------------------------------


/**
 * @brief Step 1: Pre-processes the grid. (Unchanged)
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
                option.r = r; option.c = c; option.level = current_level; option.radius = current_level;
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
    if (current_solution.size() >= min_drones_found) return;
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
        for (const auto& covered_coord : option.covered_cells) next_cells_to_cover.erase(covered_coord);
        next_available_cells.erase({option.r, option.c});
        current_solution.push_back({option.r, option.c, option.level});
        solve(next_cells_to_cover, next_available_cells, current_solution);
        current_solution.pop_back();
    }
}

// --- MODIFIED: Cost function now includes level changes ---
int placement_dist_cost(const Placement& a, const Placement& b) {
    int move_cost = abs(a.r - b.r) + abs(a.c - b.c);
    int level_cost = abs(a.level - b.level);
    return move_cost + level_cost;
}

// --- NEW: HELPER FUNCTIONS FOR DE-DUPLICATION ---

bool comparePlacements(const Placement& a, const Placement& b) {
    if (a.r != b.r) return a.r < b.r;
    if (a.c != b.c) return a.c < b.c;
    return a.level < b.level;
}

vector<vector<Placement>> getUniqueSolutions(const vector<vector<Placement>>& all_solutions) {
    vector<vector<Placement>> unique_solutions;
    
    auto set_comparator = [](const vector<Placement>& a, const vector<Placement>& b) {
        if (a.size() != b.size()) return a.size() < b.size();
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i].r != b[i].r) return a[i].r < b[i].r;
            if (a[i].c != b[i].c) return a[i].c < b[i].c;
            if (a[i].level != b[i].level) return a[i].level < b[i].level;
        }
        return false; // They are equal
    };
    
    set<vector<Placement>, decltype(set_comparator)> seen_signatures(set_comparator);

    for (const auto& sol : all_solutions) {
        vector<Placement> signature = sol;
        sort(signature.begin(), signature.end(), comparePlacements);

        if (seen_signatures.find(signature) == seen_signatures.end()) {
            unique_solutions.push_back(sol); 
            seen_signatures.insert(signature);
        }
    }
    return unique_solutions;
}


// --- KARP'S ALGORITHM (with cycle printing) ---

void findMinMeanWeightCycle(const vector<vector<Edge>>& graph) {
    int K = graph.size(); 
    if (K == 0) {
        cout << "No graph to analyze." << endl;
        return;
    }

    vector<vector<double>> dp(K + 1, vector<double>(K, INF));
    vector<vector<int>> parent(K + 1, vector<int>(K, -1));
    dp[0][0] = 0;

    for (int k = 1; k <= K; ++k) {
        for (int u = 0; u < K; ++u) {
            if (dp[k - 1][u] == INF) continue;
            for (const auto& edge : graph[u]) {
                int v = edge.to_node;
                double cost = (double)edge.cost;
                if (dp[k - 1][u] + cost < dp[k][v]) {
                    dp[k][v] = dp[k - 1][u] + cost;
                    parent[k][v] = u;
                }
            }
        }
    }

    double min_mean_weight = INF;
    int best_node = -1;

    for (int v = 0; v < K; ++v) {
        bool reachable = false;
        for(int k=0; k<=K; ++k) if(dp[k][v] != INF) reachable = true;
        if (!reachable) continue;

        double max_ratio = -INF;

        for (int k = 0; k < K; ++k) {
            if (dp[k][v] != INF) {
                // Using standard mean cost formula: sum / (num_nodes)
                if (K - k > 0) { 
                    max_ratio = max(max_ratio, (dp[K][v] - dp[k][v]) / (K - k));
                }
            }
        }

        if (max_ratio != -INF && max_ratio < min_mean_weight) {
            min_mean_weight = max_ratio;
            best_node = v;
        }
    }

    cout << "\n=================================================" << endl;
    cout << "Karp's Minimum Mean Weight Cycle Algorithm:" << endl;

    if (best_node == -1 || min_mean_weight == INF) {
        cout << "Result: NO cycle was found reachable from Node 1." << endl;
    } else {
        cout << fixed << setprecision(4);
        cout << "Result: The minimum mean weight is: " << min_mean_weight << " (standard mean)" << endl;

        vector<int> path(K + 1);
        vector<int> first_seen(K, -1);
        vector<int> cycle_nodes;

        int curr = best_node;
        for (int k = K; k >= 0; --k) {
            if (curr == -1) break; 
            path[k] = curr;
            curr = parent[k][curr];
        }

        for (int k = 0; k <= K; ++k) {
            int node = path[k];
            if(node == -1) break; 

            if (first_seen[node] != -1) {
                int start_index = first_seen[node];
                for (int j = start_index; j < k; ++j) {
                    cycle_nodes.push_back(path[j]);
                }
                break;
            }
            first_seen[node] = k;
        }

        if (!cycle_nodes.empty()) {
            cout << "   Cycle Found (" << cycle_nodes.size() << " nodes): ";
            for (size_t i = 0; i < cycle_nodes.size(); ++i) {
                cout << (cycle_nodes[i] + 1) << (i == cycle_nodes.size() - 1 ? "" : " -> ");
            }
            cout << " -> " << (cycle_nodes[0] + 1) << endl;
        } else {
             cout << "   Could not reconstruct cycle." << endl;
        }
    }
    cout << "=================================================" << endl;
}


int main() {

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


    if (min_drones_found == INT_MAX) {
        cout << "No solution found." << endl;
    } else {
        cout << "Minimum Drones Required: " << min_drones_found << endl;
        cout << endl;
        
        cout << "Found " << optimal_solutions.size() << " raw solutions. De-duplicating..." << endl;
        vector<vector<Placement>> unique_solutions = getUniqueSolutions(optimal_solutions);
        cout << "Found " << unique_solutions.size() << " Unique Optimal Solution(s):" << endl;

        for (int i = 0; i < unique_solutions.size(); ++i) {
            cout << "  Solution " << (i + 1) << ":" << endl;
            for (const auto& placement : unique_solutions[i]) {
                cout << "    - Drone at (" << placement.r << ", " << placement.c
                     << ") at level " << placement.level << endl;
            }
        }

        cout << "=================================================" << endl;
        cout << "Building Solution Graph..." << endl;
        cout << "=================================================" << endl;
        
        int K = unique_solutions.size();
        if (K == 0) {
            cout << "No solutions to build graph from." << endl;
            return 0;
        }
        vector<vector<Edge>> solution_graph(K);
        int num_drones = unique_solutions[0].size();
        
        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                vector<vector<int>> cost_matrix(num_drones, vector<int>(num_drones));
                for (int drone_i = 0; drone_i < num_drones; ++drone_i) {
                    for (int drone_j = 0; drone_j < num_drones; ++drone_j) {
                        // --- MODIFIED: Calling new cost function ---
                        int dist = placement_dist_cost(unique_solutions[i][drone_i],
                                                       unique_solutions[j][drone_j]);
                        
                        // Check if move is valid (dist <= 1 *per drone*)
                        // This rule is tricky. Our Hungarian implementation finds TOTAL cost.
                        // We must define what a "valid" move is.
                        // Let's stick to your original rule: *each* drone moves at most 1.
                        // We must check this *before* solving.
                        
                        // We'll use the original cost for the matrix
                        int manhattan_move = abs(unique_solutions[i][drone_i].r - unique_solutions[j][drone_j].r) +
                                             abs(unique_solutions[i][drone_i].c - unique_solutions[j][drone_j].c);

                        if (manhattan_move > 1) {
                             cost_matrix[drone_i][drone_j] = INF_COST; // This move is invalid
                        } else {
                             // The move is valid, so its cost is the (r,c) move + level change
                             cost_matrix[drone_i][drone_j] = placement_dist_cost(unique_solutions[i][drone_i],
                                                                                 unique_solutions[j][drone_j]);
                        }
                    }
                }

                Munkres hungarian_solver;
                int total_cost = hungarian_solver.solve(cost_matrix);

                // The total cost is valid if it's less than INF_COST
                // (meaning a valid 1-to-1 matching was found)
                if (total_cost < INF_COST) {
                    solution_graph[i].push_back({j, total_cost});
                    solution_graph[j].push_back({i, total_cost});
                }
            }
        }

        for (int i = 0; i < K; ++i) {
            cout << "\n## Node " << (i + 1) << " (Solution " << (i + 1) << ")" << endl;
            cout << "   Placements:" << endl;
            for (const auto& p : unique_solutions[i]) {
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

        findMinMeanWeightCycle(solution_graph);
    }

    return 0;
}