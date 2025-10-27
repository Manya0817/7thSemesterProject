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
}; // Removed the typo 'a' here

struct PlacementOption {
    int r, c;
    int level;
    int radius;
    set<pair<int, int>> covered_cells;
};

// --- NEW: Struct for the solution graph ---
struct Edge {
    int to_node; // The index (0 to K-1) of the solution it connects to
    int cost;
};

// --- Global Variables ---
int n, m;
int min_drones_found = INT_MAX;
vector<vector<Placement>> optimal_solutions;
vector<PlacementOption> all_options;

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
 * @brief Step 3: The core backtracking recursive function. (BUG FIXED)
 */
void solve(set<pair<int, int>> cells_to_cover, 
           set<pair<int, int>> available_placement_cells, 
           vector<Placement>& current_solution) {

    // 1. Base Case: All cells are covered.
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

// --- Manhattan distance helper ---
int manhattan_dist(const Placement& a, const Placement& b) {
    return abs(a.r - b.r) + abs(a.c - b.c);
}

/**
 * @brief Solves the assignment problem to find the min cost to transform A to B.
 */
int find_min_cost_valid_match(int drone_A_idx, 
                              const vector<Placement>& sol_A, 
                              const vector<Placement>& sol_B, 
                              vector<bool>& used_B_drones) {
    
    // Base case: All drones in A have been successfully matched.
    if (drone_A_idx == sol_A.size()) {
        return 0; 
    }

    int min_total_cost = INT_MAX;

    // Try to match drone_A_idx with every available drone_B
    for (int j = 0; j < sol_B.size(); ++j) {
        if (!used_B_drones[j]) {
            int move_cost = manhattan_dist(sol_A[drone_A_idx], sol_B[j]);

            // This is the edge condition: move must be at most 1
            if (move_cost <= 1) {
                used_B_drones[j] = true; // "Use" this drone

                int cost_of_rest = find_min_cost_valid_match(drone_A_idx + 1, sol_A, sol_B, used_B_drones);

                if (cost_of_rest != INT_MAX) {
                    min_total_cost = min(min_total_cost, move_cost + cost_of_rest);
                }

                used_B_drones[j] = false; // Backtrack
            }
        }
    }
    
    return min_total_cost; 
}


// --- NEW: CYCLE DETECTION FUNCTIONS ---

/**
 * @brief A recursive DFS helper to detect cycles.
 * @param u The current node.
 * @param p The parent node (to avoid going backward immediately).
 * @param visited The array tracking visited nodes.
 * @param graph The adjacency list of the graph.
 * @return true if a cycle is found, false otherwise.
 */
bool dfs_cycle_check(int u, int p, vector<bool>& visited, const vector<vector<Edge>>& graph) {
    visited[u] = true;

    for (const auto& edge : graph[u]) {
        int v = edge.to_node;

        // If the neighbor is the parent, skip it.
        if (v == p) {
            continue;
        }

        // If the neighbor is already visited (and not the parent), we found a back-edge.
        if (visited[v]) {
            return true; // Cycle detected!
        }

        // Recurse for the unvisited neighbor.
        if (dfs_cycle_check(v, u, visited, graph)) {
            return true; // Cycle found downstream.
        }
    }
    
    // No cycle found from this node.
    return false;
}

/**
 * @brief Checks if the given undirected graph contains a cycle.
 * @param graph The adjacency list of the graph.
 * @return true if a cycle exists, false otherwise.
 */
bool hasCycle(const vector<vector<Edge>>& graph) {
    int K = graph.size();
    if (K == 0) {
        return false;
    }
    
    vector<bool> visited(K, false);
    
    // We must loop through all nodes in case the graph is disconnected.
    for (int i = 0; i < K; ++i) {
        if (!visited[i]) {
            // Start a new DFS. -1 indicates no parent.
            if (dfs_cycle_check(i, -1, visited, graph)) {
                return true;
            }
        }
    }
    
    return false;
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
    
    // --- MODIFIED SECTION ---
    cout << "Setting all cell min/max levels to match drone levels (" 
         << drone.d_level_min << ", " << drone.d_level_max << ")..." << endl;

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            grid[r][c].r = r;
            grid[r][c].c = c;
            grid[r][c].h_level_min = drone.d_level_min;
            grid[r][c].h_level_max = drone.d_level_max;
        }
    }
    // --- END MODIFIED SECTION ---
    
    cout << "---" << endl;

    // 1. Pre-processing
    all_options = generatePlacementOptions(grid, drone);
    cout << "Generated " << all_options.size() << " total possible placement options." << endl;
    cout << "---" << endl;

    // 2. Initialize State
    set<pair<int, int>> cells_to_cover;
    set<pair<int, int>> available_placement_cells;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            cells_to_cover.insert({r, c});
            available_placement_cells.insert({r, c});
        }
    }
    vector<Placement> current_solution;

    // 3. Solve
    cout << "Running solver..." << endl;
    solve(cells_to_cover, available_placement_cells, current_solution);
    cout << "---" << endl;


    // 4. Print Results
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

        // --- 5. Build and Print the Solution Graph ---
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

        for (int i = 0; i < K; ++i) {
            for (int j = i + 1; j < K; ++j) {
                vector<bool> used_B_drones(num_drones, false);
                
                int cost = find_min_cost_valid_match(0, 
                                                     optimal_solutions[i], 
                                                     optimal_solutions[j], 
                                                     used_B_drones);

                if (cost != INT_MAX) {
                    solution_graph[i].push_back({j, cost});
                    solution_graph[j].push_back({i, cost});
                }
            }
        }

        // Print the final graph
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

        // --- 6. NEW: Check for Cycles ---
        cout << "\n=================================================" << endl;
        cout << "Cycle Detection:" << endl;
        
        bool cycleExists = hasCycle(solution_graph);
        
        if (cycleExists) {
            cout << "Result: A cycle EXISTS in the solution graph." << endl;
        } else {
            cout << "Result: NO cycle was found in the solution graph." << endl;
        }
        cout << "=================================================" << endl;
    }

    return 0;
}