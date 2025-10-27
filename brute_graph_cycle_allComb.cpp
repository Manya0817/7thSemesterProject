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

// --- NEW: Globals for cycle finding ---
set<vector<int>> unique_cycles;
map<vector<int>, int> cycle_costs; // map[normalized_cycle_nodes] -> cost

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


// --- NEW: CYCLE FINDING FUNCTIONS ---

/**
 * @brief Helper function to get the cost of a specific edge (u, v).
 */
int get_edge_cost(int u, int v, const vector<vector<Edge>>& graph) {
    for (const auto& edge : graph[u]) {
        if (edge.to_node == v) {
            return edge.cost;
        }
    }
    return 0; // Should not happen in a consistent graph
}

/**
 * @brief Recursive DFS to find all simple cycles.
 * @param u The current node.
 * @param p The parent node (to avoid going backward immediately).
 * @param graph The adjacency list of the graph.
 * @param path The list of nodes in the current recursion path.
 * @param visited_in_path A bool array to mark nodes in the current path for O(1) lookup.
 */
void dfs_find_all_cycles(int u, int p, 
                         const vector<vector<Edge>>& graph,
                         vector<int>& path, 
                         vector<bool>& visited_in_path) {
    
    // Add current node to path and mark as 'visiting'
    path.push_back(u);
    visited_in_path[u] = true;

    for (const auto& edge : graph[u]) {
        int v = edge.to_node;

        if (v == p) {
            continue; // Skip the parent
        }

        if (visited_in_path[v]) {
            // --- CYCLE DETECTED ---
            // Found a back-edge to an ancestor 'v'
            // The cycle is the path from 'v' back to 'u', plus edge (u,v)

            vector<int> cycle_nodes;
            int total_cost = edge.cost; // Cost of the closing edge (u,v)

            // Iterate backwards from the end of the path (u) until we find v
            for (auto it = path.rbegin(); it != path.rend(); ++it) {
                cycle_nodes.push_back(*it);
                if (*it == v) {
                    break; // Found the start of the cycle
                }
            }
            reverse(cycle_nodes.begin(), cycle_nodes.end()); // Now path is [v, ..., u]

            // Calculate cost of the cycle path (v -> ... -> u)
            for (size_t i = 0; i < cycle_nodes.size() - 1; ++i) {
                total_cost += get_edge_cost(cycle_nodes[i], cycle_nodes[i+1], graph);
            }

            // --- Normalize and Store ---
            vector<int> normalized_cycle = cycle_nodes;
            sort(normalized_cycle.begin(), normalized_cycle.end());

            if (unique_cycles.find(normalized_cycle) == unique_cycles.end()) {
                unique_cycles.insert(normalized_cycle);
                cycle_costs[normalized_cycle] = total_cost;
            }

        } else {
            // --- Recurse ---
            dfs_find_all_cycles(v, u, graph, path, visited_in_path);
        }
    }

    // --- Backtrack ---
    path.pop_back();
    visited_in_path[u] = false;
}

/**
 * @brief Main function to find and print all cycles.
 */
void findAllAndPrintCycles(const vector<vector<Edge>>& graph) {
    int K = graph.size();
    if (K == 0) return;

    // Clear global stores
    unique_cycles.clear();
    cycle_costs.clear();

    vector<bool> visited_in_path(K, false);
    vector<int> path;

    // We must run DFS from every node to find all cycles
    // in all disconnected components
    for (int i = 0; i < K; ++i) {
        // We only start a new DFS if the node isn't already part of a path
        // This is slightly different from standard component-finding
        // We are finding all paths, so we must start from all nodes
        
        // This logic is complex. A simpler, common way:
        // We just need to start the DFS from every node, but avoid re-exploring.
        // `visited_in_path` handles the path logic.
        // We need a *separate* `visited` to handle components.
        
        // Let's use the simplest approach:
        // Run a full DFS from each node 'i' as the root.
        // This will find duplicates, but the `set` will handle it.
        fill(visited_in_path.begin(), visited_in_path.end(), false);
        path.clear();
        dfs_find_all_cycles(i, -1, graph, path, visited_in_path);
    }

    // --- Print All Found Cycles ---
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

        // --- 6. MODIFIED: Find and Print All Cycles ---
        findAllAndPrintCycles(solution_graph);
    }

    return 0;
}