#include<bits/stdc++.h>

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
};a

struct PlacementOption {
    int r, c;
    int level;
    int radius;
    set<pair<int, int>> covered_cells;
};

// --- Global Variables ---
int n, m;
int min_drones_found = INT_MAX;
vector<vector<Placement>> optimal_solutions;
vector<PlacementOption> all_options;

/**
 * @brief Step 1: Pre-processes the grid.
 * (The validation part "max(cell..., drone...)" is here and remains unchanged)
 */
vector<PlacementOption> generatePlacementOptions(const vector<vector<Cell>>& grid, const DroneType& drone) {
    vector<PlacementOption> options;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            const Cell& cell = grid[r][c];
            
            // This validation is still active, as requested
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

    // 1. Base Case: All cells are covered. (MOVED TO THE TOP)
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

    // 2. Pruning: If we've already used more drones than the current best, stop.
    // (MOVED TO SECOND)
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
    // Automatically set all cell levels to match the drone levels
    cout << "Setting all cell min/max levels to match drone levels (" 
         << drone.d_level_min << ", " << drone.d_level_max << ")..." << endl;

    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < m; ++c) {
            grid[r][c].r = r;
            grid[r][c].c = c;
            grid[r][c].h_level_min = drone.d_level_min; // Set to drone min
            grid[r][c].h_level_max = drone.d_level_max; // Set to drone max
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
    }

    return 0;
}