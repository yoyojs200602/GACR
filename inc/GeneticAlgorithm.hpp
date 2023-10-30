#pragma once

#include "Individual.hpp"
#include "utils.hpp"


class GeneticAlgorithm
{
    private:
        vector<Individual> Population;
        vector<int> T, B, signals;
        string infile, outfile;
        Individual final_ind;
        unordered_map<int, vector<vector<int>>> signal_map;
        mt19937 gen;
        discrete_distribution<int> selection_gen;

        int Pc;
        int generation;
        int ymin;
        int nsignal;
        int max_descendant;
        int rank;
        bool to_optimize;

    public:
        GeneticAlgorithm(){};
        
        void initialize(Config&);
        void run();
        void parseInfile();
        void initial_population();
        void reduction();
        void optimize(Individual&);
        void mutation();
        void output_ith_best_layout();
        void sort_solution();

        // return true if found solution successfully
        bool crossover(Individual&, Individual&, int);
        bool random_route(Individual&, bounding_box&, vector<int>&, vector<int>&);

        double fitness_calculation();
        Individual& selection();
        Individual mutation(vector<Individual>&);

    private:
        void initialize_visit_table(Individual&, bounding_box&, vector<int>&, int, char);
        bool searchVH(Individual&, bounding_box&, vector<int>&, int, int, char);
        void ensurePath(Individual&, vector<int>&, int, int);
        void concat_individual(Individual&, Individual&);
        void backtrack(Individual&, vector<int>&, int);
        void remove_unused_rows(Individual&);
        void ripup_to_steiner(Individual&, vector<grid>&, int, int, int);
        void ripup_net(Individual&, vector<int>&);
        void plow_net_segment(Individual&, int, int, int, int);

        // return true if found solution successfully
        bool mut_1(Individual&);
        bool mut_2(Individual&);
        bool mut_3(Individual&);
        bool mut_4(Individual&);
        
        vector<int> bfs_search_shortest_path(Individual&, vector<int>&, int, char);
        sig_grid_Map ripup_and_reroute_partial(Individual&, Individual&, int, int);
        sig_grid_Map ripup_rectangle(Individual&, bounding_box&);
};  

