#pragma once

#include "utils.hpp"

class Individual
{
    public:
        vector<vector<vector<int>>> g;
        Visit_Table visit_table;
        int nrow, ncol;
        double fit_score;

        Individual() = default;
        Individual(int nrow, int ncol, vector<int>& T, vector<int>& B): nrow(nrow),
                                                                        ncol(ncol),
                                                                        g(nrow, vector<vector<int>>(ncol, vector<int>(2, 0))), 
                                                                        visit_table(nrow, ncol),
                                                                        fit_score(0.0) {
                                                                            for (int i=0; i<ncol; i++)
                                                                            {
                                                                                (*this)(i, 0, 1) = B[i];
                                                                                (*this)(i, nrow-1, 1) = T[i];
                                                                            }
                                                                        };
        Individual(Individual& in, int from, int to): nrow(in.nrow),
                                                    ncol(to-from),
                                                    visit_table(in.nrow, to-from),
                                                    fit_score(0.0) {
                                                        g.resize(in.nrow);
                                                        for (int i=0; i<in.nrow; i++)
                                                        {
                                                            auto it = in.get_row(i).begin();
                                                            g[i].reserve(in.ncol);
                                                            g[i].assign(it+from, it+to);
                                                        }
                                                    };

        int& operator()(int x, int y, int z) { return g[y][x][z]; }
        int predict_extension_row(vector<int>&, vector<int>&, int, bounding_box&);
        int calculate_net_length();
        int calculate_via_number();
        void add_channel(int);
        void add_random_channels(int);
        void delete_channel(int);
        void delete_channels(vector<int>&);
        void delete_channels(bounding_box&);
        void concat_individual(Individual&);
        void remove_if_isolated(vector<int>&);
        vector<int> find_row_index(grid&);
        vector<vector<int>>& get_row(int y) { return g[y]; }
};