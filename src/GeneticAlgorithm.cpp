#include "GeneticAlgorithm.hpp"


void GeneticAlgorithm::initialize(Config& cf)
{
    this->infile = cf.infile;
    this->outfile = cf.outfile;
    this->Pc = cf.Pc;
    this->generation = cf.generation;
    this->ymin = cf.ymin;
    this->rank = cf.rank;
    this->nsignal = 0;
    this->to_optimize = cf.to_optimize;

    this->gen = mt19937(random_device{}());
    this->max_descendant = cf.max_descendant;
    
    srand(time(NULL));

    parseInfile();
}


void GeneticAlgorithm::run()
{
    initial_population();
    fitness_calculation();
    reduction();

    for (int i=0; i<generation; i++)
    {
        for (int j=Pc; j<Pc+max_descendant; j++)
        {
            Individual& i0 = selection();
            Individual& i1 = selection();
            if (!crossover(i0, i1, j)) j--;
        }

        fitness_calculation();
        reduction();
        mutation();
    }

    sort_solution();
}


void GeneticAlgorithm::parseInfile()
{
    ifstream fin(infile);
    string line, tmp;

    getline(fin, line);
    stringstream ss(line);
    while (getline(ss, tmp, ' ') && tmp != "\r")
    {
        T.push_back(stoi(tmp));
        nsignal = max(nsignal, T.back());
    }

    getline(fin, line);
    ss.clear();
    ss.str(line);
    while (getline(ss, tmp, ' ') && tmp != "\r")
    {
        B.push_back(stoi(tmp));
        nsignal = max(nsignal, B.back());
    }

    fin.close();

    for (int i=0; i<T.size(); i++)
    {
        this->signal_map[T[i]].push_back({i, 1, 1});
        this->signal_map[B[i]].push_back({i, 0, 1});
    }
    if (this->signal_map.find(0) != this->signal_map.end()) 
        this->signal_map.erase(0);
    for (auto& it : signal_map) 
        this->signals.push_back(it.first);
}


void GeneticAlgorithm::initial_population()
{
    int ncol = T.size();
    auto randInt = uniform_int_distribution<int>(2*ymin, 4*ymin);
    unordered_map<int, vector<vector<int>>> S_, T_;
    vector<vector<int>> sig_;
    sig_.reserve(ncol);
    
    for (int j=0; j<ncol; j++)
    {
        // x, y, z = (j, 0/1, 1)
        if (T[j] != 0)
        {
            S_[T[j]].push_back(vector<int>{j, 1, 1});
            sig_.push_back(S_[T[j]].back());
        }

        if (B[j] != 0)
        {
            S_[B[j]].push_back(vector<int>{j, 0, 1});
            sig_.push_back(S_[B[j]].back());
        }
    }

    Population.resize(this->Pc + this->max_descendant);
    for (int i=0; i<Population.size(); i++)
    {
        int nrow = randInt(gen);
        bool success = true;
        vector<vector<bool>> connected(2, vector<bool>(ncol, false));
        bounding_box bb(-1, ncol, 0, nrow-1);
        shuffle(sig_.begin(), sig_.end(), gen);
        Population[i] = Individual(nrow, ncol, T, B);
        T_.clear();

        for (auto s : sig_)
        {
            if (connected[s[1]][s[0]] == true) continue;
            connected[s[1]][s[0]] = true;

            auto &in = Population[i];
            int sig = s[1] > 0 ? in(s[0], in.nrow-1, 1) : in(s[0], s[1], 1);

            if (T_.find(sig) == T_.end())
            {
                auto t = &(S_[sig][rand() % S_[sig].size()]);
                while (s == *t) { t = &(S_[sig][rand() % S_[sig].size()]); }
                T_[sig].push_back(*t);
                connected[(*t)[1]][(*t)[0]] = true;
            }

            auto t = T_[sig][rand() % T_[sig].size()];
            T_[sig].push_back(s);

            if (s[1] > 0) s[1] = in.nrow-1;
            if (t[1] > 0) t[1] = in.nrow-1;

            if (!random_route(in, bb, s, t))
            {
                success = false;
                break;
            }
        }

        if (success) 
            remove_unused_rows(Population[i]);
        else 
            i--;
    }
}


void GeneticAlgorithm::initialize_visit_table(Individual& in, 
                                            bounding_box& bb,
                                            vector<int>& s, 
                                            int target_in_sig,
                                            char target_visit_sig)
{
    if (in.visit_table(s[0], s[1], 0) != 'x') return;

    queue<point2D> q;

    auto inside_bounding_box = [&](int x, int y) { return x > bb.left && x < bb.right && y > bb.bottom && y < bb.top; };

    if (inside_bounding_box(s[0]-1, s[1]) && in(s[0]-1, s[1], 0) == target_in_sig) q.emplace(s[0]-1, s[1]);
    if (inside_bounding_box(s[0]+1, s[1]) && in(s[0]+1, s[1], 0) == target_in_sig) q.emplace(s[0]+1, s[1]);
    if (inside_bounding_box(s[0], s[1]-1) && in(s[0], s[1]-1, 1) == target_in_sig) q.emplace(s[0], s[1]-1);
    if (inside_bounding_box(s[0], s[1]+1) && in(s[0], s[1]+1, 1) == target_in_sig) q.emplace(s[0], s[1]+1);

    in.visit_table(s[0], s[1], 0) = target_visit_sig;
    in.visit_table(s[0], s[1], 1) = '0';

    while (!q.empty())
    {
        auto& pt = q.front();
        char& cur_direc = in.visit_table(pt.x, pt.y, 1);

        if (cur_direc != '0')
        {
            cur_direc = '0';
            in.visit_table(pt.x, pt.y, 0) = target_visit_sig;
            if (pt.x-1 > bb.left && in(pt.x-1, pt.y, 0) == target_in_sig) q.emplace(pt.x-1, pt.y);
            if (pt.x+1 < bb.right && in(pt.x+1, pt.y, 0) == target_in_sig) q.emplace(pt.x+1, pt.y);
            if (pt.y-1 > bb.bottom && in(pt.x, pt.y-1, 1) == target_in_sig) q.emplace(pt.x, pt.y-1);
            if (pt.y+1 < bb.top && in(pt.x, pt.y+1, 1) == target_in_sig) q.emplace(pt.x, pt.y+1);
        }

        q.pop();
    }
}


bool GeneticAlgorithm::searchVH(Individual& in, 
                                bounding_box& bb,
                                vector<int>& st, 
                                int ori_in_sig,
                                int direc,
                                char target_visit_sig)
{
    // visit_table[i][j][0]: current signal ('s' or 't')
    // visit_table[i][j][1]: direction of pervious grid (right: 1, top: 2, left: 3, bottom: 4)
    char tmp_sig = in.visit_table(st[0], st[1], 0);
    if (tmp_sig == target_visit_sig || bb.top == bb.bottom+1 || bb.right == bb.left+1) return true;

    int tmp_in_sig = in(st[0], st[1], direc);
    if (tmp_in_sig != 0 && tmp_in_sig != ori_in_sig) return false;

    int max_ = bb.top, min_ = bb.bottom, ori_start_grid = st[direc];
    if (0 == direc)
    {
        if (st[1] == bb.bottom || st[1] == bb.top)
            return false;
        else
        {
            max_ = bb.right;
            min_ = bb.left;
        }
    }

    // check grid toward up/right
    while (++st[direc] < max_)
    {
        char& sig = in.visit_table(st[0], st[1], 0);
        tmp_in_sig = in(st[0], st[1], direc);

        if (tmp_in_sig != 0 && tmp_in_sig != ori_in_sig)
        {
            max_ = st[direc];
            break;
        }

        if (sig == target_visit_sig) return true;
        if (sig == 'x') sig = tmp_sig;
    }

    // check grid toward bottom/left
    st[direc] = ori_start_grid;
    while (--st[direc] > min_)
    {
        char& sig = in.visit_table(st[0], st[1], 0);
        tmp_in_sig = in(st[0], st[1], direc);

        if (tmp_in_sig != 0 && tmp_in_sig != ori_in_sig)
        {
            min_ = st[direc];
            break;
        }

        if (sig == target_visit_sig)  return true;
        if (sig == 'x') sig = tmp_sig;
    }

    st[direc] = max_ > min_+2 ? rand() % (max_ - min_ - 1) + min_ + 1 : ori_start_grid;

    return false;
}


void GeneticAlgorithm::ensurePath(Individual& in, 
                                vector<int>& st,
                                int direc,
                                int ori_in_sig)
{
    char tmp_visit_sig = in.visit_table(st[0], st[1], 0);
    char target_visit_sig = tmp_visit_sig == 's' ? 't' : 's';

    auto st0 = bfs_search_shortest_path(in, st, ori_in_sig, tmp_visit_sig);
    auto st1 = bfs_search_shortest_path(in, st, ori_in_sig, target_visit_sig);
    in.visit_table(st[0], st[1], 1) = '0';

    backtrack(in, st0, ori_in_sig);
    backtrack(in, st1, ori_in_sig);
}


vector<int> GeneticAlgorithm::bfs_search_shortest_path(Individual& in, 
                                                        vector<int>& st,
                                                        int ori_in_sig, 
                                                        char target_visit_sig)
{
    if (in.visit_table(st[0], st[1], 0) == target_visit_sig 
        && in.visit_table(st[0], st[1], 1) == '0') return st;

    static const int directions[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
    static const int direc[4] = {0, 0, 1, 1};
    static const char symbols[4] = {'1', '3', '2', '4'};
    int max_row = in.nrow - 1, max_col = in.ncol - 1;
    queue<pair<int, int>> q;
    q.emplace(st[0], st[1]);

    while (!q.empty())
    {
        auto& pt = q.front();

        for (int i=0; i<4; i++)
        {
            int new_col = pt.first + directions[i][0];
            int new_row = pt.second + directions[i][1];

            if (new_row >= 0 && new_row <= max_row && new_col >= 0 && new_col <= max_col &&
                in.visit_table(new_col, new_row, 0) == target_visit_sig)
            {
                char& visit_direc = in.visit_table(new_col, new_row, 1);
                int in_direc = in(new_col, new_row, direc[i]);

                if (visit_direc == '0')
                {
                    visit_direc = symbols[i];
                    in.visit_table(st[0], st[1], 1) = '0';
                    return { new_col, new_row };
                }

                if (visit_direc == 'x' && (in_direc == 0 || in_direc == ori_in_sig)) 
                {
                    visit_direc = symbols[i];
                    q.emplace(new_col, new_row);
                }
            }
        }

        q.pop();
    }

    return st;
}


void GeneticAlgorithm::backtrack(Individual& in, 
                                vector<int>& st,
                                int in_target_sig)
{
    while (in.visit_table(st[0], st[1], 1) != '0')
    {
        if (in.visit_table(st[0], st[1], 1) == '1')
        {
            while (in.visit_table(st[0], st[1], 1) == '1') { in(st[0]++, st[1], 0) = in_target_sig; }
            in(st[0], st[1], 0) = in_target_sig;
        }
        if (in.visit_table(st[0], st[1], 1) == '2')
        {
            while (in.visit_table(st[0], st[1], 1) == '2') { in(st[0], st[1]++, 1) = in_target_sig; }
            in(st[0], st[1], 1) = in_target_sig;
        }
        if (in.visit_table(st[0], st[1], 1) == '3')
        {
            while (in.visit_table(st[0], st[1], 1) == '3') { in(st[0]--, st[1], 0) = in_target_sig; }
            in(st[0], st[1], 0) = in_target_sig;
        }
        if (in.visit_table(st[0], st[1], 1) == '4')
        {
            while (in.visit_table(st[0], st[1], 1) == '4') { in(st[0], st[1]--, 1) = in_target_sig; }
            in(st[0], st[1], 1) = in_target_sig;
        }
    }
}


bool GeneticAlgorithm::random_route(Individual& in, 
                                    bounding_box& bb, 
                                    vector<int>& s, 
                                    vector<int>& t)
{
    int dirs = s[2], dirt = t[2];
    int itr = 3 * abs(s[0] - t[0]) + in.nrow + 10;
    int target_sig = in(s[0], s[1], s[2]);
    int extension_cnt = 0;

    while(extension_cnt++ < 10)
    {
        in.visit_table.clear();
        initialize_visit_table(in, bb, s, target_sig, 's');
        initialize_visit_table(in, bb, t, target_sig, 't');
        auto ss = s, st = t;
        int cnt = itr;

        while (cnt > 0)
        {
            if (searchVH(in, bb, ss, target_sig, dirs, 't'))
            {
                ensurePath(in, ss, dirs, target_sig);
                break;
            }

            if (searchVH(in, bb, st, target_sig, dirt, 's'))
            {
                ensurePath(in, st, dirt, target_sig);
                break;
            }

            dirs = 1 - dirs;
            dirt = 1 - dirt;
            cnt--;
        }

        if (cnt == 0)
        {
            int new_row_idx = in.predict_extension_row(s, t, target_sig, bb);
            in.add_channel(new_row_idx);
            if (s[1] >= new_row_idx) s[1]++;
            if (t[1] >= new_row_idx) t[1]++;
            bb.top++;
        }
        else
            return true;
    }

    return false;
}


Individual& GeneticAlgorithm::selection()
{
    return Population[this->selection_gen(this->gen)];
}


void GeneticAlgorithm::ripup_to_steiner(Individual& in, 
                                        vector<grid>& existed_pt,
                                        int x, 
                                        int y, 
                                        int direc)
{
    int target_sig = in(x, y, direc);
    int max_row = in.nrow-1, max_col = in.ncol-1;

    while (true)
    {
        if (y == 0 || y == max_row) { 
            existed_pt.push_back({x, &in(x, y, direc), 1}); break; }

        if (in(x, y, 0) == target_sig && in(x, y, 1) == target_sig)
        {
            int cnt[2] = {(x > 0 && in(x-1, y, 0) == target_sig) + (x < max_col && in(x+1, y, 0) == target_sig),
                        (y > 0 && in(x, y-1, 1) == target_sig) + (y < max_row && in(x, y+1, 1) == target_sig)};

            existed_pt.erase(std::remove_if(existed_pt.begin(), existed_pt.end(), 
                            [&](const grid& pt) { return get<1>(pt) == &in(x, y, 0) || get<1>(pt) == &in(x, y, 1); }), 
                            existed_pt.end());

            if (cnt[0]+cnt[1] > 1)
            {
                existed_pt.push_back({x, &in(x, y, direc), direc});
                break;
            }
            else
            {
                in(x, y, direc) = in(x, y, 1-direc) = 0;
                if (cnt[direc] == 0) {
                    direc = 1 - direc;
                }
            }
        }
        else
            in(x, y, direc) = 0;

        if (direc == 0)
        {
            if (x > 0 && in(x-1, y, 0) == target_sig) x--;
            else if (x < max_col && in(x+1, y, 0) == target_sig) x++;
            else break;
        }
        else
        {
            if (y > 0 && in(x, y-1, 1) == target_sig) y--;
            else if (y < max_row && in(x, y+1, 1) == target_sig) y++;
            else break;
        }
    }
}


void GeneticAlgorithm::ripup_net(Individual& in, vector<int>& st)
{
    int ori_in_sig = in(st[0], st[1], st[2]);
    int max_row = in.nrow-1, max_col = in.ncol-1;
    queue<point2D> q;
    q.emplace(st[0], st[1]);

    while (!q.empty())
    {
        auto &pt = q.front();
        int &hor_ = in(pt.x, pt.y, 0), &ver_ = in(pt.x, pt.y, 1);

        if (hor_ == ori_in_sig)
        {
            if (pt.x > 0 && in(pt.x-1, pt.y, 0) == ori_in_sig) q.emplace(pt.x-1, pt.y);
            if (pt.x < max_col && in(pt.x+1, pt.y, 0) == ori_in_sig) q.emplace(pt.x+1, pt.y);
            hor_ = 0;
        }

        if (ver_ == ori_in_sig)
        {
            if (pt.y > 0 && in(pt.x, pt.y-1, 1) == ori_in_sig) q.emplace(pt.x, pt.y-1);
            if (pt.y < max_row && in(pt.x, pt.y+1, 1) == ori_in_sig) q.emplace(pt.x, pt.y+1);
            if (pt.y != 0 && pt.y != max_row) ver_ = 0;
        }
    
        q.pop();
    }
}


void GeneticAlgorithm::remove_unused_rows(Individual& in)
{
    vector<int> remove_list;
    remove_list.reserve(in.nrow);
    for (int i=1; i<in.nrow-1; i++)
        if (!any_of(in.get_row(i).begin(), in.get_row(i).end(), [&](vector<int>& grid){ return grid[0] != 0; }))
            remove_list.push_back(i);
    in.delete_channels(remove_list);
}



sig_grid_Map GeneticAlgorithm::ripup_and_reroute_partial(Individual& in, 
                                                        Individual& in_ori, 
                                                        int cut_column, 
                                                        int ripup_col)
{
    sig_grid_Map mp;
    vector<int> signal_to_reroute;
    vector<int> remove_list;
    signal_to_reroute.reserve(in.nrow);

    // rip up nets crossing the cut column
    for (int i=1; i<in.nrow-1; i++)
        if (in(ripup_col, i, 0) != 0 && in(ripup_col, i, 0) == in_ori(cut_column, i, 0))
            ripup_to_steiner(in, mp[in(ripup_col, i, 0)], ripup_col, i, 0);
    for (auto& sig_pt : mp)
        if (!sig_pt.second.empty())
            signal_to_reroute.push_back(sig_pt.first); 
    
    // scan and remove unused rows
    for (int i=1; i<in.nrow-1; i++)
    {
        if (any_of(in.get_row(i).begin(), in.get_row(i).end(), [&](vector<int>& grid){ return grid[0] != 0; }))
            continue;
        remove_list.push_back(i);
    }
    in.delete_channels(remove_list);

    // decide the order of reroute
    shuffle(signal_to_reroute.begin(), signal_to_reroute.end(), this->gen);

    // reroute the (steiner) points of each ripped net
    bounding_box bb(-1, in.ncol, 0, in.nrow-1);
    for (auto sig : signal_to_reroute)
    {
        auto& steiner_pts = mp[sig];
        auto pt0 = in.find_row_index(steiner_pts.front());

        for (int i=steiner_pts.size()-1; i>0; i--)
        {
            auto pt1 = in.find_row_index(steiner_pts.back());

            if (!random_route(in, bb, pt0, pt1))
            {
                auto failed_map = sig_grid_Map();
                failed_map[-1] = {};
                return failed_map;
            }

            in.remove_if_isolated(pt1);
            steiner_pts.pop_back();
        }
    }

    // return map of < signal, vector of the (steiner) point remained >
    return mp;
}


// use i0, i1 as parent, do crossover to generate a child at Population[idx]
bool GeneticAlgorithm::crossover(Individual& i0, Individual& i1, int idx)
{
    int cut_column = 1 + (rand()%(i0.ncol-3));

    Population[idx] = Individual(i0, 0, cut_column);
    auto &in_l = Population[idx];
    Individual in_r(i1, cut_column, i1.ncol);

    auto mp_l = ripup_and_reroute_partial(in_l, i0, cut_column, cut_column-1);
    if (mp_l.find(-1) != mp_l.end()) return false;

    auto mp_r = ripup_and_reroute_partial(in_r, i1, cut_column, 0);
    if (mp_r.find(-1) != mp_r.end()) return false;

    concat_individual(in_l, in_r);

    bounding_box bb(-1, in_l.ncol, 0, in_l.nrow-1);

    for (auto& sig : mp_r)
    {
        if (sig.second.empty()) continue;
        auto pt = in_r.find_row_index(sig.second.front());
        pt[0] += cut_column;
        mp_l[sig.first].emplace_back(pt[0], &in_l(pt[0], pt[1], pt[2]), pt[2]);
    }

    for (auto& sig : mp_l)
    {
        if (sig.second.size() < 2) continue;
        auto pt0 = in_l.find_row_index(sig.second.front());
        auto pt1 = in_l.find_row_index(sig.second.back());
        if (!random_route(in_l, bb, pt0, pt1)) return false;

        in_l.remove_if_isolated(pt0);
        in_l.remove_if_isolated(pt1);
    }

    remove_unused_rows(in_l);

    return true;
}


void GeneticAlgorithm::reduction()
{
    sort(Population.begin(), Population.end(), [&](Individual& i0, Individual& i1){
        return i0.fit_score > i1.fit_score;
    });
}


void GeneticAlgorithm::optimize(Individual& in)
{
    Individual tmp = in;

    mut_1(tmp);
    tmp.nrow < in.nrow ? in = tmp : tmp = in;

    mut_2(tmp);
    tmp.nrow < in.nrow ? in = tmp : tmp = in;

    mut_3(tmp);
    tmp.nrow < in.nrow ? in = tmp : tmp = in;

    mut_4(tmp);
    if (tmp.nrow < in.nrow) in = tmp;
}


sig_grid_Map GeneticAlgorithm::ripup_rectangle(Individual& in, bounding_box& bb)
{
    sig_grid_Map mp;

    // save top/bottom row signal to reroute map
    for (int i=bb.left+1; i<bb.right; i++)
    {
        int *top = &in(i, bb.top, 1), *bottom = &in(i, bb.bottom, 1);
        if (*top != 0 && in(i, bb.top-1, 1) == *top) mp[*top].emplace_back(i, top, 1);
        if (*bottom != 0 && in(i, bb.bottom+1, 1) == *bottom) mp[*bottom].emplace_back(i, bottom, 1);
    }

    // save right/left colummn signal to reroute map
    for (int i=bb.bottom+1; i<bb.top; i++)
    {
        int *left = &in(bb.left, i, 0), *right = &in(bb.right, i, 0);
        if (*right != 0 && in(bb.right-1, i, 0) == *right) mp[*right].emplace_back(bb.right, right, 0);
        if (*left != 0 && in(bb.left+1, i, 0) == *left) mp[*left].emplace_back(bb.left, left, 0);
    }

    in.delete_channels(bb);

    return mp;
}


void GeneticAlgorithm::mutation()
{
    uniform_real_distribution<double> dis(0.0, 1.0);

    for (int i=0; i<Pc; i++)
    {
        if (dis(this->gen) < 0.01) while(!mut_1(Population[i])){};
        if (dis(this->gen) < 0.01) while(!mut_2(Population[i])){};
        if (dis(this->gen) < 0.01) mut_3(Population[i]);
        if (dis(this->gen) < 0.05) while(!mut_4(Population[i])){};
    }
}


bool GeneticAlgorithm::mut_1(Individual& in)
{
    int x0 = rand() % in.ncol, x1 = rand() % in.ncol;
    int y0 = rand() % (in.nrow-2)+1, y1 = rand() % (in.nrow-2)+1;
    while (x0 == x1) { x1 = rand() % in.ncol; }
    while (y0 == y1) { y1 = rand() % (in.nrow-2)+1; }
    if (x0 > x1) swap(x0, x1);
    if (y0 > y1) swap(y0, y1);

    Individual backup = in;
    bounding_box bb(x0, x1, y0, y1);
    auto mp = ripup_rectangle(in, bb);

    for (auto& sig_pt : mp)
    {
        auto& pts = sig_pt.second;
        auto pt0 = in.find_row_index(pts.front());

        for (int i=pts.size()-1; i>0; i--)
        {
            auto pt1 = in.find_row_index(pts.back());

            if (!random_route(in, bb, pt0, pt1))
            {
                in = backup;
                return false;
            }

            pts.pop_back();
        }
    }

    return true;
}


bool GeneticAlgorithm::mut_2(Individual& in)
{
    int ripup_cnt = 1 + random() % (this->signals.size()-1);
    Individual backup = in;
    shuffle(signals.begin(), signals.end(), this->gen);

    for (int i=0; i<ripup_cnt; i++) 
    {
        auto st = signal_map[signals[i]].front();
        if (st[1] > 0) st[1] = in.nrow-1; 
        ripup_net(in, st);
    }

    remove_unused_rows(in);

    bounding_box bb(-1, in.ncol, 0, in.nrow-1);

    for (int i=0; i<ripup_cnt; i++)
    {
        auto& signal_pts = this->signal_map[this->signals[i]];
        auto s = signal_pts.front();

        for (int j=1; j<signal_pts.size(); j++)
        {
            auto t = signal_pts[j];
            if (s[1] > 0) s[1] = in.nrow-1;
            if (t[1] > 0) t[1] = in.nrow-1;

            if (!random_route(in, bb, s, t))
            {
                in = backup;
                return false;
            }
        }
    }

    return true;
}


bool GeneticAlgorithm::mut_3(Individual& in)
{
    int random_row = 1 + random() % (in.nrow-2), target_row;

    in.add_channel(random_row);

    for (int i=0; i<in.ncol; i++)
    {
        int above_ver = in(i, random_row+1, 1), below_ver = in(i, random_row-1, 1);
        int above_hor = in(i, random_row+1, 0), below_hor = in(i, random_row-1, 0);
        bool above = above_hor != 0 && !(above_ver != 0 && above_ver != above_hor);
        bool below = below_hor != 0 && !(below_ver != 0 && below_ver != below_hor);

        if (!above && !below) continue;
        if (above && below) target_row = random() % 2 == 0 ? random_row-1 : random_row+1;
        else target_row = below ? random_row-1 : random_row+1;

        int from = i++, target_in_sig = in(from, target_row, 0);
        while (i < in.ncol && in(i, target_row, 0) == target_in_sig) { i++; }

        plow_net_segment(in, from, --i, target_row, random_row);
    }

    return true;
}


bool GeneticAlgorithm::mut_4(Individual& in)
{
    int random_row = 1 + rand() % (in.nrow-2);
    auto& row = in.get_row(random_row);
    Individual backup = in;
    sig_grid_Map mp;

    for (int i=0; i<in.ncol; i++)
    {
        int& ver_ = in(i, random_row, 1);

        if (ver_ != 0 && in(i, random_row, 0) == ver_)
        {
            int tmp_in_sig = ver_;
            ver_ = 0;
            if (in(i, random_row-1, 1) == tmp_in_sig) ripup_to_steiner(in, mp[tmp_in_sig], i, random_row-1, 1);
            if (in(i, random_row+1, 1) == tmp_in_sig) ripup_to_steiner(in, mp[tmp_in_sig], i, random_row+1, 1);
        }
    }

    in.delete_channel(random_row);

    bounding_box bb(-1, in.ncol, 0, in.nrow-1);

    for (auto& sig_pt : mp)
    {
        auto& pts = sig_pt.second;
        int pivot = 0;
        vector<int> pt0 = in.find_row_index(pts[pivot]);
        while (pt0.empty()) { pt0 = in.find_row_index(pts[++pivot]); }

        for (int i=pts.size()-1; i>pivot; i--)
        {
            auto pt1 = in.find_row_index(pts.back());

            if (!pt1.empty() && !random_route(in, bb, pt0, pt1))
            {
                in = backup;
                return false;
            }

            pts.pop_back();
        }
    }

    remove_unused_rows(in);

    return true;
}


void GeneticAlgorithm::plow_net_segment(Individual& in, 
                                        int x0, 
                                        int x1, 
                                        int row_from, 
                                        int row_to)
{
    if (x0 >= x1) return;

    int target_in_sig = in(x0, row_from, 0);
    int check_connect_row_from, check_connect_row_to;
    auto &rf = in.get_row(row_from), &rt = in.get_row(row_to);

    if (row_to > row_from)
    {
        check_connect_row_from = row_from-1;
        check_connect_row_to = row_to+1;
    } 
    else
    {
        check_connect_row_from = row_from+1;
        check_connect_row_to = row_to-1;
    }

    for (int i=x0; i<=x1; i++)
    {
        rf[i][0] = 0;
        rt[i][0] = target_in_sig;
        rt[i][1] = rf[i][1];
        if (in(i, check_connect_row_from, 1) != target_in_sig && 
            in(i, check_connect_row_to, 1) == target_in_sig) rf[i][1] = 0;
    }

    if (x0 > 0 && in(x0-1, row_from, 0) == target_in_sig) rf[x0][0] = rf[x0][1] = rt[x0][1] = target_in_sig;
    if (x1 < in.ncol-1 && in(x1+1, row_from, 0) == target_in_sig) rf[x1][0] = rf[x1][1] = rt[x1][1] = target_in_sig;
    
}


void GeneticAlgorithm::concat_individual(Individual& in_l, Individual& in_r)
{
    if (in_l.nrow > in_r.nrow) in_r.add_random_channels(in_l.nrow - in_r.nrow);
    else if (in_r.nrow > in_l.nrow) in_l.add_random_channels(in_r.nrow - in_l.nrow);

    for (int i=0; i<in_l.nrow; i++)
        in_l.get_row(i).insert(in_l.get_row(i).end(), in_r.get_row(i).begin(), in_r.get_row(i).end());
    
    in_l.ncol += in_r.ncol;
    in_l.visit_table.addCols(in_r.ncol);
}


double GeneticAlgorithm::fitness_calculation()
{
    vector<double> F1, F2;
    vector<int> idx2in;
    int size = Population.size();
    F1.reserve(size);
    F2.reserve(size);
    idx2in.reserve(size);

    for (int i=0; i<size; i++)
    {
        F1.push_back(1.0/Population[i].nrow);
        F2.push_back(1.0/Population[i].calculate_net_length());
        idx2in.push_back(i);
    }

    sort(idx2in.begin(), idx2in.end(), [&](int a, int b){
        if (F1[a] == F1[b])
            return F2[a] > F2[b];
        return F1[a] > F1[b];
    });

    for (int i=0; i<size; i++)
    {
        int from = i++;
        while (i < size && F1[idx2in[i]] == F1[idx2in[from]]) { i++; }
        if (i < size) 
        {
            Population[idx2in[i]].fit_score = F1[idx2in[i]];
            Population[idx2in[i-1]].fit_score = Population[idx2in[i]].fit_score - (F1[idx2in[i]] - F1[idx2in[i-1]]) / (i-from);
        }
        else
            Population[idx2in[i-1]].fit_score = 1.0 / (Population[idx2in[i-1]].nrow+1);
        
        double delta_F = Population[idx2in[--i]].fit_score - Population[idx2in[from]].fit_score;
        double delta_F2 = F2[idx2in[i]] - F2[idx2in[from]];
        for (int j=from; j<i; j++)
            Population[idx2in[j]].fit_score = Population[idx2in[i]].fit_score - (delta_F * (F2[idx2in[i]] - F2[idx2in[j]])) / delta_F2;
    }

    sort(Population.begin(), Population.end(), [&](Individual& i0, Individual& i1){
        return i0.fit_score > i1.fit_score;
    });

    vector<double> fit_scores;
    fit_scores.reserve(this->Pc);
    for (auto &it : Population)
        fit_scores.push_back(it.fit_score);
    selection_gen = discrete_distribution<int>(fit_scores.begin(), fit_scores.begin()+Pc);
}


void GeneticAlgorithm::output_ith_best_layout()
{
    if (this->rank > Population.size())
    {
        cout << "The index referred is larger than number of layout generated" << endl;
        cout << "Output the worst layout instead" << endl;
        this->rank = Population.size();
    }

    Individual& in = Population[this->rank-1];
    vector<vector<string>> Net(nsignal+1, vector<string>({}));

    if (this->to_optimize) optimize(in);

    // horizontal
    for (int y=1; y<in.nrow-1; y++)
    {
        int x = 0, from = 0, sig;

        while (x < in.ncol)
        {
            if (in(x, y, 0) != 0)
            {
                sig = in(x, y, 0);
                from = x++;
                while (x < in.ncol && in(x, y, 0) == sig) { x++; }
                if (from < x-1) Net[sig].emplace_back(".H " + to_string(from) + " " + to_string(y) + " " + to_string(x-1));
            }
            else
                x++;
        }
    }

    // vertical
    for (int x=0; x<in.ncol; x++)
    {
        int y = 0, from = 0, sig;

        while (y < in.nrow)
        {
            if (in(x, y, 1) != 0)
            {
                sig = in(x, y, 1);
                from = y++;
                while (y < in.nrow && in(x, y, 1) == sig) { y++; }
                if (from < y-1) Net[sig].emplace_back(".V " + to_string(x) + " " + to_string(from) + " " + to_string(y-1));
            }
            else
                y++;
        }
    }

    ofstream fout(this->outfile);

    for (int i=0; i<Net.size(); i++)
    {
        if (Net[i].size() == 0) continue;

        fout << ".begin " << i << endl;
        
        for (int j=0; j<Net[i].size(); j++)
            fout << Net[i][j] << endl;

        fout << ".end" << endl;
    }

    fout.close();

    io::drawNets(this->infile.c_str(), this->outfile.c_str(), this->outfile.c_str());

    cout << "row number: " << in.nrow-2 << endl;
    cout << "total net length: " <<  in.calculate_net_length() << endl;
    cout << "total via number: " <<  in.calculate_via_number() << endl;
}


void GeneticAlgorithm::sort_solution()
{
    sort(Population.begin(), Population.end(), [&](const Individual& a, const Individual& b){
        return a.fit_score > b.fit_score;
    });
}