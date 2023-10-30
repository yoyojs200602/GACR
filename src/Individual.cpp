#include "Individual.hpp"


int Individual::predict_extension_row(vector<int>& s, vector<int>& t, int target_sig, bounding_box& bb)
{
    if (s[1] == nrow-1 && this->operator()(s[0], s[1]-1, 1) != 0 && this->operator()(s[0], s[1]-1, 1) != target_sig) return nrow-1;
    if (s[1] == 0 && this->operator()(s[0], s[1]+1, 1) != 0 && this->operator()(s[0], s[1]+1, 1) != target_sig) return 1;
    if (t[1] == nrow-1 && this->operator()(t[0], t[1]-1, 1) != 0 && this->operator()(t[0], t[1]-1, 1) != target_sig) return nrow-1;
    if (t[1] == 0 && this->operator()(t[0], t[1]+1, 1) != 0 && this->operator()(t[0], t[1]+1, 1) != target_sig) return 1;
    
    int random = rand()%4;
    if (random == 0)
        return bb.bottom+1;
    else if (random == 1)
        return bb.top-1;
    else
        return rand() % (bb.top-bb.bottom-1) + bb.bottom + 1;
}


int Individual::calculate_net_length()
{
    int length = 0;

    // horizontal net length
    for (int i = 1; i < nrow-1; i++) {
        int tmp_sig = 0;
        for (int j = 0; j < ncol; j++) {
            if (tmp_sig == 0) {
                tmp_sig = this->operator()(j, i, 0);
            } else if (this->operator()(j, i, 0) == tmp_sig) {
                length++;
            } else {
                tmp_sig = this->operator()(j, i, 0);
            }
        }
    }

    // vertical net length
    for (int i = 0; i < ncol; i++) {
        int tmp_sig = 0;
        for (int j = 0; j < nrow; j++) {
            if (tmp_sig == 0) {
                tmp_sig = this->operator()(i, j, 1);
            } else if (this->operator()(i, j, 1) == tmp_sig) {
                length++;
            } else {
                tmp_sig = this->operator()(i, j, 1);
            }
        }
    }

    return length;
}


int Individual::calculate_via_number()
{
    int cnt = 0;
    for (int i=1; i<nrow-1; i++)
    {
        auto &row = this->get_row(i);
        for (int j=0; j<ncol; j++)
            if (row[j][0] == row[j][1])
                cnt++;
    }

    return cnt;
}


void Individual::add_channel(int n)
{
    vector<int> zero(2, 0);
    g.emplace(g.begin()+n, ncol, zero);

    if (nrow > 2)
        for (int i=0; i<ncol; i++)
            if (this->operator()(i, n+1, 1) == this->operator()(i, n-1, 1))
                this->operator()(i, n, 1) = this->operator()(i, n+1, 1);

    this->nrow++;
    visit_table.addRows(1);
}


void Individual::add_random_channels(int n)
{
    nrow += n;
    vector<bool> channel_order(nrow, false);
    vector<int> zero(2, 0);
    for (int i=1; i<n+1; i++) channel_order[i] = true;
    shuffle(channel_order.begin()+1, channel_order.end()-1, mt19937{random_device{}()});
    
    g.reserve(nrow);

    for (int i=1; i<nrow-1; i++)
    {
        if (channel_order[i])
        {
            g.emplace(g.begin()+i, ncol, zero);
            for (int j=0; j<ncol; j++)
                if (this->operator()(j, i-1, 1) == this->operator()(j, i+1, 1))
                    this->operator()(j, i, 1) = this->operator()(j, i-1, 1);
        }
    }

    this->visit_table.addRows(n);
}


void Individual::delete_channel(int n)
{
    g.erase(g.begin() + n);
    this->visit_table.deleteRows(1);
    nrow--;
}


void Individual::delete_channels(vector<int>& remove_list)
{
    for (auto it=remove_list.rbegin(); it!=remove_list.rend(); ++it)
        g.erase(g.begin() + *it);
    this->visit_table.deleteRows(remove_list.size());
    nrow -= remove_list.size();
}


void Individual::delete_channels(bounding_box& bb)
{
    for (int i=bb.bottom+1; i<bb.top; i++)
    {
        auto& row = this->get_row(i);
        for (auto it=row.begin()+bb.left+1; it!=row.begin()+bb.right; it++)
            (*it)[0] = (*it)[1] = 0;
    }
}


void Individual::remove_if_isolated(vector<int>& pt) 
{
    int target_in_sig = this->operator()(pt[0], pt[1], pt[2]);

    if (pt[2] == 0)
        if (!((pt[0] > 0 && this->operator()(pt[0]-1, pt[1], 0) == target_in_sig) || 
            (pt[0] < ncol-1 && this->operator()(pt[0]+1, pt[1], 0) == target_in_sig)))
            this->operator()(pt[0], pt[1], 0) = 0;
    else
        if (!((pt[1] > 0 && this->operator()(pt[0], pt[1]-1, 1) == target_in_sig) || 
            (pt[1] < nrow-1 && this->operator()(pt[0], pt[1]+1, 1) == target_in_sig)))
            this->operator()(pt[0], pt[1], 1) = 0;
}


vector<int> Individual::find_row_index(grid& pt)
{
    int x = get<0>(pt), direc = get<2>(pt);
    int* grid_addr = get<1>(pt);

    for (int y=0; y<nrow; y++)
        if (&(this->operator()(x, y, direc)) == grid_addr)
            return {x, y, direc};

    return {};
}