#pragma once
#include "typedefs.h"
#include "Vector.h"
#include "chelper.h"
#include "MathFunctions.h"
#include "Grid.h"
#include "Matrix2D.h"
#include "Source.h"
#include "Vector.h"

class CGridVoronoi: public CGridBasic
{
public:

    CGridVoronoi(void)
    {
        basic_path = 0;
        buffer_size = 0;

        max_cells = 0;
        max_value = 0;
        max_data = 0;

        min_delta = 0;
        max_delta = 0;

        min_mach = 0;
        max_mach = 0;

        min_mag = 0;
        max_mag = 0;

        min_vel = 0;
        max_vel = 0;

        min_len = 0;
        max_len = 0;

        min_gas_temp = 0;
        max_gas_temp = 0;

        min_dust_temp = 0;
        max_dust_temp = 0;

        min_gas_dens = 0;
        max_gas_dens = 0;

        min_dust_dens = 0;
        max_dust_dens = 0;

        aalg_min = 1e300;
        aalg_max = 0;

        min_larm_limit = 0;
        max_larm_limit = 0;

        min_pres = 0;
        max_pres = 0;

        line_counter = 0;
        char_counter = 0;
        ru[0] = '|';
        ru[1] = '/';
        ru[2] = '-';
        ru[3] = '\\';

        conv_length_in_SI = 1;
        conv_Bfield_in_SI = 1;
        conv_Vfield_in_SI = 1;

        nrOfGnuPoints = 1000;
        nrOfGnuVectors = 1000;
        maxGridLines = 3;

        cell_list = 0;

        data_offset = 6;
        dataID = 0;
        data_len = 0;

        mu = 2.0;

        nrOfDensRatios = 0;
        nrOfOpiateIDs = 0;

        //data_pos_gd = MAX_UINT;
        //data_pos_dd = MAX_UINT;
        //data_pos_td = MAX_UINT;
        data_pos_tg = MAX_UINT;
        data_pos_mx = MAX_UINT;
        data_pos_my = MAX_UINT;
        data_pos_mz = MAX_UINT;
        data_pos_vx = MAX_UINT;
        data_pos_vy = MAX_UINT;
        data_pos_vz = MAX_UINT;
        data_pos_px = MAX_UINT;
        data_pos_py = MAX_UINT;
        data_pos_pz = MAX_UINT;
        data_pos_aalg = MAX_UINT;
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_eq = MAX_UINT;
        data_pos_ra = MAX_UINT;

        data_pos_vt = MAX_UINT;
        data_pos_pda = MAX_UINT;

        data_pos_op = MAX_UINT;

        data_pos_n_th = MAX_UINT;
        data_pos_T_e = MAX_UINT;
        data_pos_n_cr = MAX_UINT;
        data_pos_g_min = MAX_UINT;
        data_pos_g_max = MAX_UINT;
        data_pos_p = MAX_UINT;

        pos_GasSpecRatios = 0;
        pos_OpiateIDS = 0;

        plt_gas_dens = false;
        plt_dust_dens = false;
        plt_gas_temp = false;
        plt_dust_temp = false;
        plt_mag = false;
        plt_vel = false;
        plt_rat = false;
        plt_delta = false;
        plt_larm = false;
        plt_mach = false;
        plt_dust_id = false;

        total_volume = 0;
        cell_volume = 0;

        cell_list = 0;

        line_counter = 0;
        char_counter = 0;

        mu = 2.0;

        rot_angle1 = 0;
        rot_angle2 = 0;

        hull_list = 0;
        hull_size = 0;
        stree = 0;

        min_nrOfNeigbors = uint(1e6);
        max_nrOfNeigbors = 0;
        pos_counter = 0;

        turbulent_velocity = 0;
    }

    ~CGridVoronoi()
    {
        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        if(hull_list != 0)
        {
            delete [] hull_list;
            hull_list = 0;
        }

        if(stree != 0)
        {
            delete stree;
            stree = 0;
        }

        cout << CLR_LINE << flush;
    }

    bool writeGNUPlotFiles(string path, parameter & param);

    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(cell_basic * cell)
    {
        return((cell_vo *) cell)->getCenter();
    }

    Vector3D getCenter(uint id)
    {
        cell_vo * cell = ((cell_vo *) cell_list[id]);
        return cell->getCenter();
    }

    double maxLength()
    {
        return max_len;
    }

    bool next(photon_package *pp)
    {
        if(!positionPhotonInGrid(pp))
            return false;

        if(!goToNextCellBorder(pp))
            return false;

        return true;
    };

    /*uint getXIndex(cell_vo * cell, vector<cell_vo*> & list)
    {
    uint N = list.size();
    uint min = 1, max = N;

    double x_n = cell->getCenter().X();
    double y_n = cell->getCenter().Y();
    double z_n = cell->getCenter().Z();

    if(v < x[0] || v > x[N])
    return 0;

    while(max - min > 1)
    {
    const uint i = min + (max - min) / 2;
    if(x[i] > v)
    max = i;
    else
    min = i;
    }

    return min;
    }*/



    /*void getBoundingPoints(Vector3D & p_min, Vector3D & p_max)
    {
    p_min.set(cell_oc_root->x_min, cell_oc_root->y_min,
    cell_oc_root->z_min);
    p_max.set(cell_oc_root->x_max, cell_oc_root->y_max,
    cell_oc_root->z_max);
    }

    void getBoundingPoints(cell_basic * cell, Vector3D & p_min,
    Vector3D & p_max)
    {
    cell_oc * curr_cell = (cell_oc*) cell;
    p_min.set(curr_cell->x_min, curr_cell->y_min, curr_cell->z_min);
    p_max.set(curr_cell->x_max, curr_cell->y_max, curr_cell->z_max);
    }*/

    bool findStartingPoint(photon_package *pp);

    void getLengths(uint bins, double & step_xy, double & off_xy)
    {
        step_xy = 2 * max_len / double(bins);
        off_xy = step_xy / 2.0;
    }

    bool createCellList()
    {
        cout << CLR_LINE;
        cout << "- Creating of cell list         : done          \n" << flush;
        return true;
    }

    double getVolume(cell_basic * cell)
    {
        return((cell_vo*) cell)->getVolume();
    }

    double getVolume(photon_package *pp)
    {
        cell_basic * cell_pos = pp->getPositionCell();

        return getVolume(cell_pos);
    }

    /*double getMinArea(photon_package *pp)
    {
    //tbd
    cell_oc * cell_pos = (cell_oc *) pp->getPositionCell();
    return 0;
    }*/

    bool positionPhotonInGrid(photon_package * pp);
    bool positionPhotonInGridTest(photon_package * pp);

    double getMaxLength()
    {
        return max_len;
    }

    void getAmountOfCells(uint & N_x, uint & N_y, uint & N_z)
    {
        uint nr_pixel = uint(pow(max_cells, 1.0 / 3.0));
        N_x = nr_pixel;
        N_y = nr_pixel;
        N_z = nr_pixel;
    }

    bool createArtificialGrid(string path);

    bool saveBinaryGridFile(string filename)
    {
        return saveBinaryGridFile(filename, GRID_ID_VOR, data_offset);
    }

    bool loadGridFromBinrayFile(parameter & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinrayFile(parameter & param)
    {
        return loadGridFromBinrayFile(param, 0);
    };

    void clear()
    {
        line_counter = 0;
        char_counter = 0;


        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        if(hull_list != 0)
        {
            delete [] hull_list;
            hull_list = 0;
        }

        cout << CLR_LINE;
        cout << "Final cleanup                                :  done     \n";
    }

    void printParameter();

private:
    uint pos_counter;

    class h_list
    {
    public:

        h_list()
        {
            id = 0;
        }

        h_list(double _x, double _y, double _z, uint _id)
        {
            pos = Vector3D(_x, _y, _z);
            id = _id;
        }

        Vector3D pos;
        uint id;
    };

    class tree_node
    {
    public:

        class list_element
        {
        public:

            list_element()
            {
                cell = 0;
                next = 0;
            }

            cell_vo * cell;
            list_element * next;
        };

        tree_node()
        {
            first = 0;
            last = 0;
            size = 0;

            x_min = 0;
            y_min = 0;
            z_min = 0;
            length = 0;

            level = 0;
            branch = 0;
            leafs = 0;
        }

        ~tree_node()
        {
            clear();
        }

        bool isInNode(Vector3D point)
        {
            if(point.X() < getXMin() || point.Y() < getYMin() || point.Z() < getZMin())
                return false;

            if(point.X() > getXMax() || point.Y() > getYMax() || point.Z() > getZMax())
                return false;

            return true;
        }

        void add_cell(cell_vo * cell)
        {
            list_element * new_element = new list_element();
            new_element->cell = cell;

            if(first == 0)
            {
                first = new_element;
                last = new_element;
            }
            else
            {
                last->next = new_element;
                last = new_element;
            }

            size++;
        }

        void clear()
        {
            if(first == 0)
                return;

            list_element * pos = first;

            while(pos != 0)
            {
                list_element * tmp_element = pos;
                pos = pos->next;

                delete tmp_element;
            }

            first = 0;
            last = 0;
            size = 0;
        }

        cell_vo * findClosestCell(Vector3D point, double & _min_distance)
        {
            list_element * pos = first;
            double min_distance = 1e200;
            cell_vo * res = 0;

            double X = point.X();
            double Y = point.Y();
            double Z = point.Z();

            while(pos != 0)
            {
                cell_vo * cell = pos->cell;
                Vector3D center = cell->getCenter();

                double X1 = center.X();
                double Y1 = center.Y();
                double Z1 = center.Z();

                double sq_distance = (X - X1)*(X - X1) + (Y - Y1)*(Y - Y1) + (Z - Z1)*(Z - Z1);

                if(sq_distance<min_distance)
                {
                    min_distance = sq_distance;
                    res = cell;
                }

                pos = pos->next;
            }

            _min_distance = min_distance;
            return res;
        }

        uint get_size()
        {
            return size;
        }

        bool is_emty()
        {
            return size == 0;
        }

        double getXMin()
        {
            return x_min;
        }

        double getYMin()
        {
            return y_min;
        }

        double getZMin()
        {
            return z_min;
        }

        double getXMax()
        {
            return x_min + length;
        }

        double getYMax()
        {
            return y_min + length;
        }

        double getZMax()
        {
            return z_min + length;
        }

        double getXCenter()
        {
            return x_min + 0.5 * length;
        }

        double getYCenter()
        {
            return y_min + 0.5 * length;
        }

        double getZCenter()
        {
            return z_min + 0.5 * length;
        }

        void setXMin(double x)
        {
            x_min = x;
        }

        void setYMin(double y)
        {
            y_min = y;
        }

        void setZMin(double y)
        {
            z_min = y;
        }

        double getLength()
        {
            return length;
        }

        void setLength(double l)
        {
            length = l;
        }

        void setLevel(uint l)
        {
            level = l;
        }

        uint getLevel()
        {
            return level;
        }

        tree_node * getLeafs()
        {
            return leafs;
        }

        tree_node * getLeaf(uint index)
        {
            return &leafs[index];
        }

        tree_node * getBranch()
        {
            return branch;
        }

        void setLeafs(tree_node * l)
        {
            leafs = l;
        }

        void setBranch(tree_node * b)
        {
            branch = b;
        }

    private:
        list_element * first;
        list_element * last;

        uint size;

        uint level;

        double x_min;
        double y_min;
        double z_min;
        double length;

        tree_node * branch;
        tree_node * leafs;
    };

    class search_tree
    {
    public:

        search_tree()
        {
            root = 0;
            max_level = 0;
        };

        ~search_tree()
        {
            clear(root);
        };

        bool addCell(cell_vo* cell)
        {
            Vector3D center = cell->getCenter();

            if(center.X() < root->getXMin() || center.Y() < root->getYMin()
                    || center.Z() < root->getZMin())
                return false;

            if(center.X() > root->getXMax() || center.Y() > root->getYMax()
                    || center.Z() > root->getZMax())
                return false;

            return addCell(root, cell);
        }

        void initTree(ulong nrOfCells, double _side_length)
        {
            side_length = _side_length;
            max_level = calcMaxLevel(nrOfCells);

            root = new tree_node();
            root->setLength(side_length);
            root->setXMin(-0.5 * side_length);
            root->setYMin(-0.5 * side_length);
            root->setZMin(-0.5 * side_length);

            createLeafNodes(root, max_level);
        }

        cell_vo * findClosestCell(Vector3D pos, cell_basic ** cell_list)
        {
            //find current node
            tree_node * node = findMatchingNode(pos, max_level);
            double min_distance = 0;

            //no node found
            if(node == 0)
                return 0;

            //when node has not enough cells go to next lower level
            while(node->get_size() < 20)
            {
                node = node->getBranch();

                if(node->getBranch() == 0)
                    break;
            }

            //find closest cell in current node
            cell_vo * cell = node->findClosestCell(pos, min_distance);
            uint level = node->getLevel();
            vector<tree_node*> node_list;
            tree_node * neigbour_node = 0;

            Vector3D center = cell->getCenter();

            double px=pos.X();
            double py=pos.Y();
            double pz=pos.Z();

            //find possible closer cells in neighboring nodes
            double len = min_distance * 1.001;

            for(int i_x=-1;i_x<=1;i_x++)
            {
                for(int i_y=-1;i_y<=1;i_y++)
                {
                    for(int i_z=-1;i_z<=1;i_z++)
                    {
                        if(i_x+i_y+i_z==0)
                            continue;

                        Vector3D tmp_pos=Vector3D(px + double(i_x)*len, py + double(i_y)*len, pz + double(i_z)*len);

                        if(!node->isInNode(tmp_pos))
                        {
                            neigbour_node = findMatchingNode(tmp_pos, level);

                            if(neigbour_node != 0)
                            {
                                if(node!=neigbour_node)
                                    node_list.push_back(neigbour_node);
                            }
                        }
                    }
                }
            }

            //eliminate double entries
            if(node_list.size() > 2)
            {
                sort(node_list.begin(), node_list.end());
                node_list.erase(unique(node_list.begin(), node_list.end()), node_list.end());
            }

            //check all neighboring cells for shortest distance
            double n_distance = 2e200;
            cell_vo * n_cell = 0;

            for(uint i = 0; i < node_list.size(); i++)
            {
                tree_node * tmp_node = node_list[i];

                if(tmp_node->is_emty())
                    continue;

                double sq_distance = 1e200;
                cell_vo * tmp_cell = tmp_node->findClosestCell(pos, sq_distance);

                if(sq_distance < n_distance)
                {
                    n_distance = sq_distance;
                    n_cell = tmp_cell;
                }
            }

            if(n_distance < min_distance)
            {
                cell = n_cell;
                min_distance = n_distance;
            }

            uint n_size = cell->getNrOfNeighbors();

            center = cell->getCenter();
            bool found = false;

            for(uint i = 0; i < n_size; i++)
            {
                if(cell->getNeighborID(i)>-1)
                {
                    uint id = cell->getNeighborID(i);
                    cell_vo * tmp_cell = ((cell_vo*) cell_list[id]);
                    Vector3D diff_pos = pos - tmp_cell->getCenter();
                    double len = diff_pos.sq_length();

                    if(len < min_distance)
                    {
                        n_cell = tmp_cell;
                        min_distance=len;
                        found = true;
                    }
                }
            }

            if(found)
                cell = n_cell;/**/

            return cell;
        }

    private:

        tree_node * findMatchingNode(Vector3D point)
        {
            return findMatchingNode(point, max_level);
        }

        tree_node * findMatchingNode(Vector3D point, uint _level)
        {
            if(point.X() < root->getXMin() || point.Y() < root->getYMin()
                    || point.Z() < root->getZMin())
                return 0;

            if(point.X() > root->getXMax() || point.Y() > root->getYMax()
                    || point.Z() > root->getZMax())
                return 0;

            return goNextLevelDown(root, point, _level);
        }

        uint calcMaxLevel(ulong nrOfCells)
        {
            uint level = 1;

            double nc = double(nrOfCells) / pow(8.0, double(level));

            while(nc > 30.0)
            {
                level++;
                nc = double(nrOfCells) / pow(8.0, double(level));
            }

            return level+2;
        }

        bool addCell(tree_node * node, cell_vo* cell)
        {
            Vector3D center = cell->getCenter();
            double X = center.X();
            double Y = center.Y();
            double Z = center.Z();

            node->add_cell(cell);

            if(node->getLeafs() == 0)
                return true;

            double xmid = node->getXCenter();
            double ymid = node->getYCenter();
            double zmid = node->getZCenter();

            if(Z < zmid) //z 0 1 2 3
            {
                if(Y < ymid) //y 0 1
                {
                    if(X < xmid) //x 0
                        node = node->getLeaf(0);
                    else
                        //x 1
                        node = node->getLeaf(1);
                }
                else //y 2 3
                {
                    if(X < xmid) //x 2
                        node = node->getLeaf(2);
                    else //x 3
                        node = node->getLeaf(3);
                }
            }
            else //z 4 5 6 7
            {
                if(Y < ymid) //y 4 5
                {
                    if(X < xmid) //x 4
                        node = node->getLeaf(4);
                    else //x 5
                        node = node->getLeaf(5);
                }
                else //y 6 7
                {
                    if(X < xmid) //x 6
                        node = node->getLeaf(6);
                    else //x 7
                        node = node->getLeaf(7);
                }
            }

            return addCell(node, cell);
        }

        tree_node * goNextLevelDown(tree_node * node, Vector3D point, uint _level)
        {
            if(node->getLeafs() == 0)
                return node;

            if(node->getLevel() == _level)
                return node;


            double xmid = node->getXCenter();
            double ymid = node->getYCenter();
            double zmid = node->getZCenter();

            if(point.Z() < zmid) //z 0 1 2 3
            {
                if(point.Y() < ymid) //y 0 1
                {
                    if(point.X() < xmid) //x 0
                        node = node->getLeaf(0);
                    else
                        //x 1
                        node = node->getLeaf(1);
                }
                else //y 2 3
                {
                    if(point.X() < xmid) //x 2
                        node = node->getLeaf(2);
                    else //x 3
                        node = node->getLeaf(3);
                }
            }
            else //z 4 5 6 7
            {
                if(point.Y() < ymid) //y 4 5
                {
                    if(point.X() < xmid) //x 4
                        node = node->getLeaf(4);
                    else //x 5
                        node = node->getLeaf(5);
                }
                else //y 6 7
                {
                    if(point.X() < xmid) //x 6
                        node = node->getLeaf(6);
                    else //x 7
                        node = node->getLeaf(7);
                }
            }

            return goNextLevelDown(node, point, _level);
        }

        void clear(tree_node * cell)
        {
            tree_node * leafs = cell->getLeafs();

            if(leafs == 0)
                return;

            for(uint i = 0; i < 8; i++)
            {
                clear(&leafs[i]);
            }

            delete[] leafs;
            cell->setLeafs(0);
        }

        void createLeafNodes(tree_node * node, uint max_level)
        {
            uint next_level = 1 + node->getLevel();

            if(max_level < next_level)
                return;

            double ox = node->getXMin();
            double oy = node->getYMin();
            double oz = node->getZMin();

            double length = 0.5 * node->getLength();
            tree_node * leafs = new tree_node[8];

            leafs[0].setBranch(node);
            leafs[0].setXMin(ox);
            leafs[0].setYMin(oy);
            leafs[0].setZMin(oz);

            leafs[0].setLength(length);
            leafs[0].setLevel(next_level);


            leafs[1].setBranch(node);
            leafs[1].setXMin(ox + length);
            leafs[1].setYMin(oy);
            leafs[1].setZMin(oz);

            leafs[1].setLength(length);
            leafs[1].setLevel(next_level);


            leafs[2].setBranch(node);
            leafs[2].setXMin(ox);
            leafs[2].setYMin(oy + length);
            leafs[2].setZMin(oz);

            leafs[2].setLength(length);
            leafs[2].setLevel(next_level);


            leafs[3].setBranch(node);
            leafs[3].setXMin(ox + length);
            leafs[3].setYMin(oy + length);
            leafs[3].setZMin(oz);

            leafs[3].setLength(length);
            leafs[3].setLevel(next_level);


            leafs[4].setBranch(node);
            leafs[4].setXMin(ox);
            leafs[4].setYMin(oy);
            leafs[4].setZMin(oz + length);

            leafs[4].setLength(length);
            leafs[4].setLevel(next_level);


            leafs[5].setBranch(node);
            leafs[5].setXMin(ox + length);
            leafs[5].setYMin(oy);
            leafs[5].setZMin(oz + length);

            leafs[5].setLength(length);
            leafs[5].setLevel(next_level);


            leafs[6].setBranch(node);
            leafs[6].setXMin(ox);
            leafs[6].setYMin(oy + length);
            leafs[6].setZMin(oz + length);

            leafs[6].setLength(length);
            leafs[6].setLevel(next_level);


            leafs[7].setBranch(node);
            leafs[7].setXMin(ox + length);
            leafs[7].setYMin(oy + length);
            leafs[7].setZMin(oz + length);

            leafs[7].setLength(length);
            leafs[7].setLevel(next_level);


            node->setLeafs(leafs);

            for(uint i = 0; i < 8; i++)
                createLeafNodes(&leafs[i], max_level);
        }

        double side_length;
        uint max_level;

        tree_node * root;
    };

    uint min_nrOfNeigbors;
    uint max_nrOfNeigbors;

    uint hull_size;
    h_list * hull_list;
    search_tree * stree;

    static inline bool cell_compID(cell_vo * v1, cell_vo * v2)
    {
        if(v1->getID() < v2->getID())
            return true;

        return false;
    }

    double abs_min(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
    {
        double res = 1e300;

        if(res > abs(x_min))
            res = abs(x_min);

        if(res > abs(x_max))
            res = abs(x_max);

        if(res > abs(y_min))
            res = abs(y_min);

        if(res > abs(y_max))
            res = abs(y_max);

        if(res > abs(z_min))
            res = abs(z_min);

        if(res > abs(z_max))
            res = abs(z_max);

        return res;
    }

    double abs_max(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
    {
        double res = 0;

        if(res < abs(x_min))
            res = abs(x_min);

        if(res < abs(x_max))
            res = abs(x_max);

        if(res < abs(y_min))
            res = abs(y_min);

        if(res < abs(y_max))
            res = abs(y_max);

        if(res < abs(z_min))
            res = abs(z_min);

        if(res < abs(z_max))
            res = abs(z_max);

        return res;
    }

    /*bool find_hull_point(photon_package * pp)
    {
        bool found = false;
        double max_radius = 2e300;
        Vector3D pos = pp->getPosition();

        for(uint i = 0; i < hull_size; i++)
        {
            uint id = hull_list[i].id;
            cell_vo * tmp_cell = ((cell_vo*) cell_list[id]);
            Vector3D tmp_pos = pos - tmp_cell->getCenter();
            double len = tmp_pos.sq_length();

            if(len < max_radius)
            {
                max_radius = len;
                pp->setPositionCell(tmp_cell);
                found = true;
            }

            find_neighboring_cell(pp);
        }

        return found;
    }*/

    /*bool find_neighboring_cell(photon_package * pp)
    {
        cell_vo * center_cell = (cell_vo*) pp->getPositionCell();

        if(center_cell == 0)
            return false;

        Vector3D pos = pp->getPosition();

        Vector3D cell_pos = center_cell->getCenter();

        Vector3D diff_pos = pos - cell_pos;

        double max_radius = diff_pos.sq_length();

        uint c_size = center_cell->getNrOfNeighbors();
        uilist id_list;

        for(uint i = 0; i < c_size; i++)
        {
            if(isNeigboringVoroCell(center_cell, i))
            {
                uint id = center_cell->getNeighborID(i);
                cell_vo * tmp_cell = ((cell_vo*) cell_list[id]);
                diff_pos = pos - tmp_cell->getCenter();
                double len = diff_pos.sq_length();

                uint n_size = tmp_cell->getNrOfNeighbors();

                for(uint j = 0; j < n_size; j++)
                {
                    if(isNeigboringVoroCell(tmp_cell, j))
                    {
                        uint id = tmp_cell->getNeighborID(j);
                        insertInList(id, id_list);
                    }
                }

                if(len <= max_radius)
                {
                    max_radius = len;
                    pp->setPositionCell(tmp_cell);
                }
            }
        }

        for(uint i = 0; i < id_list.size(); i++)
        {
            uint id = id_list[i];
            cell_vo * tmp_cell = ((cell_vo*) cell_list[id]);
            diff_pos = pos - tmp_cell->getCenter();
            double len = diff_pos.sq_length();

            if(len <= max_radius)
            {
                max_radius = len;
                pp->setPositionCell(tmp_cell);
            }
        }

        return true;
    }*/

    bool isHullPoint(uint id)
    {
        uint N = hull_size;
        uint min = 0, max = N - 1;

        if(id < hull_list[min].id)
            return false;

        if(id > hull_list[max].id)
            return false;

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(hull_list[i].id > id)
                max = i;
            else
                min = i;
        }

        if(id == hull_list[min].id)
            return true;

        uint upper = min + 1;

        if(upper > N)
            upper = N;

        if(id == hull_list[upper].id)
            return true;

        uint lower = min - 1;

        if(lower == uint(-1))
            lower = 0;

        if(id == hull_list[lower].id)
            return true;

        return false;
    }

    uint biListIndexSearch(uint val, uilist & list)
    {
        uint N = uint(list.size());
        uint min = 0, max = N - 1;

        if(val < list[0] || val > list[max])
            return uint(-1);

        while(max - min > 1)
        {
            uint i = min + (max - min) / 2;
            if(list[i] >= val)
                max = i;
            else
                min = i;
        }

        return min;
    }

    void insertInList(uint val, uilist& list)
    {
        uint N = uint(list.size());

        if(N == 0)
        {
            list.push_back(val);
            return;
        }

        if(N == 1)
        {
            if(list[0] == val)
                return;
        }

        if(val < list[0])
        {
            uilist::iterator it = list.begin();
            list.insert(it, val);
            return;
        }

        if(val > list[N - 1])
        {
            list.push_back(val);
            return;
        }

        uilist::iterator it = list.begin();
        uint index = biListIndexSearch(val, list);

        if(list[index] == val)
            return;

        if(list[index + 1] == val)
            return;

        list.insert(it + index + 1, val);
    }

    bool isInside(Vector3D & pos, cell_basic * _cell)
    {
        cell_vo * cell = (cell_vo *) _cell;
        cout << "WARNING: This function needs to be implimented if needed!     \n";

        /*if(cell->getRID() == uint(-1))
        {
        if(pos.sq_length() > Rmin * Rmin)
        return false;
        }

        double r1 = listR[cell->getRID()];
        double r2 = listR[cell->getRID() + 1];
        double ph1 = listPhi[cell->getPhID()];
        double ph2 = listPhi[cell->getPhID() + 1];
        double th1 = listTheta[cell->getThID()];
        double th2 = listTheta[cell->getThID() + 1];

        Vector3D tmp_pos = pos.getSphericalCoord();

        if(tmp_pos.R() < r1)
        return false;

        if(tmp_pos.R() > r2)
        return false;

        if(tmp_pos.Phi() < ph1)
        return false;

        if(tmp_pos.Phi() > ph2)
        return false;

        if(tmp_pos.Theta() < th1)
        return false;

        if(tmp_pos.Theta() > th2)
        return false;*/

        return true;
    }

    bool isInside(Vector3D & pos)
    {
        double l_min = -0.5 * max_len;
        double l_max = 0.5 * max_len;

        if(pos.X() < l_min)
            return false;
        if(pos.Y() < l_min)
            return false;
        if(pos.Z() < l_min)
            return false;

        if(pos.X() > l_max)
            return false;
        if(pos.Y() > l_max)
            return false;
        if(pos.Z() > l_max)
            return false;

        return true;
    }

    bool isInside(photon_package * pp, Vector3D & pos)
    {
        return isInside(pos, pp->getPositionCell());
    }

    void setRndPositionInCell(photon_package * pp)
    {
        cell_vo * tmp_cell = (cell_vo*) pp->getPositionCell();
        Vector3D pos = getCenter(tmp_cell);
        pp->setPosition(pos);
    }

    void addGNULines(uint cID, stringstream & str)
    {
        cell_vo * tmp_cell = (cell_vo*) cell_list[cID];
        uint nr_neighbors = tmp_cell->getNrOfNeighbors();

        Vector3D p1 = tmp_cell->getCenter();

        for(uint i = 0; i < nr_neighbors; i++)
        {
            if(isNeigboringVoroCell(tmp_cell, i))
            {
                Vector3D p2 = getNeighborCenter(tmp_cell, i);
                Vector3D tmp_len = p2 - p1;
                str << p1.X() << " " << p1.Y() << " " << p1.Z() << " " << tmp_len.X() << " " << tmp_len.Y() << " " << tmp_len.Z() << "\n";
            }
        }

    }

    bool isNeigboringVoroCell(cell_vo* cell, uint nID)
    {
        int id = cell->getNeighborID(nID);
        
        if(id>max_cells)
            return false;
            
        return id>-1;
    }

    Vector3D getNeighborCenter(cell_vo* cell, uint nID)
    {
        int id = cell->getNeighborID(nID);
        cell_vo* n_cell = ((cell_vo *) cell_list[id]);
        return n_cell->getCenter();
    }

};
