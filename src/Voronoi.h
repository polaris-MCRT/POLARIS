#pragma once
#include "Grid.h"
#include "Vector.h"

// search tree parameters
#define MAX_CELLS 31 // max. cells per tree node
#define MAX_LEVEL 30 // max. tree level

class CGridVoronoi : public CGridBasic
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

        // data_pos_gd = MAX_UINT;
        // data_pos_dd = MAX_UINT;
        // data_pos_td = MAX_UINT;
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
        // data_pos_aalg = MAX_UINT;
        data_pos_amin = MAX_UINT;
        data_pos_amax = MAX_UINT;
        data_pos_size_param = MAX_UINT;
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
        plt_mol_dens = false;
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
            delete[] hull_list;
            hull_list = 0;
        }

        if(stree != 0)
        {
            delete stree;
            stree = 0;
        }

        cout << CLR_LINE << flush;
    }

    bool isInside(const Vector3D & pos)const
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

    bool writeGNUPlotFiles(string path, parameters & param);

    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(const cell_basic & cell) const
    {
        return ((const cell_vo *)&cell)->getCenter();
    }

    Vector3D getCenter(uint id) const
    {
        const cell_vo * cell = ((const cell_vo *)cell_list[id]);
        return cell->getCenter();
    }

    double maxLength()
    {
        return max_len;
    }

    bool next(photon_package * pp)
    {
        if(!positionPhotonInGrid(pp))
            return false;

        if(!goToNextCellBorder(pp))
            return false;

        return true;
    };

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy)
    {
        step_xy = 2 * max_len / double(bins);
        off_xy = step_xy / 2.0;
    }

    bool createCellList()
    {
        // cout << CLR_LINE;
        // cout << "- Creating of cell list                : done          \n" << flush;
        return true;
    }

    double getVolume(const cell_basic & cell) const
    {
        const cell_vo * cell_pos = (const cell_vo *)&cell;
        return cell_pos->getVolume();
    }

    bool positionPhotonInGrid(photon_package * pp);

    // for debugging only
    bool positionPhotonInGridTest(photon_package * pp);
    bool createArtificialGrid(string path);

    double getMaxLength()
    {
        return max_len;
    }

    bool saveBinaryGridFile(string filename)
    {
        return saveBinaryGridFile(filename, GRID_ID_VOR, data_offset);
    }

    bool loadGridFromBinrayFile(parameters & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinrayFile(parameters & param)
    {
        return loadGridFromBinrayFile(param, 0);
    };

    // final cleanup
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
            delete[] hull_list;
            hull_list = 0;
        }

        cout << CLR_LINE;
        cout << "Final cleanup                                : done     \n";
    }

    void printParameters();

  private:
    uint pos_counter;

    // list of convex hull points
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

    // search tree class
    class search_tree
    {
      public:
        // list element of for linked Voronoi cell lists
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

        // node object for the search tree
        class tree_node
        {
          public:
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

            // finds all nodes within search radius
            bool nodeIntersection(Vector3D point, double _sq_distance)
            {
                double X = point.X();
                double Y = point.Y();
                double Z = point.Z();

                // is point within node?
                if(isPointInNode(point))
                    return true;

                // edge intersection
                for(uint ix = 0; ix <= 1; ix++)
                    for(uint iy = 0; iy <= 1; iy++)
                        for(uint iz = 0; iz <= 1; iz++)
                        {
                            double X1 = x_min + double(ix) * length;
                            double Y1 = y_min + double(iy) * length;
                            double Z1 = z_min + double(iz) * length;

                            double sq_distance =
                                (X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) + (Z - Z1) * (Z - Z1);

                            if(sq_distance <= _sq_distance)
                                return true;
                        }

                // surface intersection
                double distance = sqrt(_sq_distance);
                for(int ix = -1; ix <= 1; ix++)
                    for(int iy = -1; iy <= 1; iy++)
                        for(int iz = -1; iz <= 1; iz++)
                        {
                            if(ix + iy + iz == 0)
                                continue;

                            double x = point.X() + double(ix) * distance;
                            double y = point.Y() + double(iy) * distance;
                            double z = point.Z() + double(iz) * distance;

                            if(isPointInNode(Vector3D(x, y, z)))
                                return true;
                        } /**/

                return false;
            }

            bool isPointInNode(Vector3D point)
            {
                if(point.X() < getXMin() || point.Y() < getYMin() || point.Z() < getZMin())
                    return false;

                if(point.X() > getXMax() || point.Y() > getYMax() || point.Z() > getZMax())
                    return false;

                return true;
            }

            // keep track of the amount of cells in each level
            void increment()
            {
                size++;
            }

            // add cell pointer to linked list
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

            // final cleanup
            void clear()
            {
                if(first == 0)
                    return;

                list_element * pos = first;

                while(pos != 0)
                {
                    list_element * tmp_element = pos;
                    pos = pos->next;

                    // only for debugging
                    // final cell cleanup happens in Voronoi grid itself
                    // delete tmp_element->cell;
                    // tmp_element->cell=0;

                    delete tmp_element;
                    tmp_element = 0;
                }

                first = 0;
                last = 0;
                size = 0;
            }

            // empty node and refine tree if to many cell populate the node
            void push_down()
            {
                list_element * pos = first;

                while(pos != 0)
                {
                    list_element * tmp_element = pos;
                    pos = pos->next;

                    Vector3D center = tmp_element->cell->getCenter();
                    double X = center.X();
                    double Y = center.Y();
                    double Z = center.Z();

                    double xmid = this->getXCenter();
                    double ymid = this->getYCenter();
                    double zmid = this->getZCenter();

                    tree_node * leaf;

                    if(Z < zmid) // z 0 1 2 3
                    {
                        if(Y < ymid) // y 0 1
                        {
                            if(X < xmid) // x 0
                                leaf = this->getLeaf(0);
                            else
                                // x 1
                                leaf = this->getLeaf(1);
                        }
                        else // y 2 3
                        {
                            if(X < xmid) // x 2
                                leaf = this->getLeaf(2);
                            else // x 3
                                leaf = this->getLeaf(3);
                        }
                    }
                    else // z 4 5 6 7
                    {
                        if(Y < ymid) // y 4 5
                        {
                            if(X < xmid) // x 4
                                leaf = this->getLeaf(4);
                            else // x 5
                                leaf = this->getLeaf(5);
                        }
                        else // y 6 7
                        {
                            if(X < xmid) // x 6
                                leaf = this->getLeaf(6);
                            else // x 7
                                leaf = this->getLeaf(7);
                        }
                    }

                    leaf->add_cell(tmp_element->cell);

                    delete tmp_element;
                    tmp_element = 0;
                }

                first = 0;
                last = 0;
            }

            // find closest cell for given point
            cell_vo * findClosestCell(Vector3D point, double & _min_distance, tree_node *& f_node)
            {
                double min_distance = 1e200;
                cell_vo * res = 0;

                // search in highest level
                if(this->getLeafs() != 0)
                {
                    for(uint i = 0; i < 8; i++)
                    {
                        double tmp_distance = 0;
                        cell_vo * cell = leafs[i].findClosestCell(point, tmp_distance, f_node);

                        if(tmp_distance < min_distance)
                        {
                            min_distance = tmp_distance;
                            res = cell;
                        }
                    }
                }
                else // search in level plus all sub-levels within search radius
                {
                    list_element * pos = first;

                    while(pos != 0)
                    {
                        cell_vo * cell = pos->cell;
                        Vector3D center = cell->getCenter();

                        double X = point.X();
                        double Y = point.Y();
                        double Z = point.Z();

                        double X1 = center.X();
                        double Y1 = center.Y();
                        double Z1 = center.Z();

                        double sq_distance = (X - X1) * (X - X1) + (Y - Y1) * (Y - Y1) + (Z - Z1) * (Z - Z1);

                        if(sq_distance < min_distance)
                        {
                            min_distance = sq_distance;
                            res = cell;

                            f_node = this;
                        }

                        pos = pos->next;
                    }
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

            double getLength() const
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

            list_element * get_first()
            {
                return first;
            }

            list_element * get_last()
            {
                return last;
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

        search_tree()
        {
            root = 0;
            max_level = 0;
            max_nodes = 0;
        };

        ~search_tree()
        {
            if(root != 0)
            {
                clear(root);
                root = 0;
            }
        };

        bool addCell(cell_vo * cell)
        {
            Vector3D center = cell->getCenter();

            if(center.X() < root->getXMin() || center.Y() < root->getYMin() || center.Z() < root->getZMin())
                return false;

            if(center.X() > root->getXMax() || center.Y() > root->getYMax() || center.Z() > root->getZMax())
                return false;

            return addCell(root, cell);
        }

        // init. tree and first level refinement
        void initTree(double _side_length)
        {
            side_length = _side_length;

            root = new tree_node();
            root->setLength(side_length);
            root->setXMin(-0.5 * side_length);
            root->setYMin(-0.5 * side_length);
            root->setZMin(-0.5 * side_length);

            max_nodes = 1;

            createLeafNodes(root);
        }

        uint getMaxLevel()
        {
            return max_level;
        }

        uint getMaxNodes()
        {
            return max_nodes;
        }

        // find closest cell for given point in the entire tree
        cell_vo * findClosestCell(Vector3D point, cell_basic ** cell_list)
        {
            // find current node that contains the point
            tree_node * p_node = findMatchingNode(point, MAX_LEVEL);
            double min_distance = 0;

            // no node found
            if(p_node == 0)
                return 0;

            // when node has not enough cells go to next lower level
            while(p_node->is_emty())
            {
                p_node = p_node->getBranch();

                if(p_node->getBranch() == 0)
                    break;
            }

            // node that contains the closest cell
            // not necessarily identical with p_node
            tree_node * f_node = 0;

            // find closest cell in current node
            cell_vo * cell = p_node->findClosestCell(point, min_distance, f_node);

            // find possible closer cells in neighboring nodes
            double len = min_distance * 1.00001;
            double n_distance = 1e200;
            cell_vo * n_cell = checkNeighboringNodes(root, p_node, f_node, point, len, n_distance);

            if(n_distance < min_distance)
            {
                cell = n_cell;
                min_distance = n_distance;
            }

            return cell;
        }

      private:
        // check neighboring nodes for shortest distance
        cell_vo * checkNeighboringNodes(tree_node * node,
                                        tree_node * f_node,
                                        tree_node * p_node,
                                        Vector3D point,
                                        double distance,
                                        double & min_distance)
        {
            cell_vo * res_cell = 0;
            if(node->getLeafs() == 0)
            {
                // do not search f_node again
                if(node != f_node)
                {
                    // do not search p_node again
                    if(node != p_node)
                    {
                        tree_node * dummy;
                        res_cell = node->findClosestCell(point, min_distance, dummy);
                    }
                }
            }
            else
            {
                for(uint i = 0; i < 8; i++)
                {
                    if(node->getLeaf(i)->nodeIntersection(point, distance))
                    {
                        double tmp_distance = 1e200;
                        cell_vo * tmp_cell = checkNeighboringNodes(
                            node->getLeaf(i), f_node, p_node, point, distance, tmp_distance);

                        if(tmp_distance < min_distance)
                        {
                            min_distance = tmp_distance;
                            res_cell = tmp_cell;
                        }
                    }
                }
            }

            return res_cell;
        }

        tree_node * findMatchingNode(Vector3D point)
        {
            return findMatchingNode(point, MAX_LEVEL);
        }

        // go to maximal level for starting the point search
        tree_node * findMatchingNode(Vector3D point, uint _level)
        {
            if(point.X() < root->getXMin() || point.Y() < root->getYMin() || point.Z() < root->getZMin())
                return 0;

            if(point.X() > root->getXMax() || point.Y() > root->getYMax() || point.Z() > root->getZMax())
                return 0;

            return goNextLevelDown(root, point, _level);
        }

        // add cell to tree
        bool addCell(tree_node * node, cell_vo * cell)
        {
            Vector3D center = cell->getCenter();
            double X = center.X();
            double Y = center.Y();
            double Z = center.Z();

            if(node->getLeafs() == 0)
            {
                node->add_cell(cell);

                // check for maximal cells in node
                if(node->get_size() >= MAX_CELLS)
                {
                    // check for maximal level
                    if(node->getLevel() < MAX_LEVEL)
                    {
                        // do refinement
                        createLeafNodes(node);
                        node->push_down();
                    }
                    return true;
                }

                return true;
            }

            node->increment();

            double xmid = node->getXCenter();
            double ymid = node->getYCenter();
            double zmid = node->getZCenter();

            if(Z < zmid) // z 0 1 2 3
            {
                if(Y < ymid) // y 0 1
                {
                    if(X < xmid) // x 0
                        node = node->getLeaf(0);
                    else
                        // x 1
                        node = node->getLeaf(1);
                }
                else // y 2 3
                {
                    if(X < xmid) // x 2
                        node = node->getLeaf(2);
                    else // x 3
                        node = node->getLeaf(3);
                }
            }
            else // z 4 5 6 7
            {
                if(Y < ymid) // y 4 5
                {
                    if(X < xmid) // x 4
                        node = node->getLeaf(4);
                    else // x 5
                        node = node->getLeaf(5);
                }
                else // y 6 7
                {
                    if(X < xmid) // x 6
                        node = node->getLeaf(6);
                    else // x 7
                        node = node->getLeaf(7);
                }
            }

            return addCell(node, cell);
        }

        // find connected nodes in the next level
        tree_node * goNextLevelDown(tree_node * node, Vector3D point, uint _level)
        {
            if(node->getLeafs() == 0)
                return node;

            if(node->getLevel() == _level)
                return node;

            double xmid = node->getXCenter();
            double ymid = node->getYCenter();
            double zmid = node->getZCenter();

            if(point.Z() < zmid) // z 0 1 2 3
            {
                if(point.Y() < ymid) // y 0 1
                {
                    if(point.X() < xmid) // x 0
                        node = node->getLeaf(0);
                    else
                        // x 1
                        node = node->getLeaf(1);
                }
                else // y 2 3
                {
                    if(point.X() < xmid) // x 2
                        node = node->getLeaf(2);
                    else // x 3
                        node = node->getLeaf(3);
                }
            }
            else // z 4 5 6 7
            {
                if(point.Y() < ymid) // y 4 5
                {
                    if(point.X() < xmid) // x 4
                        node = node->getLeaf(4);
                    else // x 5
                        node = node->getLeaf(5);
                }
                else // y 6 7
                {
                    if(point.X() < xmid) // x 6
                        node = node->getLeaf(6);
                    else // x 7
                        node = node->getLeaf(7);
                }
            }

            return goNextLevelDown(node, point, _level);
        }

        // final tree cleanup
        void clear(tree_node * node)
        {
            tree_node * leafs = node->getLeafs();

            node->clear();

            if(leafs == 0)
                return;

            for(uint i = 0; i < 8; i++)
                clear(&leafs[i]);

            delete[] leafs;
            node->setLeafs(0);
        }

        // create eight sub-nodes for node pointer
        bool createLeafNodes(tree_node * node)
        {
            uint next_level = 1 + node->getLevel();

            if(MAX_LEVEL < next_level)
                return false;

            max_level = next_level;
            max_nodes += 4;

            double ox = node->getXMin();
            double oy = node->getYMin();
            double oz = node->getZMin();

            double tmp_length = 0.5 * node->getLength();
            tree_node * leafs = new tree_node[8];

            leafs[0].setBranch(node);
            leafs[0].setXMin(ox);
            leafs[0].setYMin(oy);
            leafs[0].setZMin(oz);

            leafs[0].setLength(tmp_length);
            leafs[0].setLevel(next_level);

            leafs[1].setBranch(node);
            leafs[1].setXMin(ox + tmp_length);
            leafs[1].setYMin(oy);
            leafs[1].setZMin(oz);

            leafs[1].setLength(tmp_length);
            leafs[1].setLevel(next_level);

            leafs[2].setBranch(node);
            leafs[2].setXMin(ox);
            leafs[2].setYMin(oy + tmp_length);
            leafs[2].setZMin(oz);

            leafs[2].setLength(tmp_length);
            leafs[2].setLevel(next_level);

            leafs[3].setBranch(node);
            leafs[3].setXMin(ox + tmp_length);
            leafs[3].setYMin(oy + tmp_length);
            leafs[3].setZMin(oz);

            leafs[3].setLength(tmp_length);
            leafs[3].setLevel(next_level);

            leafs[4].setBranch(node);
            leafs[4].setXMin(ox);
            leafs[4].setYMin(oy);
            leafs[4].setZMin(oz + tmp_length);

            leafs[4].setLength(tmp_length);
            leafs[4].setLevel(next_level);

            leafs[5].setBranch(node);
            leafs[5].setXMin(ox + tmp_length);
            leafs[5].setYMin(oy);
            leafs[5].setZMin(oz + tmp_length);

            leafs[5].setLength(tmp_length);
            leafs[5].setLevel(next_level);

            leafs[6].setBranch(node);
            leafs[6].setXMin(ox);
            leafs[6].setYMin(oy + tmp_length);
            leafs[6].setZMin(oz + tmp_length);

            leafs[6].setLength(tmp_length);
            leafs[6].setLevel(next_level);

            leafs[7].setBranch(node);
            leafs[7].setXMin(ox + tmp_length);
            leafs[7].setYMin(oy + tmp_length);
            leafs[7].setZMin(oz + tmp_length);

            leafs[7].setLength(tmp_length);
            leafs[7].setLevel(next_level);

            node->setLeafs(leafs);

            return true;
        }

        double side_length;
        tree_node * root;
        uint max_level;
        uint max_nodes;
    };

    uint min_nrOfNeigbors;
    uint max_nrOfNeigbors;

    uint hull_size;
    h_list * hull_list;
    search_tree * stree;

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

    // obsolete function
    // for testing purposes only
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

    void addGNULines(uint cID, stringstream & str)
    {
        cell_vo * tmp_cell = (cell_vo *)cell_list[cID];
        uint nr_neighbors = tmp_cell->getNrOfNeighbors();

        Vector3D p1 = tmp_cell->getCenter();

        for(uint i = 0; i < nr_neighbors; i++)
        {
            if(isNeigboringVoroCell(tmp_cell, i))
            {
                Vector3D p2 = getNeighborCenter(tmp_cell, i);
                Vector3D tmp_len = p2 - p1;
                str << p1.X() << " " << p1.Y() << " " << p1.Z() << " " << tmp_len.X() << " " << tmp_len.Y()
                    << " " << tmp_len.Z() << "\n";
            }
        }
    }

    bool isNeigboringVoroCell(cell_vo * cell, uint nID)
    {
        int id = cell->getNeighborID(nID);

        if(id > int(max_cells))
            return false;

        return id > -1;
    }

    Vector3D getNeighborCenter(cell_vo * cell, uint nID)
    {
        int id = cell->getNeighborID(nID);
        cell_vo * n_cell = ((cell_vo *)cell_list[id]);
        return n_cell->getCenter();
    }
};
