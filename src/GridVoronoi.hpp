/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_VORONOI_H
#define CGRID_VORONOI_H

#include "GridBasic.hpp"
#include "Vector3D.hpp"
#include "CellVoronoi.hpp"
#include "Typedefs.hpp"

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

        /*min_delta = 0;
        max_delta = 0;*/

        min_mach = 1e300;
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

        nrOfPlotPoints = 1000;
        nrOfPlotVectors = 1000;
        maxPlotLines = 3;

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
        //plt_delta = false;
        plt_larm = false;
        //plt_mach = false;
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

        // cout << CLR_LINE << flush;
    }

    bool isInside(const Vector3D & pos) const;

    bool writePlotFiles(string path, parameters & param);

    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    Vector3D getCenter(const cell_basic & cell) const;

    Vector3D getCenter(uint id) const;

    bool next(photon_package * pp);

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy);

    bool createCellList();

    double getVolume(const cell_basic & cell) const;

    bool positionPhotonInGrid(photon_package * pp);

    // for debugging only
    bool positionPhotonInGridTest(photon_package * pp);
    bool createArtificialGrid(string path);

    bool saveBinaryGridFile(string filename);

    bool loadGridFromBinaryFile(parameters & param, uint _data_len);
    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinaryFile(parameters & param);

    // final cleanup
    void clear();

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
            bool nodeIntersection(Vector3D point, double _sq_distance);

            bool isPointInNode(Vector3D point);

            // keep track of the amount of cells in each level
            void increment();

            // add cell pointer to linked list
            void add_cell(cell_vo * cell);

            // final cleanup
            void clear();

            // empty node and refine tree if to many cell populate the node
            void push_down();

            // find closest cell for given point
            cell_vo * findClosestCell(Vector3D point, double & _min_distance, tree_node *& f_node);

            uint get_size();

            bool is_emty();

            double getXMin();

            double getYMin();

            double getZMin();

            double getXMax();

            double getYMax();

            double getZMax();

            double getXCenter();

            double getYCenter();

            double getZCenter();

            void setXMin(double x);

            void setYMin(double y);

            void setZMin(double y);

            double getLength() const;

            void setLength(double l);

            void setLevel(uint l);

            list_element * get_first();

            list_element * get_last();

            uint getLevel();

            tree_node * getLeafs();

            tree_node * getLeaf(uint index);

            tree_node * getBranch();

            void setLeafs(tree_node * l);

            void setBranch(tree_node * b);

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

        bool addCell(cell_vo * cell);

        // init. tree and first level refinement
        void initTree(double _side_length);

        uint getMaxLevel();

        uint getMaxNodes();

        // find closest cell for given point in the entire tree
        cell_vo * findClosestCell(Vector3D point, cell_basic ** cell_list);

    private:
        // check neighboring nodes for shortest distance
        cell_vo * checkNeighboringNodes(tree_node * node,
                                        tree_node * f_node,
                                        tree_node * p_node,
                                        Vector3D point,
                                        double distance,
                                        double & min_distance);

        tree_node * findMatchingNode(Vector3D point);

        // go to maximal level for starting the point search
        tree_node * findMatchingNode(Vector3D point, uint _level);

        // add cell to tree
        bool addCell(tree_node * node, cell_vo * cell);

        // find connected nodes in the next level
        tree_node * goNextLevelDown(tree_node * node, Vector3D point, uint _level);

        // final tree cleanup
        void clear(tree_node * node);

        // create eight sub-nodes for node pointer
        bool createLeafNodes(tree_node * node);

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

    double abs_min(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

    // obsolete function (for testing purposes only)
    bool isHullPoint(uint id);

    void addPlotLines(uint cID, stringstream & str);

    bool isNeigboringVoroCell(cell_vo * cell, uint nID);

    Vector3D getNeighborCenter(cell_vo * cell, uint nID);
};

#endif /* CGRID_VORONOI_H */
