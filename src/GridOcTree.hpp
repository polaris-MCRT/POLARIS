/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CGRID_OCTREE_H
#define CGRID_OCTREE_H

#include "GridBasic.hpp"
#include "RandomGenerator.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"
#include "CellOcTree.hpp"
#include "Photon.hpp"

class CGridOcTree : public CGridBasic
{
public:
    CGridOcTree(void)
    {
        // octree internals
        cell_oc_root = 0;
        cell_oc_pos  = 0;
        rec_counter  = 0;
        max_level    = 0;
    }

    ~CGridOcTree()
    {
        if(cell_oc_root == 0)
            return;

        clear(cell_oc_root);

        if(cell_list != 0)
        {
            delete[] cell_list;
            cell_list = 0;
        }

        cout << CLR_LINE;
    }

    

    void goToRoot();

    bool nextLowLevelCell();
    bool nextLowLevelCell(cell_basic * cell);
    // end   IO functions
    bool goToNextCellBorder(photon_package * pp);
    bool updateShortestDistance(photon_package * pp);

    // void deleteSubCells(slist & source);

    bool reduceBinaryFile(string in_filename, string out_filename, uint tr_level);
    bool reduceLevelOfBinaryFile(cell_oc * cell, uint tr_level);

    Vector3D getCenter(const cell_basic & cell) const;

    Vector3D getMidplaneCenter(cell_basic * cell);

    bool createCellList();

    bool findMatchingCell(photon_package * pp);

    bool next(photon_package * pp);

    bool findStartingPoint(photon_package * pp);

    void getLengths(uint bins, double & step_xy, double & off_xy);

    double getVolume(const cell_basic & cell) const;

    bool positionPhotonInGrid(photon_package * pp);

    const cell_oc * getTopLevelCell() const;

    const cell_oc * getCurrentCell() const;

    void printParameters();

    bool createArtificialGrid(string path);

    void createNextLevel(cell_oc * cell);

    bool saveBinaryGridFile(string filename);

    bool loadGridFromBinaryFile(parameters & param, uint data_len);

    bool saveBinaryGridFile(string filename, ushort id, ushort data_size);

    bool loadGridFromBinaryFile(parameters & param);

    void nextBinaryDataCell(ofstream & file_stream, cell_oc * cell, uint data_size);

    void clear();

    void goNextLevelDown(photon_package * pp);

    bool createTree(cell_oc * parent,
                    double _x_min,
                    double _y_min,
                    double _z_min,
                    double _length,
                    uint _level);

    bool initiateTreeFromFile(uint _nx,
                              uint _max_level,
                              double _fa,
                              double _length,
                              string str_dens,
                              string str_temp,
                              string str_magx,
                              string str_magy,
                              string str_magz);

private:
    void clear(cell_oc * cell);
    cell_oc * cell_oc_root;
    cell_oc * cell_oc_pos;

    Matrix3D datdens;
    Matrix3D dattemp;
    Matrix3D datmx;
    Matrix3D datmy;
    Matrix3D datmz;

    uint nx;
    // uint ny;
    // uint nz;
    double f_min;
    double f_max;
    double factor;
    uint treelevel_counter;
    uint tagged_cells;

    uint rec_counter;
    double max_level;

    void createBoundingCell();
    void createBoundingCell(cell_oc * cell);

    void goNextLevelUp(photon_package * pp);

    bool isInside(const Vector3D & pos, const cell_basic & _cell) const;

    bool isInside(const Vector3D & pos) const;

    // bool isInside(photon_package * pp, Vector3D & pos);

    void setRndPositionInCell(photon_package * pp, CRandomGenerator * rand_gen);
};

#endif /* CGRID_OCTREE_H */
