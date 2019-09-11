#include "Typedefs.h"
#include "Vector.h"

#ifndef CELL
#define CELL

class cell_basic
{
  public:
    cell_basic()
    {
        data = 0;
    }

    virtual ~cell_basic()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    bool isValid()
    {
        return data != 0;
    }

    double getData(uint i) const
    {
        if(data == 0)
            return 0;

        if(i == MAX_UINT)
            return 0;

        if(i > 4138637280)
            return 0;

        return data[i];
    }

    void setData(uint i, double d)
    {
#pragma omp atomic write
        data[i] = d;
    }

    void updateData(uint i, double d)
    {
#pragma omp atomic update
        data[i] += d;
    }

    void convertData(uint i, double c)
    {
#pragma omp atomic update
        data[i] *= c;
    }

    void resize(uint size)
    {
        if(data != 0)
            delete[] data;

        data = new double[size];
        for(uint i = 0; i < size; i++)
            data[i] = 0;
    }

    void setID(uint _id)
    {
        id = _id;
    }

    uint getID()
    {
        return id;
    }

    virtual ulong getUniqueID()
    {
        return ulong(id);
    }

    void updateID(uint _id)
    {
        id += _id;
    }

  protected:
    double * data;
    uint id;
};

class cell_oc : public cell_basic
{
  public:
    cell_oc()
    {
        parent = 0;
        children = 0;
        level = 0;
        id = 0;
        x_min = 0;
        y_min = 0;
        z_min = 0;
        length = 0;
        data = 0;
    }

    ~cell_oc()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setXmin(double _x_min)
    {
        x_min = _x_min;
    }

    void setYmin(double _y_min)
    {
        y_min = _y_min;
    }

    void setZmin(double _z_min)
    {
        z_min = _z_min;
    }

    void setLength(double _length)
    {
        length = _length;
    }

    void setLevel(uchar _level)
    {
        level = _level;
    }

    void setParent(cell_oc * _parent)
    {
        parent = _parent;
    }

    void setChildren(cell_oc * _children)
    {
        children = _children;
    }

    double getXmin() const
    {
        return x_min;
    }

    double getYmin() const
    {
        return y_min;
    }

    double getZmin() const
    {
        return z_min;
    }

    double getXmax() const
    {
        return x_min + length;
    }

    double getYmax() const
    {
        return y_min + length;
    }

    double getZmax() const
    {
        return z_min + length;
    }

    double getLength() const
    {
        return length;
    }

    uchar getLevel()
    {
        return level;
    }

    ulong getUniqueID()
    {
        uint parent_id = parent->getID();
        ulong res = 8 ^ level + 8 * parent_id + id;
        return res;
    }

    cell_oc * getParent()
    {
        return parent;
    }

    cell_oc * getChildren()
    {
        return children;
    }

    cell_oc * getChild(uint i)
    {
        return &children[i];
    }

  private:
    double x_min, y_min, z_min, length;
    cell_oc * children;
    cell_oc * parent;
    uchar level;
};

class cell_vo : public cell_basic
{
  public:
    cell_vo()
    {
        id = 0;
        data = 0;
        nr_neighbors = 0;
        neighbors = 0;
        volume = -1;
    }

    ~cell_vo()
    {
        if(data != 0)
            delete[] data;

        data = 0;

        if(neighbors != 0)
            delete[] neighbors;

        neighbors = 0;
    }

    void initNeighbors(short nr)
    {
        nr_neighbors = nr;
        neighbors = new int[nr];

        for(ushort i = 0; i < nr; i++)
            neighbors[i] = 0;
    }

    void setNeighbor(uint pos, int id)
    {
        if(pos > 500)
            return;

        neighbors[pos] = id;
    }

    void setCenter(double cx, double cy, double cz)
    {
        center.set(cx, cy, cz);
    };

    void setVolume(double v)
    {
        volume = v;
    };

    Vector3D getCenter()
    {
        return center;
    };

    double getX()
    {
        return center.X();
    }

    double getY()
    {
        return center.Y();
    }

    double getZ()
    {
        return center.Z();
    }

    double getVolume() const
    {
        return volume;
    };

    int getNeighborID(uint pos)
    {
        return neighbors[pos];
    }

    ushort getNrOfNeighbors()
    {
        return nr_neighbors;
    }

  private:
    Vector3D center;
    ushort nr_neighbors;
    int * neighbors;
    double volume;
};

class cell_sp : public cell_basic
{
  public:
    cell_sp()
    {
        rID = 0;
        phID = 0;
        thID = 0;
        data = 0;
        id = 0;
    }

    cell_sp(uint _rID, uint _phID, uint _thID)
    {
        rID = _rID;
        phID = _phID;
        thID = _thID;
        data = 0;
        id = 0;
    }

    ~cell_sp()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id)
    {
        rID = id;
    }

    void setPhID(uint id)
    {
        phID = id;
    }

    void setThID(uint id)
    {
        thID = id;
    }

    uint getRID() const
    {
        return rID;
    }

    uint getPhID() const
    {
        return phID;
    }

    uint getThID() const
    {
        return thID;
    }

  private:
    uint rID, phID, thID;
};

class cell_cyl : public cell_basic
{
  public:
    cell_cyl()
    {
        rID = 0;
        phID = 0;
        zID = 0;
        data = 0;
        id = 0;
    }

    ~cell_cyl()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id)
    {
        rID = id;
    }

    void setPhID(uint id)
    {
        phID = id;
    }

    void setZID(uint id)
    {
        zID = id;
    }

    uint getRID() const
    {
        return rID;
    }

    uint getPhID() const
    {
        return phID;
    }

    uint getZID() const
    {
        return zID;
    }

  private:
    uint rID, phID, zID;
};
#endif