/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) <year>  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#ifndef HDF5IO_H
#define HDF5IO_H

#include <iostream>
#include <cstdlib>

#include <hdf5.h>


#include "VectorMap.h"




class HDF5IO
{
    hid_t   fid, RawGID, StatGID;
    hid_t   currentSnapshot;
    int nSnapshots;


    hid_t GetH5T(double x) {
        return H5T_NATIVE_DOUBLE;
    }


    hid_t GetH5T(int x) {
        return H5T_NATIVE_INT;
    }

public:
    HDF5IO() : fid(-1) {}
    HDF5IO(string fname) {
        open(fname);
    }
    
    void open(string fname) {
        nSnapshots = 0;
        fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        if (fid < 0) {
            std::cerr << "File " << fname << " could not be created\n.";
            std::exit(EXIT_FAILURE);
        }
        RawGID = H5Gcreate(fid, "Snapshots", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (RawGID < 0) {
            std::cerr << "Group /Snapshots could not be created\n.";
            std::exit(EXIT_FAILURE);
        }
        StatGID = H5Gcreate(fid, "Stats", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (StatGID < 0) {
            std::cerr << "Group /Stats could not be created\n.";
            std::exit(EXIT_FAILURE);
        }
    }

    ~HDF5IO() {
        if (fid < 0) return;
        H5Gclose(StatGID);
        H5Fclose(fid);
    }
    
    hid_t getFID() {
        return fid;
    }

    void newSnapshot() {
        stringstream gname;

        gname << nSnapshots;

        currentSnapshot = H5Gcreate(RawGID, gname.str().c_str(), H5P_DEFAULT,
                                      H5P_DEFAULT, H5P_DEFAULT);
    }

    void closeSnapshot() {
        H5Gclose(currentSnapshot);
        nSnapshots++;
    }

    template <typename T>
    void addSnapshotFields(VectorMap<T>& f, hsize_t N) {

        hid_t data_space, data_set_id, data_type;
        typename VectorMap<T>::iterator i;
        T x;

        hsize_t dims[1] = {N};
        data_space = H5Screate_simple(1, dims, NULL);

        data_type = GetH5T(x);

        for (i = f.begin(); i != f.end(); i++) {

            data_set_id = H5Dcreate(currentSnapshot, i->first.c_str(), data_type,
                                    data_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            H5Dwrite(data_set_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     i->second->memptr());

            H5Dclose(data_set_id);
        }

        H5Sclose(data_space);
    }

    template <typename T>
    void addSnapshotField(const char *name, Col<T> &v, hsize_t N) {
        hid_t data_space, data_set_id, data_type;
        T x;

        data_space = H5Screate_simple(1, &N, NULL);
        data_type = GetH5T(x);

        data_set_id = H5Dcreate(currentSnapshot, name, data_type, data_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(data_set_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 v.memptr());

        H5Dclose(data_set_id);

        H5Sclose(data_space);
    }

    void addStatsField(const char *name, mat &m) {
        hid_t data_space, data_set_id;
        hsize_t dims[2];

        dims[1] = m.n_rows;
        dims[0] = m.n_cols;

        data_space = H5Screate_simple(2, dims, NULL);

        data_set_id = H5Dcreate(StatGID, name, H5T_NATIVE_DOUBLE, data_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(data_set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 m.memptr());

        H5Dclose(data_set_id);

        H5Sclose(data_space);
    }

    void addStatsField(const char *name, cube &m) {
        hid_t data_space, data_set_id;
        hsize_t dims[3];

        dims[2] = m.n_rows;
        dims[1] = m.n_cols;
        dims[0] = m.n_slices;

        data_space = H5Screate_simple(3, dims, NULL);

        data_set_id = H5Dcreate(StatGID, name, H5T_NATIVE_DOUBLE, data_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(data_set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 m.memptr());

        H5Dclose(data_set_id);

        H5Sclose(data_space);
    }
    
        void addStatsCubeSeries(const char *name, cube &m, hsize_t nCubes) {
        hid_t data_space, data_set_id;
        hsize_t dims[4];

        dims[3] = m.n_rows;
        dims[2] = m.n_cols;
        dims[1] = m.n_slices/nCubes;
        dims[0] = nCubes;

        data_space = H5Screate_simple(4, dims, NULL);

        data_set_id = H5Dcreate(StatGID, name, H5T_NATIVE_DOUBLE, data_space,
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite(data_set_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 m.memptr());

        H5Dclose(data_set_id);

        H5Sclose(data_space);
    }
};

#endif // HDF5IO_H
