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


#include <hdf5.h>

#include "VectorMap.h"


class HDF5IO
{
    hid_t   fid;
    hid_t   current_group_id;
    int ngroups;


    hid_t GetH5T(double x) {
        return H5T_NATIVE_DOUBLE;
    }


    hid_t GetH5T(int x) {
        return H5T_NATIVE_INT;
    }

public:
    HDF5IO(string fname) : ngroups(0) {
        fid = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    ~HDF5IO() {
        H5Fclose(fid);
    }

    void newGroup() {
        stringstream gname;

        ngroups++;

        gname << ngroups;

        current_group_id = H5Gcreate(fid, gname.str().c_str(), H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
    }

    void closeGroup() {
        H5Gclose(current_group_id);
    }

    template <typename T>
    void add(VectorMap<T>& f, hsize_t N) {

        hid_t data_space, data_set_id, data_type;
        typename VectorMap<T>::iterator i;
        T x;

        hsize_t dims[1] = {N};
        data_space = H5Screate_simple(1, dims, NULL);

        data_type = GetH5T(x);

        for (i = f.begin(); i != f.end(); i++) {

            data_set_id = H5Dcreate(current_group_id, i->first.c_str(), data_type,
                                    data_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            H5Dwrite(data_set_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     i->second->memptr());

            H5Dclose(data_set_id);
        }

        H5Sclose(data_space);
    }

    template <typename T>
    void add(const char *name, Col<T> &v, hsize_t N) {
        hid_t data_space, data_set_id, data_type;
        T x;

        data_space = H5Screate_simple(1, &N, NULL);
        data_type = GetH5T(x);

	data_set_id = H5Dcreate(current_group_id, name, data_type, data_space,
				H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	H5Dwrite(data_set_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     v.memptr());

	H5Dclose(data_set_id);

        H5Sclose(data_space);
    }

};

#endif // HDF5IO_H
