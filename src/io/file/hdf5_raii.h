/***
Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced
under U.S. Government contract DE-AC02-06CH11357 for Argonne National
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
modified to produce derivative works, such modified software should
be clearly marked, so as not to confuse it with the version available
from ANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
***/

/// Initial Author <2026>: Arthur Glowacki

#ifndef HDF5_RAII_H
#define HDF5_RAII_H

//-----------------------------------------------------------------------------
//
// Scope-bound ownership for HDF5 identifiers.
//
// The older pattern in hdf5_io.h tracks handles in a
// std::stack<std::pair<hid_t, H5_OBJECTS>> and relies on every exit path
// remembering to drain it. Because the handle and the tag naming its closer
// are stored as two independent values, a copy/paste slip can pair a handle
// with the wrong closer. That is not hypothetical: commit fc68201 fixed a
// Linux segfault caused by pushing a dataset id tagged H5O_DATATYPE, so
// H5Tclose() ran against a dataset.
//
// Here the closer is part of the type, so it cannot disagree with the handle,
// and release happens at scope exit on every path including exceptions.
// Obtain handles through the h5_* factories below and bind them with `auto`;
// then the correct closer is selected by construction rather than by hand.
//
//-----------------------------------------------------------------------------

#include <string>

#include "hdf5.h"
#include "core/defines.h"

namespace io
{
namespace file
{

//-----------------------------------------------------------------------------

/// Move-only owner of an hid_t. CloseFn is baked into the type, so the
/// release call can never be mismatched to the kind of handle held.
template<herr_t (*CloseFn)(hid_t)>
class H5Handle
{
public:

    H5Handle() noexcept = default;

    explicit H5Handle(hid_t id) noexcept : _id(id) {}

    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;

    H5Handle(H5Handle&& other) noexcept : _id(other._id)
    {
        other._id = -1;
    }

    H5Handle& operator=(H5Handle&& other) noexcept
    {
        if (this != &other)
        {
            close();
            _id = other._id;
            other._id = -1;
        }
        return *this;
    }

    ~H5Handle()
    {
        close();
    }

    /// Release early. Safe to call more than once; the handle is cleared first
    /// so the destructor will not close it again.
    void close() noexcept
    {
        if (_id > -1)
        {
            CloseFn(_id);
            _id = -1;
        }
    }

    /// Implicit so an owner can be passed straight to the HDF5 C API.
    operator hid_t() const noexcept { return _id; }

    hid_t get() const noexcept { return _id; }

    bool valid() const noexcept { return _id > -1; }

    explicit operator bool() const noexcept { return _id > -1; }

    /// Hand ownership to the caller; this object no longer closes the handle.
    hid_t release() noexcept
    {
        hid_t tmp = _id;
        _id = -1;
        return tmp;
    }

    /// Take ownership of a new handle, closing whatever was held before.
    void reset(hid_t id = -1) noexcept
    {
        if (id != _id)
        {
            close();
            _id = id;
        }
    }

private:

    hid_t _id = -1;
};

//-----------------------------------------------------------------------------

using H5File    = H5Handle<&H5Fclose>;
using H5Group   = H5Handle<&H5Gclose>;
using H5Dataset = H5Handle<&H5Dclose>;
using H5Space   = H5Handle<&H5Sclose>;
using H5Type    = H5Handle<&H5Tclose>;
using H5Attr    = H5Handle<&H5Aclose>;
using H5Plist   = H5Handle<&H5Pclose>;

//-----------------------------------------------------------------------------
// Factories. Each returns the owner type matching the handle it produced, so
// `auto h = h5_get_type(dset);` cannot pick up the wrong closer. A failed open
// yields an invalid handle rather than throwing; check with valid().
//-----------------------------------------------------------------------------

inline H5File h5_open_file(const std::string& path, unsigned flags = H5F_ACC_RDONLY, bool log_error = true)
{
    hid_t id = H5Fopen(path.c_str(), flags, H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to open file " << path << "\n";
    }
    return H5File(id);
}

//-----------------------------------------------------------------------------

inline H5File h5_create_file(const std::string& path, bool log_error = true)
{
    hid_t id = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to create file " << path << "\n";
    }
    return H5File(id);
}

//-----------------------------------------------------------------------------

inline H5Group h5_open_group(hid_t parent, const std::string& name, bool log_error = true)
{
    hid_t id = H5Gopen(parent, name.c_str(), H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to open group " << name << "\n";
    }
    return H5Group(id);
}

//-----------------------------------------------------------------------------

/// Open name under parent, creating it when absent.
inline H5Group h5_open_or_create_group(hid_t parent, const std::string& name, bool log_error = true)
{
    hid_t id = H5Gopen(parent, name.c_str(), H5P_DEFAULT);
    if (id < 0)
    {
        id = H5Gcreate(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (id < 0 && log_error)
    {
        logE << "creating group " << name << "\n";
    }
    return H5Group(id);
}

//-----------------------------------------------------------------------------

inline H5Dataset h5_open_dataset(hid_t parent, const std::string& name, bool log_error = true)
{
    hid_t id = H5Dopen2(parent, name.c_str(), H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to open dataset " << name << "\n";
    }
    return H5Dataset(id);
}

//-----------------------------------------------------------------------------

inline H5Dataset h5_create_dataset(hid_t parent, const std::string& name, hid_t type_id, hid_t space_id, hid_t dcpl_id = H5P_DEFAULT, bool log_error = true)
{
    hid_t id = H5Dcreate(parent, name.c_str(), type_id, space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to create dataset " << name << "\n";
    }
    return H5Dataset(id);
}

//-----------------------------------------------------------------------------

inline H5Space h5_create_space(int rank, const hsize_t* dims, const hsize_t* max_dims = nullptr)
{
    return H5Space(H5Screate_simple(rank, dims, max_dims));
}

//-----------------------------------------------------------------------------

inline H5Space h5_get_space(hid_t dset_id)
{
    return H5Space(H5Dget_space(dset_id));
}

//-----------------------------------------------------------------------------

inline H5Space h5_get_attr_space(hid_t attr_id)
{
    return H5Space(H5Aget_space(attr_id));
}

//-----------------------------------------------------------------------------

inline H5Type h5_get_type(hid_t dset_id)
{
    return H5Type(H5Dget_type(dset_id));
}

//-----------------------------------------------------------------------------

inline H5Type h5_get_attr_type(hid_t attr_id)
{
    return H5Type(H5Aget_type(attr_id));
}

//-----------------------------------------------------------------------------

inline H5Type h5_copy_type(hid_t type_id)
{
    return H5Type(H5Tcopy(type_id));
}

//-----------------------------------------------------------------------------

/// H5Tcopy followed by H5Tset_size, the fixed-length-string idiom used
/// throughout hdf5_io.h.
inline H5Type h5_string_type(hid_t base_type_id, size_t size)
{
    H5Type t(H5Tcopy(base_type_id));
    if (t.valid())
    {
        H5Tset_size(t, size);
    }
    return t;
}

//-----------------------------------------------------------------------------

inline H5Type h5_native_type(hid_t type_id, H5T_direction_t direction = H5T_DIR_ASCEND)
{
    return H5Type(H5Tget_native_type(type_id, direction));
}

//-----------------------------------------------------------------------------

inline H5Attr h5_open_attr(hid_t obj_id, const std::string& name, bool log_error = true)
{
    hid_t id = H5Aopen(obj_id, name.c_str(), H5P_DEFAULT);
    if (id < 0 && log_error)
    {
        logW << "Failed to open attribute " << name << "\n";
    }
    return H5Attr(id);
}

//-----------------------------------------------------------------------------

inline H5Plist h5_create_plist(hid_t cls_id)
{
    return H5Plist(H5Pcreate(cls_id));
}

//-----------------------------------------------------------------------------

inline H5Plist h5_get_create_plist(hid_t dset_id)
{
    return H5Plist(H5Dget_create_plist(dset_id));
}

//-----------------------------------------------------------------------------

} // end namespace file
} // end namespace io

#endif // HDF5_RAII_H
