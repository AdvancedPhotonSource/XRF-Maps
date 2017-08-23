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

/// Initial Author <2016>: Arthur Glowacki


#ifndef BASE_DATASET_H
#define BASE_DATASET_H

#include "core/defines.h"



//FOR TESTING NOT USED

//#include <array>
//#include <vector>
//#include <algorithm>
//#include <memory>
/**
 * @brief The base_dataset class:
 */

/*
template <class T>
class multi_vec
{
public:
    using param=std::vector<size_t>;

    explicit multi_vec(const param& dimensions)
        : dim{dimensions}, prod {1}
    {
        std::for_each(dim.begin(), dim.end(), [this] (std::size_t val)
        {
            mult.emplace_back(prod);
            prod *= val;
        } );
        ptr.reset(new T[prod]);
    }
    std::size_t capacity() const { return prod; }

    // undefined if elements in lookup != elemenets in dim
    // undefined if any element in lookup
        // is greater than or equal to corresponding dim element
    T& operator[](const param& lookup)
    {
        return ptr[get_offset(lookup)];
    }
    const T operator[](const param& lookup) const
    {
        return ptr[get_offset(lookup)];
    }
private:
    std::size_t get_offset(const param& lookup) const
    {
        std::size_t offset=0;
        auto mit=mult.begin();
        std::for_each(lookup.begin(), lookup.end(), [&offset, &mit] (std::size_t val)
        {
            offset+=*mit * val;
           ++mit;
        } );
        return offset;
    }
    param dim;
    param mult;
    std::size_t prod;
    std::unique_ptr<T[]> ptr;
};
*/
class Base_Dataset
{
public:

    //enum DATA_TYPE {_int8, _int16, _int32, _int64, _float32, _float64 };

    /**
     * @brief Base_Dataset
     */
    Base_Dataset() {}

    /**
     * @brief Base_Dataset: destructor
     */
    ~Base_Dataset() {}

    /**
     * @brief get_data_ptr: get a pointer to the data object
     * @return
     */
    ///NOTE: change return from void* to ndarray*
    virtual real_t* get_buffer_ptr() = 0;

    /**
     * @brief get_data_type
     * @return
     */
    //virtual DATA_TYPE get_data_type()=0;

    /**
     * @brief alloc: allocate the dataset
     * @param rank
     * @param dims
     * @param dtype
     */
    //virtual void alloc(int rank, unsigned long long* dims)=0;

    virtual void alloc(unsigned short rank, int *dims)=0;

    /**
     * @brief rank: get the dataset rank
     * @return
     */
    virtual int rank() = 0;

    /**
     * @brief dims: get the dataset dims
     * @return
     */
    virtual void dims(int &dims_out) = 0;

};

#endif // BASE_DATASET_H
