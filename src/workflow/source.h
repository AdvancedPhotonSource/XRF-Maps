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

/// Initial Author <2017>: Arthur Glowacki



#ifndef Source_H
#define Source_H

#include "core/defines.h"
#include <functional>
#include "workflow/distributor.h"

namespace workflow
{

//-----------------------------------------------------------------------------
template<typename T_OUT>
class DLL_EXPORT Source
{

   typedef std::function<void (T_OUT)> Callback_Func_Def;
 //  typedef std::function<void ( void*, Callback_Func_Def )> Source_Func_Def;

public:

    Source()
    {
        _output_callback_func = nullptr;
    }

    ~Source()
    {

    }

    template<typename _T>
    void connect(Distributor<T_OUT, _T> *distributor)
    {
        _output_callback_func = std::bind(&Distributor<T_OUT, _T>::distribute, distributor, std::placeholders::_1);
    }

/*
    void connect( Callback_Func_Def out_callback_func)
    {
        _output_callback_func = out_callback_func;
    }

    void set_function(Source_Func_Def func)
    {
        _prod_func = func;
    }
*/
    virtual void run() = 0;

/*
    template<class... Args>
    virtual void run(Args&&... args)
    {

        auto task = std::make_shared< std::packaged_task<void> >(
                std::bind(std::forward<Source_Func_Def>(_prod_func), std::forward<Args>(args)...)
            );

        task();
    }
*/
protected:

    Callback_Func_Def _output_callback_func;

};

} //namespace workflow

#endif // Source_H
