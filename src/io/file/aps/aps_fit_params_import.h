/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

#ifndef APS_FIT_PARAMS_IMPORT_H
#define APS_FIT_PARAMS_IMPORT_H

#include "defines.h"
#include "fit_parameters.h"
#include "fit_element_map.h"
#include "element_info.h"

namespace io
{
namespace file
{
namespace aps
{

class DLL_EXPORT APS_Fit_Params_Import
{
public:

    APS_Fit_Params_Import();

    ~APS_Fit_Params_Import();

    bool load(std::string path,
              data_struct::xrf::Element_Info_Map *element_info_map,
              data_struct::xrf::Fit_Parameters* out_fit_params,
              std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*>* out_elements_to_fit,
              std::unordered_map<std::string, std::string>* out_values);

private:


};

}// end namespace aps
}// end namespace file
}// end namespace io

#endif // APS_FIT_PARAMS_IMPORT
