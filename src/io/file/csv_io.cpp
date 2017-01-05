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



#include "csv_io.h"

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <iostream>

namespace io
{
namespace file
{

CSV_IO::CSV_IO()
{

}

CSV_IO::~CSV_IO()
{

}

bool CSV_IO::save_fit_parameters(std::string fullpath, data_struct::xrf::Fit_Parameters fit_params)
{
    std::ofstream file_stream(fullpath);
    if(file_stream.is_open())
    {
        file_stream<<"Optimized Fit Parameters"<<std::endl;
        file_stream<<"Name,Value"<<std::endl;

        for (auto itr = fit_params.begin(); itr != fit_params.end(); itr++)
        {
            file_stream<<itr->first<<","<<itr->second.value<<std::endl;
        }
        file_stream.close();
    }
    else
    {
        return false;
    }
    return true;
}

bool load_element_info_from_csv(std::string filename, data_struct::xrf::Element_Info_Map *element_map)
{
    std::ifstream file_stream(filename);
    try
    {
        std::string value;
        for (std::string line; std::getline(file_stream, line); )
        {
            std::stringstream strstream(line);
            std::getline(strstream, value, ',');
            /*
            if (value == "version:")
            {
                std::getline(strstream, value, ',');
                if (value != "1.2")
                    std::cout<<"Warning: non expected version. Loader if for 1.2"<<std::endl;
            }*/
            //if( std::stoi(value) > 0)
            if (value[0] >= 48 && value[0] <= 57) // 0 - 9
            {

                //std::cout<< "value = "<< value<<std::endl;
                data_struct::xrf::Element_Info* element = nullptr;
                int element_number = std::stoi(value);
                element = element_map->get_element(element_number);
                std::string el_name;
                std::getline(strstream, el_name, ',');
                if (element == nullptr)
                {
                     element = new data_struct::xrf::Element_Info();
                     element->number = element_number;
                     element->name = el_name;
                     element_map->add_element(element);
                }
                else
                {
                    element->number = element_number;
                    element->name = el_name;
                }
                //std::cout<< element->number << " : "<< element->name <<std::endl;
                std::getline(strstream, value, ',');
                element->xrf["ka1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ka2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["kb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["kb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["la1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["la2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ll"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ln"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ma1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ma2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["mb"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["mg"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["k"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["m"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ka1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ka2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["kb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["kb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["la1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["la2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ll"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ln"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ma1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ma2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["mb"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["mg"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->density = std::stof(value);
                std::getline(strstream, value, ',');
                element->mass = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["K"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["L1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["L2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["L3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["M1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["N1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N5"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N6"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N7"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["O1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["P1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["P2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["P3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["K"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["L1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["L2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["L3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["M1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["N1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["O1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["O2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["O3"] = std::stof(value);

                //element_information.emplace(std::make_pair(element->name, element));
            }
        }
    }
    catch(std::exception e)
    {
        if (file_stream.eof() == 0 && (file_stream.bad() || file_stream.fail()) )
        {
            std::cerr << "ios Exception happened: " << e.what() << "\n"
                << "Error bits are: "
                << "\nfailbit: " << file_stream.fail()
                << "\neofbit: " << file_stream.eof()
                << "\nbadbit: " << file_stream.bad() << std::endl;
        }
        return false;
    }

    return true;
}


} //end namespace file
}// end namespace io
