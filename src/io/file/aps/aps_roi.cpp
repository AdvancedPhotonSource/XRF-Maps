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



#include "aps_roi.h"
#include <fstream>

namespace io
{
namespace file
{
namespace aps
{

union Num
{
    char buffer[4];
    unsigned int num;
} num;

void swapChars(char* pChar1, char* pChar2)
{
    char temp = *pChar1;
    *pChar1 = *pChar2;
    *pChar2 = temp;
}

int swapOrder(Num num)
{
    swapChars(&num.buffer[0], &num.buffer[3]);
    swapChars(&num.buffer[1], &num.buffer[2]);

    return num.num;
}

bool load_v9_rois(std::string path, std::map<int, std::vector<int_point>>& rois)
{
    std::ifstream fileStream(path, std::ios::in | std::ios::binary);

    logI << "Loading:  " << path << "\n";
    Num val;
    unsigned int width;
    unsigned int height;
    unsigned int mask;
    std::vector<unsigned int> myData;
    if (fileStream.is_open())
    {
        fileStream.read((char*)&val, sizeof(val));
        width = swapOrder(val);
        fileStream.read((char*)&val, sizeof(val));
        height = swapOrder(val);

        myData.resize(width * height);
        fileStream.read((char*)myData.data(), (width * height) * sizeof(unsigned int));
        fileStream.close();

        unsigned long i = 0;
        for (unsigned int y = 0; y < height; y++)
        {
            for (unsigned int x = 0; x < width; x++)
            {
                if (myData[i])
                {
                    val.num = myData[i];
                    mask = swapOrder(val);
                    for (int idx = 1; idx < 11; idx++)
                    {
                        if ((idx & mask) == idx)
                        {
                            rois[idx-1].push_back(int_point(x, y));
                        }
                    }
                    
                }
                i++;
            }
        }
        fileStream.close();
    }
    else
    {
        return false;
    }


    return true;

}





} //end namespace aps
} //end namespace file
}// end namespace io
