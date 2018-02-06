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


#include <valarray>
#include "core/defines.h"
#include "data_struct/spectra.h"

#define vtkDomainsChemistry_AUTOINIT 1(vtkDomainsChemistryOpenGL2)
#define vtkRenderingContext2D_AUTOINIT 1(vtkRenderingContextOpenGL2)
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

///VTK
#include <vtkVersion.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkAxis.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkGraphicsFactory.h>

namespace visual
{

void PlotSpectras(data_struct::ArrayXr spectra1, data_struct::ArrayXr spectra2, bool log_them=false);
void SavePlotSpectras(std::string path, data_struct::ArrayXr xAxis, data_struct::ArrayXr spectra1, data_struct::ArrayXr spectra2, bool log_them);

}
