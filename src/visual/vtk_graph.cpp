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

#include "vtk_graph.h"
#include <algorithm>

namespace visual
{


void PlotSpectras(data_struct::xrf::ArrayXr spectra1, data_struct::xrf::ArrayXr spectra2, bool log_them)
{
    // Create a table with some points in it
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
    arrC->SetName("Integrated Spectra");
    table->AddColumn(arrC);

    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
    arrS->SetName("Spectra Model");
    table->AddColumn(arrS);

    if(log_them)
    {
        spectra1 = spectra1.log10();
        spectra2 = spectra2.log10();
		spectra1 = spectra1.unaryExpr([](real_t v) { return std::max(v, (real_t)0.0); });
		spectra2 = spectra2.unaryExpr([](real_t v) { return std::max(v, (real_t)0.0); });
    }


    // Fill in the table with some example values
    int numPoints = spectra2.size();
    table->SetNumberOfRows(numPoints);
    real_t en = 0.0;
    for (int i = 0; i < numPoints; ++i)
    {
        table->SetValue(i, 0, en);
        table->SetValue(i, 1, spectra1[i]);
        table->SetValue(i, 2, spectra2[i]);

        //table->SetValue(i, 1, spectra1->buffer()->operator [](i));
        //table->SetValue(i, 2, spectra2->buffer()->operator [](i));
        en += 1.0;
    }

    // Set up the view
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Add multiple line plots, setting the colors etc
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 2);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(1.0);

    // For dotted line, the line type can be from 2 to 5 for different dash/dot
    // patterns (see enum in vtkPen containing DASH_LINE, value 2):
    #ifndef WIN32
    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
    #endif
    // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
    //  machines with built-in Intel HD graphics card...)

    //view->GetRenderWindow()->SetMultiSamples(0);

    // Start interactor
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();

}


void SavePlotSpectras(std::string path, data_struct::xrf::ArrayXr spectra1, data_struct::xrf::ArrayXr spectra2, bool log_them)
{
    // Setup offscreen rendering
    vtkSmartPointer<vtkGraphicsFactory> graphics_factory = vtkSmartPointer<vtkGraphicsFactory>::New();
    graphics_factory->SetOffScreenOnlyMode( 1);
    graphics_factory->SetUseMesaClasses( 1 );

//    vtkSmartPointer<vtkImagingFactory> imaging_factory = vtkSmartPointer<vtkImagingFactory>::New();
//    imaging_factory->SetUseMesaClasses( 1 );


    // Create a table with some points in it
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
    arrC->SetName("Integrated Spectra");
    table->AddColumn(arrC);

    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
    arrS->SetName("Spectra Model");
    table->AddColumn(arrS);

    if(log_them)
    {
		spectra1 = spectra1.log10();
		spectra2 = spectra2.log10();
		spectra1 = spectra1.unaryExpr([](real_t v) { return std::max(v, (real_t)0.0); });
		spectra2 = spectra2.unaryExpr([](real_t v) { return std::max(v, (real_t)0.0); });
    }


    // Fill in the table with some example values
    int numPoints = spectra2.size();
    table->SetNumberOfRows(numPoints);
    real_t en = 0.0;
    for (int i = 0; i < numPoints; ++i)
    {
        table->SetValue(i, 0, en);
        table->SetValue(i, 1, spectra1[i]);
        table->SetValue(i, 2, spectra2[i]);

        //table->SetValue(i, 1, spectra1->buffer()->operator [](i));
        //table->SetValue(i, 2, spectra2->buffer()->operator [](i));
        en += 1.0;
    }

//    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//    renderer->SetBackground(1.0, 1.0, 1.0);
//    vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->SetOffScreenRendering( 1 );
//    renderWindow->AddRenderer(renderer);

    // Set up the view
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
//    view->SetRenderer(renderer);
//    view->SetRenderWindow(renderWindow);
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
    view->GetRenderWindow()->SetOffScreenRendering( 1 );

    // Add multiple line plots, setting the colors etc
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 2);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(1.0);

    // For dotted line, the line type can be from 2 to 5 for different dash/dot
    // patterns (see enum in vtkPen containing DASH_LINE, value 2):
    #ifndef WIN32
    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
    #endif
    // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
    //  machines with built-in Intel HD graphics card...)

    //view->GetRenderWindow()->SetMultiSamples(0);

    view->GetRenderWindow()->SetSize(1024, 768);

    // Start interactor
    //view->GetInteractor()->Initialize();
    //view->GetInteractor()->Start();
    view->Render();

    //renderer->Render();

    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(view->GetRenderWindow());
    windowToImageFilter->Update();

    vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(path.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

}

/*
void PlotSpectras_carray(std::valarray<real_t> spectra1, std::valarray<real_t> spectra2, data_struct::xrf::Fit_Counts_Array* counts_arr)
{

    // Create a table with some points in it
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
    arrC->SetName("Integrated Spectra");
    table->AddColumn(arrC);

    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
    arrS->SetName("Spectra Model");
    table->AddColumn(arrS);

    vtkSmartPointer<vtkFloatArray> arr01 = vtkSmartPointer<vtkFloatArray>::New();
    arr01->SetName("Background");
    table->AddColumn(arr01);
    vtkSmartPointer<vtkFloatArray> arr02 = vtkSmartPointer<vtkFloatArray>::New();
    arr02->SetName("Ka");
    table->AddColumn(arr02);
    vtkSmartPointer<vtkFloatArray> arr03 = vtkSmartPointer<vtkFloatArray>::New();
    arr03->SetName("Kb");
    table->AddColumn(arr03);
    vtkSmartPointer<vtkFloatArray> arr04 = vtkSmartPointer<vtkFloatArray>::New();
    arr04->SetName("L");
    table->AddColumn(arr04);
    vtkSmartPointer<vtkFloatArray> arr05 = vtkSmartPointer<vtkFloatArray>::New();
    arr05->SetName("M");
    table->AddColumn(arr05);
    vtkSmartPointer<vtkFloatArray> arr06 = vtkSmartPointer<vtkFloatArray>::New();
    arr06->SetName("Elastic");
    table->AddColumn(arr06);
    vtkSmartPointer<vtkFloatArray> arr07 = vtkSmartPointer<vtkFloatArray>::New();
    arr07->SetName("Compton");
    table->AddColumn(arr07);

    // Fill in the table with some example values
    int numPoints = spectra2.size();
    table->SetNumberOfRows(numPoints);
    real_t en = 0.0;
    for (int i = 0; i < numPoints; ++i)
    {
        table->SetValue(i, 0, en);
        table->SetValue(i, 1, spectra1[i]);
        table->SetValue(i, 2, spectra2[i]);

        table->SetValue(i, 3, counts_arr->background[i]);
        table->SetValue(i, 4, counts_arr->ka[i]);
        table->SetValue(i, 5, counts_arr->kb[i]);
        table->SetValue(i, 6, counts_arr->l[i]);
        table->SetValue(i, 7, counts_arr->m[i]);
        table->SetValue(i, 8, counts_arr->elastic[i]);
        table->SetValue(i, 9, counts_arr->compton[i]);





   //     step
   //     tail
   //     pileup
   //     escape


        en += 1.0;
    }

    // Set up the view
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Add multiple line plots, setting the colors etc
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 0, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 2);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 3);
    line->SetColor(0, 0, 255, 255);
    line->SetWidth(1.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 4);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(1.0);

    // For dotted line, the line type can be from 2 to 5 for different dash/dot
    // patterns (see enum in vtkPen containing DASH_LINE, value 2):
    #ifndef WIN32
    line->GetPen()->SetLineType(vtkPen::DASH_LINE);
    #endif
    // (ifdef-ed out on Windows because DASH_LINE does not work on Windows
    //  machines with built-in Intel HD graphics card...)

    //view->GetRenderWindow()->SetMultiSamples(0);

    // Start interactor
    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();
}
*/

} //namespace visual
