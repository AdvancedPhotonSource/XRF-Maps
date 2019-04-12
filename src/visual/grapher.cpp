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

#include "grapher.h"
#include <algorithm>

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QCategoryAxis>
#include <QtCharts/QLogValueAxis>
#include <QtCharts/QValueAxis>
#include <QApplication>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QBarCategoryAxis>

namespace visual
{

void SavePlotSpectras(std::string path, data_struct::ArrayXr *energy, data_struct::ArrayXr *spectra, data_struct::ArrayXr *model, data_struct::ArrayXr *background, bool log_them)
{
    int argc = 0;
    char ** argv = nullptr;
    QApplication app(argc, argv);

    QtCharts::QLogValueAxis *axisYLog10 = new QtCharts::QLogValueAxis();
    axisYLog10->setTitleText("Counts Log10");
    axisYLog10->setLabelFormat("%.1e");
    //axisYLog10->setRange(1.0, 10000.0);
    axisYLog10->setBase(10.0);

    QtCharts::QValueAxis *axisX = new QtCharts::QValueAxis();
    axisX->setTitleText("Energy (keV)");
    axisX->setLabelFormat("%1.0f");
    //axisX->setTickCount(series->count());
    //axisX->setRange(0, 2048);
    axisX->setTickCount(11);

    QtCharts::QValueAxis *axisY = new QtCharts::QValueAxis();
    axisY->setTitleText("Counts");
    axisY->setLabelFormat("%i");
    //axisY->setTickCount(series->count());

    QtCharts::QChartView *chartView = new QtCharts::QChartView();
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->resize(1024, 768);

    QtCharts::QChart *chart = chartView->chart();
    chart->addAxis(axisX, Qt::AlignBottom);

    if(log_them)
    {
        chart->addAxis(axisYLog10, Qt::AlignLeft);
    }
    else
    {
        chart->addAxis(axisY, Qt::AlignLeft);
    }

    QtCharts::QLineSeries *series_spectra = new QtCharts::QLineSeries();
    QtCharts::QLineSeries *series_model = new QtCharts::QLineSeries();
    QtCharts::QLineSeries *series_background = new QtCharts::QLineSeries();

    series_spectra->setName("Integrated Spectra");
    series_spectra->setColor(QColor(Qt::green));
    series_model->setName("Model Spectra");
    series_model->setColor(QColor(Qt::red));
    series_background->setName("Background");
    series_background->setColor(QColor(Qt::blue));

    for(unsigned int i =0; i < spectra->rows(); i++)
    {

        float val_spec = (*spectra)[i];
        float val_mod = (*model)[i];
        float val_back = (*background)[i];

        bool isFine = std::isfinite(val_spec);
        if (false == isFine || val_spec <= 0.0f)
        {
            val_spec = 1.0;
        }
        isFine = std::isfinite(val_mod);
        if (false == isFine || val_mod <= 0.0f)
        {
            val_mod = 1.0;
        }
        isFine = std::isfinite(val_back);
        if (false == isFine || val_back <= 0.0f)
        {
            val_back = 1.0;
        }

        if(energy !=nullptr)
        {
            series_spectra->append((*energy)[i], val_spec);
            series_model->append((*energy)[i], val_mod);
            series_background->append((*energy)[i], val_back);
        }
        else
        {
            series_spectra->append(i, val_spec);
            series_model->append(i, val_mod);
            series_background->append(i, val_back);
        }
    }

    chart->addSeries(series_spectra);
    chart->addSeries(series_model);
    chart->addSeries(series_background);

    series_spectra->attachAxis(axisX);
    series_model->attachAxis(axisX);
    series_background->attachAxis(axisX);


    if(log_them)
    {
        series_spectra->attachAxis(axisYLog10);
        series_model->attachAxis(axisYLog10);
        series_background->attachAxis(axisYLog10);
    }
    else
    {
        series_spectra->attachAxis(axisY);
        series_model->attachAxis(axisY);
        series_background->attachAxis(axisY);
    }

    QPixmap pix = chartView->grab();
    QPainter painter(&pix);
    int h = 768;
    int w = 1024;
    int x = 0;
    int y = 0;
    painter.drawPixmap(x, y, w, h, pix);

    painter.end();
    pix.save(QString(path.c_str()), "png");

}

// ----------------------------------------------------------------------------

bool contains_shell(quantification::models::Electron_Shell shell_idx, string proc_type, int quant_id, map<string, data_struct::Quantification_Standard *> *standards)
{
    for(auto& s_itr: *standards)
    {
        data_struct::Quantification_Standard* quant_standard = s_itr.second;
        QtCharts::QScatterSeries *e_series = new QtCharts::QScatterSeries();
        e_series->setName(QString::fromStdString(s_itr.first));
        for(auto& itr : quant_standard->fitted_e_cal_ratio.at(proc_type).at(quant_id))
        {
            quantification::models::Electron_Shell shell = quantification::models::get_shell_by_name(itr.first);
            if(shell == shell_idx)
            {
                return true;
            }
        }
    }
    return false;
}

// ----------------------------------------------------------------------------

void SavePlotQuantification(std::string path, map<string, data_struct::Quantification_Standard *> *standards, int detector_num)
{
    const auto&itr = standards->begin();
    data_struct::Quantification_Standard* standard = itr->second;

    //iterate through proc_type {roi, nnls, fitted}
    for(auto& itr1 : standard->quantifier_map)
    {
        //iterate through quantifier {sr_current, us_ic, ds_ic}
        for(auto& itr2 : itr1.second.calib_curves)
        {
            if(contains_shell(quantification::models::K_SHELL, itr1.first, itr2.quant_id, standards))
            {
                std::string str_path_full = path + "calib_"+itr1.first+"_"+itr2.quantifier_name+"_K_det"+std::to_string(detector_num)+".png";
                SavePlotCalibrationCurve(str_path_full, standards, &itr2, itr1.first, quantification::models::K_SHELL, 13, 30);
            }
            if(contains_shell(quantification::models::L_SHELL, itr1.first, itr2.quant_id, standards))
            {
                std::string str_path_full = path + "calib_"+itr1.first+"_"+itr2.quantifier_name+"_L_det"+std::to_string(detector_num)+".png";
                SavePlotCalibrationCurve(str_path_full, standards, &itr2, itr1.first, quantification::models::L_SHELL, 39, 58);
            }
        }
    }
}

// ----------------------------------------------------------------------------

void SavePlotCalibrationCurve(std::string path,
                              map<string, data_struct::Quantification_Standard *> *standards,
                              data_struct::Calibration_Curve *calib_curve,
                              string proc_type,
                              quantification::models::Electron_Shell shell_idx,
                              int zstart, int zstop)
{
    int width_res = 1920;
    int height_res = 1080;
    //index starts at 0 so subtract one so we have correct z number
    zstart--;
    zstop--;
    int argc = 0;
    char ** argv = nullptr;
    QApplication app(argc, argv);

    QtCharts::QLogValueAxis *axisYLog10 = new QtCharts::QLogValueAxis();
    axisYLog10->setTitleText("Log10");
    axisYLog10->setLabelFormat("%.1e");
    axisYLog10->setBase(10.0);

    //QtCharts::QValueAxis *axisY = new QtCharts::QValueAxis();


    QtCharts::QChartView *chartView = new QtCharts::QChartView();
    chartView->setRenderHint(QPainter::Antialiasing);
    chartView->resize(width_res, height_res);

    QtCharts::QChart *chart = chartView->chart();
    chart->addAxis(axisYLog10, Qt::AlignLeft);
    //chart->addAxis(axisY, Qt::AlignLeft);


    QtCharts::QBarCategoryAxis *axisX = new QtCharts::QBarCategoryAxis();
    QStringList categories;
    QtCharts::QBarSet *set0 = new QtCharts::QBarSet("Element");

    QtCharts::QBarSeries *series = new QtCharts::QBarSeries();
    series->setName(QString::fromStdString(calib_curve->quantifier_name));

    real_t min_y = 9999999.0;
    real_t max_y = -9999999.0;
    for(int i=zstart; i <= zstop; i++)
    {
        categories << QString::fromStdString(calib_curve->shell_curves_labels[shell_idx][i]);
        *set0 << calib_curve->shell_curves[shell_idx][i];
        min_y = std::min(min_y, calib_curve->shell_curves[shell_idx][i]);
        max_y = std::max(max_y, calib_curve->shell_curves[shell_idx][i]);
    }
    series->append(set0);
    chart->addSeries(series);

    axisX->append(categories);
    chart->addAxis(axisX, Qt::AlignBottom);
    series->attachAxis(axisX);
    series->attachAxis(axisYLog10);
    //series->attachAxis(axisY);


    for(auto& s_itr: *standards)
    {
        data_struct::Quantification_Standard* quant_standard = s_itr.second;
        QtCharts::QScatterSeries *e_series = new QtCharts::QScatterSeries();
        e_series->setName(QString::fromStdString(s_itr.first));
        for(auto& itr : quant_standard->fitted_e_cal_ratio.at(proc_type).at(calib_curve->quant_id))
        {
            data_struct::Element_Info* element_info = data_struct::Element_Info_Map::inst()->get_element(itr.first);
            quantification::models::Electron_Shell shell = quantification::models::get_shell_by_name(itr.first);
            if(element_info != nullptr && shell == shell_idx)
            {
                e_series->append(((element_info->number -1) - zstart), itr.second);
                min_y = std::min(min_y, itr.second);
                max_y = std::max(max_y, itr.second);
            }
        }
        chart->addSeries(e_series);
        e_series->attachAxis(axisX);
        e_series->attachAxis(axisYLog10);
        //e_series->attachAxis(axisY);
    }
    //min_y -= 1.0;
    //max_y += 0.01;

    axisYLog10->setMin(min_y);
    axisYLog10->setMax(max_y);


    QPixmap pix = chartView->grab();
    QPainter painter(&pix);
    int h = height_res;
    int w = width_res;
    int x = 0;
    int y = 0;
    painter.drawPixmap(x, y, w, h, pix);

    painter.end();
    pix.save(QString(path.c_str()), "png");

}


} //namespace visual
