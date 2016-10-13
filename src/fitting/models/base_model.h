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

#ifndef Base_Model_H
#define Base_Model_H

#include "spectra_volume.h"
#include "fit_parameters.h"
#include "fit_element_map.h"
#include "calibration_standard.h"

namespace fitting
{
namespace models
{

using namespace data_struct::xrf;
using namespace std;

/**
 * @brief The Range struct to determine size of spectra we want to fit or model
 */
struct Range
{
    size_t count() const  {return (max - min) + 1;}
    int min;
    int max;
};

/**
 * @brief get_energy_range: genereates a range which consists of min and max. This represents the min energy and max enegry of the spectra to fit.
 * @param min_energy
 * @param max_energy
 * @param spectra_size
 * @param calibration: energy calibration
 * @return Range structure with the min energy and max enegry of the spectra to fit.
 */
DLL_EXPORT Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, const Calibration_Standard * const  calibration);

/**
 * @brief The Base_Model class: base class for modeling spectra and fitting elements
 */
class DLL_EXPORT Base_Model
{
public:
    /**
     * @brief Base_Model : Constructor
     */
    Base_Model();

    /**
     * @brief ~Base_Model : Destructor
     */
    ~Base_Model();

    /**
     * @brief fit_spectra : Fit a single specra ( typically 2048 in size )
     * @param fit_params : Fitting parameters required by the routine
     * @param spectra : Pointer to the spectra we are fitting to
     * @param calibration : Energy calibration
     * @param elements_to_fit : List of elemetns to fit to the spectra. This is an out variable also. Must be allocated to saved fitted value to using row_idx and col_idx
     * @param row_idx : row index used to save the fitted value back into elements_to_fit class
     * @param col_idx : column index used to save the fitted value back into elements_to_fit class
     */
    Fit_Parameters fit_spectra(const Fit_Parameters fit_params,
                              const Spectra * const spectra,
                              const Calibration_Standard * const calibration,
                              const Fit_Element_Map_Dict * const elements_to_fit,
                              Fit_Count_Dict *out_counts_dic,
                              size_t row_idx=0,
                              size_t col_idx=0);

    /**
     * @brief fit_spectra_volume : fit a row/col volume of spectra
     * @param fit_params : Fitting parameters required by the routine: Also the last spectras fitted values
     * @param spectra_volume : volume of spectra to be fitted
     * @param calibration : Energy calibration
     * @param elements_to_fit : List of elemetns to fit to the spectra. This is an out variable also. Must be allocated to saved fitted value to using row_idx and col_idx
     * @param col_idx_start : Start index of the column we want to fit
     * @param col_idx_end : End index of the column we want to fit
     * @param row_idx_start : Start index of the rows we want to fit
     * @param row_idx_end : End index of the rows we want to fit
     */
    void fit_spectra_volume(const Fit_Parameters fit_params,
                            const Spectra_Volume * const spectra_volume,
                            const Calibration_Standard * const calibration,
                            const Fit_Element_Map_Dict * const elements_to_fit,
                            Fit_Count_Dict *out_counts_dic,
                            size_t row_idx,
                            size_t col_idx);

    /**
     * @brief get_fit_parameters : returns Fit_Parameters class of the required fit parameters to run a fitting
     * @return
     */
    virtual Fit_Parameters get_fit_parameters() = 0;

    /**
     * @brief model_spectrum : Model a spectra based on the fit parameters passed in.
     * @param fit_params : Fitting parameters required to model the spectra.
     * @param spectra : Might be depricated when I remove the snip background function
     * @param calibration : Energy calibration
     * @param elements_to_fit : List of elemetns to use in modeling the spectra.
     * @param energy_range : Spectra model energy range. Basically the size of the spectra model returned;
     * @return
     */
    virtual Spectra model_spectrum(const Fit_Parameters * const fit_params,
                                   const Spectra * const spectra,
                                   const Calibration_Standard * const calibration,
                                   const Fit_Element_Map_Dict * const elements_to_fit,
                                   const struct Range energy_range) = 0;

    /**
     * @brief initialize : Initialize the model
     * @param fit_params
     * @param calibration
     * @param elements_to_fit
     * @param energy_range
     */
    virtual void initialize(Fit_Parameters *fit_params,
                            const Calibration_Standard * const calibration,
                            const Fit_Element_Map_Dict * const elements_to_fit,
                            const struct Range energy_range);

protected:

    /**
     * @brief _pre_process : Preprocessing function derived classes can use before fitting spectra. All arguments can be updated before calling fitting.
     * @param fit_params : Fitting parameters that will be passed to the fitting routine.
     * @param calibration : Energy calibration
     * @param elements_to_fit : Elements to fit
     */
    virtual void _pre_process(Fit_Parameters *fit_params,
                              const Spectra * const spectra,
                              const Calibration_Standard * const calibration,
                              const Fit_Element_Map_Dict * const elements_to_fit);

    /**
     * @brief _fit_spectra : Fitting routine created by derived class
     * @param fit_params : Fitting parameters that will be passed to the fitting routine.
     * @param spectra : Spectra to compare spectra model with
     * @param calibration : Energy Calibration
     * @param elements_to_fit : Elements to fit in spectra
     */
    virtual void _fit_spectra(Fit_Parameters *fit_params,
                              const Spectra * const spectra,
                              const Calibration_Standard * const calibration,
                              const Fit_Element_Map_Dict * const elements_to_fit) = 0;

    /**
     * @brief _post_process : Post processing function called after _fit_spectra
     * @param fit_params
     * @param spectra
     * @param calibration
     * @param elements_to_fit
     */
    virtual void _post_process(Fit_Parameters *fit_params,
                               const Spectra * const spectra,
                               const Calibration_Standard * const calibration,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               Fit_Count_Dict *out_counts_dic,
                               size_t row_idx,
                               size_t col_idx);

    /**
     * @brief _add_elements_to_fit_parameters : Adds each element to fit as a fit parameter if they don't already exist in fit_params
     * @param fit_params : elements we will be fitting
     * @param spectra
     * @param calibration : Energy calibration
     * @param elements_to_fit : elements to add to fit parameters
     */
    void _add_elements_to_fit_parameters(Fit_Parameters *fit_params,
                                         const Spectra * const spectra,
                                         const Calibration_Standard * const calibration,
                                         const Fit_Element_Map_Dict * const elements_to_fit);

    /**
     * @brief _update_elements_guess : Updates the guess values for each element to fit based on spectra details
     * @param fit_params : Stores the guess values for each element
     * @param spectra : used to create guess values
     * @param calibration : Energy calibration
     * @param elements_to_fit : elements to update
     */
    void _update_elements_guess(Fit_Parameters *fit_params,
                                const Spectra * const spectra,
                                const Calibration_Standard * const calibration,
                                const Fit_Element_Map_Dict * const elements_to_fit);

    /**
     * @brief _update_element_guess_value : Protected variable used to determine if function _update_elements_guess is called for each spectra in
     * fit_spectra_line or fit_spectra_volume
     */
    bool _update_element_guess_value;

    /**
     * @brief _counts_log_10 : tell if the model expects counts in log10
     */
    bool _counts_log_10;

    /**
     * @brief _save_counts_per_sec: should we divide counts by elapsed life time
     */
    bool _save_counts_per_sec;

private:


};

} //namespace models

} //namespace fitting

#endif // Base_Model_H
