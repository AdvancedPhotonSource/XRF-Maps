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



#ifndef Element_Info_H
#define Element_Info_H

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include "core/defines.h"

namespace data_struct
{

template<typename T_real>
struct DLL_EXPORT Element_Info
{
	Element_Info();

	~Element_Info();

	void init_f_energies(int len);
	void init_extra_energies(int len);
    void get_energies_between(T_real energy, T_real* out_low, T_real* out_high, size_t* out_low_idx, size_t* out_high_idx);

    T_real calc_beta(T_real density_val, T_real energy);

    T_real get_f2(T_real energy);

    int number;
    std::string name;
    T_real density;
    T_real mass;
    std::unordered_map<std::string, T_real> xrf;
    std::unordered_map<std::string, T_real> xrf_abs_yield;
    std::unordered_map<std::string, T_real> yieldD;
    std::unordered_map<std::string, T_real> bindingE;
    std::unordered_map<std::string, T_real> jump;

    //int atomic_scattering_len;
    //atomic scattering factors real
    std::vector<float> f1_atomic_scattering_real; //has to be float because that is how it is saved in henke.xdr
    //atomic scattering factors imaginary
    std::vector<float> f2_atomic_scattering_imaginary; //has to be float because that is how it is saved in henke.xdr

    //int extra_energies_len;
    std::vector<float> *energies;
    std::vector<float> extra_energies;
    std::vector<float> extra_f1;
    std::vector<float> extra_f2;
};

//TEMPLATE_STRUCT_DLL_EXPORT Element_Info<float>;
//TEMPLATE_STRUCT_DLL_EXPORT Element_Info<double>;


//singleton
template<typename T_real>
class DLL_EXPORT Element_Info_Map
{
public:

    static Element_Info_Map* inst();

	~Element_Info_Map();

	void clear();

	void generate_default_elements(int start_element, int end_element);

	void add_element(Element_Info<T_real>* element);

    T_real calc_beta(std::string element_name, T_real density, T_real energy);

    T_real calc_compound_beta(std::string compound_name, T_real density, T_real energy);

    Element_Info<T_real>* get_element(int element_number);

    Element_Info<T_real>* get_element(std::string element_name);

    bool is_element(std::string element_name);

    //void set_energies(float* energy_arr, int num_energies);

    bool contains(std::string element_name) {return _name_element_info_map.count(element_name) > 0 ? true : false; }

    std::vector<float> _energies;

private:

    Element_Info_Map();

    static Element_Info_Map *_this_inst;

    std::unordered_map<std::string, Element_Info<T_real>*> _name_element_info_map;
    std::map<int, Element_Info<T_real>*>  _number_element_info_map;


};
   
TEMPLATE_CLASS_DLL_EXPORT Element_Info_Map<float>;
TEMPLATE_CLASS_DLL_EXPORT Element_Info_Map<double>;


// placeholder so H starts at 1
const std::string Element_Symbols[] = {" ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
						                        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
				                                "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
												"Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
												"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
												"Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
												"Cf", "Es", "Fm", "Md", "No", "Lr", "Unq", "Unp", "Unh", "Uns", "Uno", "Une"};

template<typename T_real>
const std::unordered_map<std::string, T_real> Henke_Compound_Density_Map = {
//Compounds
    std::pair<std::string, T_real>{"water",       (T_real)1.0},
    std::pair<std::string, T_real>{"protein",     (T_real)1.35},
    std::pair<std::string, T_real>{"lipid",       (T_real)1.0},
    std::pair<std::string, T_real>{"nucleosome",  (T_real)1.5},
    std::pair<std::string, T_real>{"dna",         (T_real)1.7},
    std::pair<std::string, T_real>{"helium",      (T_real)1.66E-04},
    std::pair<std::string, T_real>{"chromatin",   (T_real)1.527},
    std::pair<std::string, T_real>{"air",         (T_real)1.20E-03},
    std::pair<std::string, T_real>{"pmma",        (T_real)1.18},
    std::pair<std::string, T_real>{"nitride",     (T_real)3.44},
    std::pair<std::string, T_real>{"graphite",    (T_real)2.26},
    std::pair<std::string, T_real>{"nickel",      (T_real)8.876},
    std::pair<std::string, T_real>{"beryl",       (T_real)1.845},
    std::pair<std::string, T_real>{"copper",      (T_real)8.96},
    std::pair<std::string, T_real>{"quartz",      (T_real)2.2},
    std::pair<std::string, T_real>{"aluminum",    (T_real)2.7},
    std::pair<std::string, T_real>{"gold",        (T_real)19.3},
    std::pair<std::string, T_real>{"ice",         (T_real)0.92},
    std::pair<std::string, T_real>{"carbon",      (T_real)1.0},
    std::pair<std::string, T_real>{"polystyrene", (T_real)1.06},
    std::pair<std::string, T_real>{"silicon",     (T_real)2.33},
    std::pair<std::string, T_real>{"germanium",   (T_real)5.323},
//Formulas
    std::pair<std::string, T_real>{"H:2,O",                            (T_real)1.0},
    std::pair<std::string, T_real>{"H:48.6,C:32.9,N:8.9,O:8.9,S:0.6",         (T_real)1.35},
    std::pair<std::string, T_real>{"H:62.5C31.5O6.3",                 (T_real)1.0},
    std::pair<std::string, T_real>{"H:42.1C31.9N10.3O13.9P1.6S0.3",   (T_real)1.5},
    std::pair<std::string, T_real>{"H:35.5C30.8N11.7O18.9P3.1", (T_real)1.7},
    std::pair<std::string, T_real>{"He", (T_real)1.66E-04},
    std::pair<std::string, T_real>{"H:49.95,C:24.64,N:8.66,O:15.57,P:1.07,S:0.03", (T_real)1.527},
    std::pair<std::string, T_real>{"N:78.08,O:20.95,Ar:0.93", (T_real)1.20E-03},
    std::pair<std::string, T_real>{"C:5,H:8,O:2", (T_real)1.18},
    std::pair<std::string, T_real>{"Si:3,N:4", (T_real)3.44},
    std::pair<std::string, T_real>{"Ni", (T_real)8.876},
    std::pair<std::string, T_real>{"Be", (T_real)1.845},
    std::pair<std::string, T_real>{"Cu", (T_real)8.96},
    std::pair<std::string, T_real>{"SiO2", (T_real)2.2},
    std::pair<std::string, T_real>{"Al", (T_real)2.7},
    std::pair<std::string, T_real>{"Au", (T_real)19.3},
    std::pair<std::string, T_real>{"C", (T_real)1.0},
    std::pair<std::string, T_real>{"C:8,H:8", (T_real)1.06},
    std::pair<std::string, T_real>{"Si", (T_real)2.33},
    std::pair<std::string, T_real>{"Ge", (T_real)5.323}
};

template<typename T_real>
const std::unordered_map<int, T_real> Element_Weight = {
    std::pair<int, T_real>(1,  (T_real)1.00794),
    std::pair<int, T_real>(2,  (T_real)4.0026),
    std::pair<int, T_real>(3,  (T_real)6.941),
    std::pair<int, T_real>(4,  (T_real)9.01218),
    std::pair<int, T_real>(5,  (T_real)10.81),
    std::pair<int, T_real>(6,  (T_real)12.011),
    std::pair<int, T_real>(7,  (T_real)14.0067),
    std::pair<int, T_real>(8,  (T_real)15.9994),
    std::pair<int, T_real>(9,  (T_real)18.9984),
    std::pair<int, T_real>(10, (T_real)21.179),
    std::pair<int, T_real>(11, (T_real)22.98977),
    std::pair<int, T_real>(12, (T_real)24.305),
    std::pair<int, T_real>(13, (T_real)26.98154),
    std::pair<int, T_real>(14, (T_real)28.0855),
    std::pair<int, T_real>(15, (T_real)30.97376),
    std::pair<int, T_real>(16, (T_real)32.06),
    std::pair<int, T_real>(17, (T_real)35.453),
    std::pair<int, T_real>(18, (T_real)39.948),
    std::pair<int, T_real>(19, (T_real)39.0983),
    std::pair<int, T_real>(20, (T_real)40.08),
    std::pair<int, T_real>(21, (T_real)44.9559),
    std::pair<int, T_real>(22, (T_real)47.88),
    std::pair<int, T_real>(23, (T_real)50.9415),
    std::pair<int, T_real>(24, (T_real)51.996),
    std::pair<int, T_real>(25, (T_real)54.9380),
    std::pair<int, T_real>(26, (T_real)55.847),
    std::pair<int, T_real>(27, (T_real)58.9332),
    std::pair<int, T_real>(28, (T_real)58.69),
    std::pair<int, T_real>(29, (T_real)63.546),
    std::pair<int, T_real>(30, (T_real)65.38),
    std::pair<int, T_real>(31, (T_real)69.72),
    std::pair<int, T_real>(32, (T_real)72.59),
    std::pair<int, T_real>(33, (T_real)74.9216),
    std::pair<int, T_real>(34, (T_real)78.96),
    std::pair<int, T_real>(35, (T_real)79.904),
    std::pair<int, T_real>(36, (T_real)83.80),
    std::pair<int, T_real>(37, (T_real)85.4678),
    std::pair<int, T_real>(38, (T_real)87.62),
    std::pair<int, T_real>(39, (T_real)88.9059),
    std::pair<int, T_real>(40, (T_real)91.22),
    std::pair<int, T_real>(41, (T_real)92.9064),
    std::pair<int, T_real>(42, (T_real)95.94),
    std::pair<int, T_real>(43, (T_real)98.0),
    std::pair<int, T_real>(44, (T_real)101.07),
    std::pair<int, T_real>(45, (T_real)102.9055),
    std::pair<int, T_real>(46, (T_real)106.42),
    std::pair<int, T_real>(47, (T_real)107.8682),
    std::pair<int, T_real>(48, (T_real)112.41),
    std::pair<int, T_real>(49, (T_real)114.82),
    std::pair<int, T_real>(50, (T_real)118.69),
    std::pair<int, T_real>(51, (T_real)121.75),
    std::pair<int, T_real>(52, (T_real)127.60),
    std::pair<int, T_real>(53, (T_real)126.9054),
    std::pair<int, T_real>(54, (T_real)131.29),
    std::pair<int, T_real>(55, (T_real)132.9054),
    std::pair<int, T_real>(56, (T_real)137.33),
    std::pair<int, T_real>(57, (T_real)138.9055),
    std::pair<int, T_real>(58, (T_real)140.12),
    std::pair<int, T_real>(59, (T_real)140.9077),
    std::pair<int, T_real>(60, (T_real)144.24),
    std::pair<int, T_real>(61, (T_real)145.0),
    std::pair<int, T_real>(62, (T_real)150.36),
    std::pair<int, T_real>(63, (T_real)151.96),
    std::pair<int, T_real>(64, (T_real)157.25),
    std::pair<int, T_real>(65, (T_real)158.9254),
    std::pair<int, T_real>(66, (T_real)162.5),
    std::pair<int, T_real>(67, (T_real)164.9304),
    std::pair<int, T_real>(68, (T_real)167.26),
    std::pair<int, T_real>(69, (T_real)168.9342),
    std::pair<int, T_real>(70, (T_real)173.04),
    std::pair<int, T_real>(71, (T_real)174.967),
    std::pair<int, T_real>(72, (T_real)178.49),
    std::pair<int, T_real>(73, (T_real)180.9479),
    std::pair<int, T_real>(74, (T_real)183.85),
    std::pair<int, T_real>(75, (T_real)186.207),
    std::pair<int, T_real>(76, (T_real)190.2),
    std::pair<int, T_real>(77, (T_real)192.22),
    std::pair<int, T_real>(78, (T_real)195.08),
    std::pair<int, T_real>(79, (T_real)196.9665),
    std::pair<int, T_real>(80, (T_real)200.59),
    std::pair<int, T_real>(81, (T_real)204.383),
    std::pair<int, T_real>(82, (T_real)207.2),
    std::pair<int, T_real>(83, (T_real)208.9804),
    std::pair<int, T_real>(84, (T_real)209.0),
    std::pair<int, T_real>(85, (T_real)210.0),
    std::pair<int, T_real>(86, (T_real)222.0),
    std::pair<int, T_real>(87, (T_real)223.0),
    std::pair<int, T_real>(88, (T_real)226.0254),
    std::pair<int, T_real>(89, (T_real)227.0278),
    std::pair<int, T_real>(90, (T_real)232.0381),
    std::pair<int, T_real>(91, (T_real)231.0359),
    std::pair<int, T_real>(92, (T_real)238.0289),
    std::pair<int, T_real>(93, (T_real)237.0),
    std::pair<int, T_real>(94, (T_real)244.0),
    std::pair<int, T_real>(95, (T_real)243.0),
    std::pair<int, T_real>(96, (T_real)247.0),
    std::pair<int, T_real>(97, (T_real)247.0),
    std::pair<int, T_real>(98, (T_real)251.0),
    std::pair<int, T_real>(99, (T_real)252.0),
    std::pair<int, T_real>(100, (T_real)257.0),
    std::pair<int, T_real>(101, (T_real)258.0),
    std::pair<int, T_real>(102, (T_real)259.0),
    std::pair<int, T_real>(103, (T_real)262.0),
    std::pair<int, T_real>(104, (T_real)261.0),
    std::pair<int, T_real>(105, (T_real)262.0),
    std::pair<int, T_real>(106, (T_real)266.0),
    std::pair<int, T_real>(107, (T_real)264.0),
    std::pair<int, T_real>(108, (T_real)277.0),
    std::pair<int, T_real>(109, (T_real)268.0)
};

} //namespace data_struct

#endif // Element_Info_H
