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

const real_t AVOGADRO = (real_t)6.02204531e23;
const real_t HC_ANGSTROMS = (real_t)12398.52;
const real_t RE = (real_t)2.817938070e-13;		// in cm

struct DLL_EXPORT Element_Info
{
	Element_Info();

	~Element_Info();

	void init_f_energies(int len);
	void init_extra_energies(int len);
    void get_energies_between(real_t energy, real_t* out_low, real_t* out_high, size_t* out_low_idx, size_t* out_high_idx);

    real_t calc_beta(real_t density_val, real_t energy);

    int number;
    std::string name;
    real_t density;
    real_t mass;
    std::unordered_map<std::string, real_t> xrf;
    std::unordered_map<std::string, real_t> xrf_abs_yield;
    std::unordered_map<std::string, real_t> yieldD;
    std::unordered_map<std::string, real_t> bindingE;
    std::unordered_map<std::string, real_t> jump;

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

//singleton
class DLL_EXPORT Element_Info_Map
{
public:

    static Element_Info_Map* inst();

	~Element_Info_Map();

	void clear();

	void generate_default_elements(int start_element, int end_element);

	void add_element(Element_Info* element);

    real_t calc_beta(std::string element_name, real_t density, real_t energy);

    Element_Info* get_element(int element_number);

    Element_Info* get_element(std::string element_name);

    //void set_energies(float* energy_arr, int num_energies);

    bool contains(std::string element_name) {return _name_element_info_map.count(element_name) > 0 ? true : false; }

    std::vector<float> _energies;

private:

    Element_Info_Map();

    static Element_Info_Map *_this_inst;

    std::unordered_map<std::string, Element_Info*> _name_element_info_map;
    std::map<int, Element_Info*>  _number_element_info_map;


};
                              // placeholder so H starts at 1
const std::string Element_Symbols[] = {" ", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
						                        "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
				                                "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
												"Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
												"Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
												"Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
												"Cf", "Es", "Fm", "Md", "No", "Lr", "Unq", "Unp", "Unh", "Uns", "Uno", "Une"};

const std::unordered_map<std::string, real_t> Henke_Compound_Density_Map = {
//Compounds
    std::pair<std::string, real_t>{"water",       (real_t)1.0},
    std::pair<std::string, real_t>{"protein",     (real_t)1.35},
    std::pair<std::string, real_t>{"lipid",       (real_t)1.0},
    std::pair<std::string, real_t>{"nucleosome",  (real_t)1.5},
    std::pair<std::string, real_t>{"dna",         (real_t)1.7},
    std::pair<std::string, real_t>{"helium",      (real_t)1.66E-04},
    std::pair<std::string, real_t>{"chromatin",   (real_t)1.527},
    std::pair<std::string, real_t>{"air",         (real_t)1.20E-03},
    std::pair<std::string, real_t>{"pmma",        (real_t)1.18},
    std::pair<std::string, real_t>{"nitride",     (real_t)3.44},
    std::pair<std::string, real_t>{"graphite",    (real_t)2.26},
    std::pair<std::string, real_t>{"nickel",      (real_t)8.876},
    std::pair<std::string, real_t>{"beryl",       (real_t)1.845},
    std::pair<std::string, real_t>{"copper",      (real_t)8.96},
    std::pair<std::string, real_t>{"quartz",      (real_t)2.2},
    std::pair<std::string, real_t>{"aluminum",    (real_t)2.7},
    std::pair<std::string, real_t>{"gold",        (real_t)19.3},
    std::pair<std::string, real_t>{"ice",         (real_t)0.92},
    std::pair<std::string, real_t>{"carbon",      (real_t)1.0},
    std::pair<std::string, real_t>{"polystyrene", (real_t)1.06},
    std::pair<std::string, real_t>{"silicon",     (real_t)2.33},
    std::pair<std::string, real_t>{"germanium",   (real_t)5.323},
//Formulas
    std::pair<std::string, real_t>{"H2O",                            (real_t)1.0},
    std::pair<std::string, real_t>{"H48.6C32.9N8.9O8.9S0.6",         (real_t)1.35},
    std::pair<std::string, real_t>{"H62.5C31.5O6.3",                 (real_t)1.0},
    std::pair<std::string, real_t>{"H42.1C31.9N10.3O13.9P1.6S0.3",   (real_t)1.5},
    std::pair<std::string, real_t>{"H35.5C30.8N11.7O18.9P3.1", (real_t)1.7},
    std::pair<std::string, real_t>{"He", (real_t)1.66E-04},
    std::pair<std::string, real_t>{"H49.95C24.64N8.66O15.57P1.07S0.03", (real_t)1.527},
    std::pair<std::string, real_t>{"N78.08O20.95Ar0.93", (real_t)1.20E-03},
    std::pair<std::string, real_t>{"C5H8O2", (real_t)1.18},
    std::pair<std::string, real_t>{"Si3N4", (real_t)3.44},
    std::pair<std::string, real_t>{"Ni", (real_t)8.876},
    std::pair<std::string, real_t>{"Be", (real_t)1.845},
    std::pair<std::string, real_t>{"Cu", (real_t)8.96},
    std::pair<std::string, real_t>{"SiO2", (real_t)2.2},
    std::pair<std::string, real_t>{"Al", (real_t)2.7},
    std::pair<std::string, real_t>{"Au", (real_t)19.3},
    std::pair<std::string, real_t>{"C", (real_t)1.0},
    std::pair<std::string, real_t>{"C8H8", (real_t)1.06},
    std::pair<std::string, real_t>{"Si", (real_t)2.33},
    std::pair<std::string, real_t>{"Ge", (real_t)5.323}
};


const std::unordered_map<int, real_t> Element_Weight = {
    std::pair<int, real_t>(1,  (real_t)1.00794),
    std::pair<int, real_t>(2,  (real_t)4.0026),
    std::pair<int, real_t>(3,  (real_t)6.941),
    std::pair<int, real_t>(4,  (real_t)9.01218),
    std::pair<int, real_t>(5,  (real_t)10.81),
    std::pair<int, real_t>(6,  (real_t)12.011),
    std::pair<int, real_t>(7,  (real_t)14.0067),
    std::pair<int, real_t>(8,  (real_t)15.9994),
    std::pair<int, real_t>(9,  (real_t)18.9984),
    std::pair<int, real_t>(10, (real_t)21.179),
    std::pair<int, real_t>(11, (real_t)22.98977),
    std::pair<int, real_t>(12, (real_t)24.305),
    std::pair<int, real_t>(13, (real_t)26.98154),
    std::pair<int, real_t>(14, (real_t)28.0855),
    std::pair<int, real_t>(15, (real_t)30.97376),
    std::pair<int, real_t>(16, (real_t)32.06),
    std::pair<int, real_t>(17, (real_t)35.453),
    std::pair<int, real_t>(18, (real_t)39.948),
    std::pair<int, real_t>(19, (real_t)39.0983),
    std::pair<int, real_t>(20, (real_t)40.08),
    std::pair<int, real_t>(21, (real_t)44.9559),
    std::pair<int, real_t>(22, (real_t)47.88),
    std::pair<int, real_t>(23, (real_t)50.9415),
    std::pair<int, real_t>(24, (real_t)51.996),
    std::pair<int, real_t>(25, (real_t)54.9380),
    std::pair<int, real_t>(26, (real_t)55.847),
    std::pair<int, real_t>(27, (real_t)58.9332),
    std::pair<int, real_t>(28, (real_t)58.69),
    std::pair<int, real_t>(29, (real_t)63.546),
    std::pair<int, real_t>(30, (real_t)65.38),
    std::pair<int, real_t>(31, (real_t)69.72),
    std::pair<int, real_t>(32, (real_t)72.59),
    std::pair<int, real_t>(33, (real_t)74.9216),
    std::pair<int, real_t>(34, (real_t)78.96),
    std::pair<int, real_t>(35, (real_t)79.904),
    std::pair<int, real_t>(36, (real_t)83.80),
    std::pair<int, real_t>(37, (real_t)85.4678),
    std::pair<int, real_t>(38, (real_t)87.62),
    std::pair<int, real_t>(39, (real_t)88.9059),
    std::pair<int, real_t>(40, (real_t)91.22),
    std::pair<int, real_t>(41, (real_t)92.9064),
    std::pair<int, real_t>(42, (real_t)95.94),
    std::pair<int, real_t>(43, (real_t)98.0),
    std::pair<int, real_t>(44, (real_t)101.07),
    std::pair<int, real_t>(45, (real_t)102.9055),
    std::pair<int, real_t>(46, (real_t)106.42),
    std::pair<int, real_t>(47, (real_t)107.8682),
    std::pair<int, real_t>(48, (real_t)112.41),
    std::pair<int, real_t>(49, (real_t)114.82),
    std::pair<int, real_t>(50, (real_t)118.69),
    std::pair<int, real_t>(51, (real_t)121.75),
    std::pair<int, real_t>(52, (real_t)127.60),
    std::pair<int, real_t>(53, (real_t)126.9054),
    std::pair<int, real_t>(54, (real_t)131.29),
    std::pair<int, real_t>(55, (real_t)132.9054),
    std::pair<int, real_t>(56, (real_t)137.33),
    std::pair<int, real_t>(57, (real_t)138.9055),
    std::pair<int, real_t>(58, (real_t)140.12),
    std::pair<int, real_t>(59, (real_t)140.9077),
    std::pair<int, real_t>(60, (real_t)144.24),
    std::pair<int, real_t>(61, (real_t)145.0),
    std::pair<int, real_t>(62, (real_t)150.36),
    std::pair<int, real_t>(63, (real_t)151.96),
    std::pair<int, real_t>(64, (real_t)157.25),
    std::pair<int, real_t>(65, (real_t)158.9254),
    std::pair<int, real_t>(66, (real_t)162.5),
    std::pair<int, real_t>(67, (real_t)164.9304),
    std::pair<int, real_t>(68, (real_t)167.26),
    std::pair<int, real_t>(69, (real_t)168.9342),
    std::pair<int, real_t>(70, (real_t)173.04),
    std::pair<int, real_t>(71, (real_t)174.967),
    std::pair<int, real_t>(72, (real_t)178.49),
    std::pair<int, real_t>(73, (real_t)180.9479),
    std::pair<int, real_t>(74, (real_t)183.85),
    std::pair<int, real_t>(75, (real_t)186.207),
    std::pair<int, real_t>(76, (real_t)190.2),
    std::pair<int, real_t>(77, (real_t)192.22),
    std::pair<int, real_t>(78, (real_t)195.08),
    std::pair<int, real_t>(79, (real_t)196.9665),
    std::pair<int, real_t>(80, (real_t)200.59),
    std::pair<int, real_t>(81, (real_t)204.383),
    std::pair<int, real_t>(82, (real_t)207.2),
    std::pair<int, real_t>(83, (real_t)208.9804),
    std::pair<int, real_t>(84, (real_t)209.0),
    std::pair<int, real_t>(85, (real_t)210.0),
    std::pair<int, real_t>(86, (real_t)222.0),
    std::pair<int, real_t>(87, (real_t)223.0),
    std::pair<int, real_t>(88, (real_t)226.0254),
    std::pair<int, real_t>(89, (real_t)227.0278),
    std::pair<int, real_t>(90, (real_t)232.0381),
    std::pair<int, real_t>(91, (real_t)231.0359),
    std::pair<int, real_t>(92, (real_t)238.0289),
    std::pair<int, real_t>(93, (real_t)237.0),
    std::pair<int, real_t>(94, (real_t)244.0),
    std::pair<int, real_t>(95, (real_t)243.0),
    std::pair<int, real_t>(96, (real_t)247.0),
    std::pair<int, real_t>(97, (real_t)247.0),
    std::pair<int, real_t>(98, (real_t)251.0),
    std::pair<int, real_t>(99, (real_t)252.0),
    std::pair<int, real_t>(100, (real_t)257.0),
    std::pair<int, real_t>(101, (real_t)258.0),
    std::pair<int, real_t>(102, (real_t)259.0),
    std::pair<int, real_t>(103, (real_t)262.0),
    std::pair<int, real_t>(104, (real_t)261.0),
    std::pair<int, real_t>(105, (real_t)262.0),
    std::pair<int, real_t>(106, (real_t)266.0),
    std::pair<int, real_t>(107, (real_t)264.0),
    std::pair<int, real_t>(108, (real_t)277.0),
    std::pair<int, real_t>(109, (real_t)268.0)
};

} //namespace data_struct

#endif // Element_Info_H
