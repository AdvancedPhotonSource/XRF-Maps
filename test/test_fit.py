from numpy.lib.shape_base import column_stack
import pyxrfmaps as px
import h5py
import matplotlib.pyplot as plt
import numpy as np

element_csv_filename = "../reference/xrf_library.csv"
element_henke_filename = "../reference/henke.xdr"
detector_num = 0
dataset_dir = 'D:/bnp/MinghettiTi_fly49_XRFMaps/'
dataset = 'bnp_fly0049.mda.h50'
full_path = dataset_dir + 'img.dat/' + dataset

def load_dataset(filename):
    # Load integreated spectra from dataset
    df = h5py.File(filename)
    if not df == None:
        if '/MAPS/Spectra/Integrated_Spectra/Spectra' in df:
            return df['/MAPS/Spectra/Integrated_Spectra/Spectra'][...]
        elif '/MAPS/int_spec' in df:
            return df['/MAPS/int_spec'][...]
    else:
        print('Failed to open')
    return None

def plot_results(int_spec, fit_spec):
    i_ax = np.linspace(0,int_spec.size-1, int_spec.size)
    f_ax = np.linspace(0,fit_spec.size-1, fit_spec.size)
    fig, axs = plt.subplots(2,3)
    axs[0,0].plot(i_ax, int_spec)
    axs[0,1].plot(i_ax, int_spec)
    axs[0,1].set_yscale('log')
    
    fft_spec = np.fft.fft(int_spec)
    freq = np.fft.fftfreq(fft_spec.size)
    axs[0,2].plot(freq, fft_spec.real**2 + fft_spec.imag**2)
    
    axs[0,0].plot(f_ax, fit_spec)
    axs[0,1].plot(f_ax, fit_spec)
    #axs[0,1].set_yscale('log')
    
    fft_fit_spec = np.fft.fft(fit_spec)
    axs[0,2].plot(freq, fft_fit_spec.real**2 + fft_fit_spec.imag**2)
    
    diff_spec = np.abs(int_spec - fit_spec)
    axs[1,0].plot(i_ax, diff_spec)
    axs[1,1].plot(i_ax, diff_spec)
    axs[1,1].set_yscale('log')
    #print(fft_fit_spec.imag)
    ffdiff = fft_spec - fft_fit_spec
    axs[1,2].plot(freq, ffdiff.real**2+ ffdiff.imag**2)
    

    plt.show()

if __name__ == '__main__':
    # initialize element info
    px.load_element_info(element_henke_filename, element_csv_filename)

    # Select fitting routine
    #fit_rout = px.fitting.routines.nnls()
    #fit_rout = px.fitting.routines.svd()
  ###  opt = px.fitting.optimizers.mpfit()
    fit_rout = px.fitting.routines.matrix()
    fit_rout.set_optimizer(opt)
    #fit_rout = px.fitting.routines.roi()
    
    # Use Gausian Model
    model = px.fitting.models.GaussModel()
    
    # Load fit parameters 
    po = px.load_override_params(dataset_dir, detector_num, True)
    
    # Load dataset
    int_spec = load_dataset(full_path)
    
    # Initialize model and fit routine with fit parameters
    energy_range = px.get_energy_range(int_spec.size, po.fit_params)
    model.update_fit_params_values(po.fit_params)
    fit_rout.initialize(model, po.elements_to_fit, energy_range)
    
    # Fit element counts
    counts = fit_rout.fit_counts(model, int_spec, po.elements_to_fit)
    print(counts)
    
    # Get Fit Spectra 
    fit_spec = fit_rout.fit_spectra(model, int_spec, po.elements_to_fit)

    # Resize int_spec to match fit_spec
    int_spec = int_spec[energy_range.min:energy_range.max+1]

    # Plot both
    plot_results(int_spec, fit_spec)