"""
Inject a NR waveform into GW strain,
assume zero noise,
save to .gwf file and plot the waveform.
"""

import numpy as np
import h5py
# lal imports
import lal
import lalsimulation as lalsim
# pycbc imports
from pycbc.frame.frame import write_frame
from pycbc.waveform import get_td_waveform
from pycbc.waveform.utils import taper_timeseries
from pycbc.types import TimeSeries
from pycbc.types import float64, float32, zeros
injection_func_map = {
    np.dtype(float32): lambda *args: lalsim.SimAddInjectionREAL4TimeSeries(*args),
    np.dtype(float64): lambda *args: lalsim.SimAddInjectionREAL8TimeSeries(*args),
}
# bilby imports
import bilby
from gwpy.timeseries import TimeSeries as gwpy_timeseries


Msun_s = 4.925491025e-6  # 1 solar mass in seconds

def get_time_domain_data(ifo, waveform_arguments):
    strain = TimeSeries(zeros(ifo.duration * ifo.sampling_frequency), 
                        delta_t = 1./ifo.sampling_frequency,
                        epoch=ifo.start_time)
    
    lalstrain = strain.lal()
    add_injection = injection_func_map[strain.dtype]

    hplus, hcross = get_td_waveform(approximant = waveform_arguments['waveform_approximant'],
                                    numrel_data=waveform_arguments['numerical_relativity_file'],
                                    mass1 = waveform_arguments['mass_1'],
                                    mass2 = waveform_arguments['mass_2'],
                                    spin1x = waveform_arguments['spin_1x'],
                                    spin1y = waveform_arguments['spin_1y'],
                                    spin1z = waveform_arguments['spin_1z'],
                                    spin2x = waveform_arguments['spin_2x'],
                                    spin2y = waveform_arguments['spin_2y'],
                                    spin2z = waveform_arguments['spin_2z'],
                                    distance = waveform_arguments['luminosity_distance'],
                                    inclination = waveform_arguments['iota'],
                                    coa_phase = waveform_arguments['phase'],
                                    f_lower = waveform_arguments['minimum_frequency'],
                                    f_ref = waveform_arguments['reference_frequency'],
                                    delta_t = 1/ waveform_arguments['sampling_frequency'])
    i_mrg = np.argmax(hplus**2 + hcross**2)
    t_mrg = hplus.sample_times[i_mrg]
    # have merger (more or less) at geocent time
    hplus.start_time  += waveform_arguments['geocent_time'] - t_mrg
    hcross.start_time += waveform_arguments['geocent_time'] - t_mrg

    hplus_tapered  = taper_timeseries(hplus, 'TAPER_START',)
    hcross_tapered = taper_timeseries(hcross, 'TAPER_START',)
    
    waveform = {'plus': hplus_tapered, 'cross': hcross_tapered}
    
    time_delay = ifo.time_delay_from_geocenter(ra=waveform_arguments['ra'],
                                               dec=waveform_arguments['dec'],
                                               time=waveform_arguments['geocent_time'])
    signal = 0
    for key in waveform.keys():
        antenna_response = ifo.antenna_response(ra=waveform_arguments['ra'], 
                                                dec=waveform_arguments['dec'], 
                                                psi=waveform_arguments['psi'],
                                                time=waveform_arguments['geocent_time'],
                                                mode=key)
        signal += antenna_response * waveform[key]
    signal.start_time += time_delay    
    signal = signal.astype(strain.dtype)
    signal_lal = signal.lal()
    add_injection(lalstrain, signal_lal, None)

    strain.data[:] = lalstrain.data.data[:]
    return strain, waveform

def load_nr_metadata(filename, M=100):

    f = h5py.File(filename, 'r')

    metadata = {}

    Mt    = f.attrs['mass1'] + f.attrs['mass2']
    m1    = f.attrs['mass1']*M/Mt
    m2    = f.attrs['mass2']*M/Mt
    q     = m1/m2
    nu    = q/(1+q)**2
    flow  = f.attrs['f_lower_at_1MSUN']/M
    fref  = -1
    spins = lalsim.SimInspiralNRWaveformGetSpinsFromHDF5File(fref, M, filename)
    c1x, c1y, c1z = spins[0], spins[1], spins[2]
    c2x, c2y, c2z = spins[3], spins[4], spins[5]

    metadata['numerical_relativity_file'] = filename
    metadata['waveform_approximant'] = 'NR_hdf5'

    # overwrite some of the metadata
    metadata['mass_1']    = m1
    metadata['mass_2']    = m2
    metadata['M']     = m1+m2
    metadata['q']     = m1/m2
    metadata['nu']    = nu
    metadata['spin_1x'] = c1x; metadata['spin_1y'] = c1y;  metadata['spin_1z'] = c1z
    metadata['spin_2x'] = c2x; metadata['spin_2y'] = c2y;  metadata['spin_2z'] = c2z
    metadata['initial_frequency_geo'] = flow*Msun_s*M
    metadata['eccentricity']   = f.attrs['eccentricity']

    metadata['m1SI']   = m1 * lal.MSUN_SI
    metadata['m2SI']   = m2 * lal.MSUN_SI
    metadata['reference_frequency'] = flow

    return metadata

if __name__ == "__main__":

    M        = 80  # total mass in solar masses

    iota     = 0. # inclination angle in radians
    phase    = 0.0  # coalescence phase in radians
    distance = 500.0  # luminosity distance in Mpc
    sampling_frequency = 4096.0  # sampling frequency in Hz
    minimum_frequency  = 20.0  # minimum frequency in Hz
    duration = 4.0  # duration of the waveform in seconds
    gpstime  = 1000000000  # GPS time of the event
    dec      = 0.0
    ra       = 0.0
    psi      = 0.0  # polarization angle in radians
    # mode_array for the waveform
    mode_array = [[2, 2], [2, -2], 
                  [3, 3], [3, -3], 
                  [4, 4], [4, -4],
                ]
    
    gaussian_noise = False  # whether to add gaussian noise to the detector

    # open the file
    nr_fname = '../dat/Boson-star-waveforms/GRChombo/'
    #nr_fname += 'GRChombo_BBSsol02_A147A147q100d12p000_Res40.h5' 
    nr_fname +='GRChombo_BBSsol02_A17A17q100d17p180_Res40.h5'
    
    waveform_arguments = load_nr_metadata(nr_fname, M=M)
    waveform_arguments['iota'] = iota
    waveform_arguments['phase'] = phase
    waveform_arguments['luminosity_distance'] = distance
    waveform_arguments['sampling_frequency'] = sampling_frequency
    waveform_arguments['minimum_frequency'] = minimum_frequency

    waveform_arguments['duration'] = duration
    waveform_arguments['geocent_time'] = gpstime
    waveform_arguments['ra'] = ra
    waveform_arguments['dec'] = dec
    waveform_arguments['psi'] = psi

    print("Injection parameters:")
    for key, value in waveform_arguments.items():
        print(f"\t{key}: {value}")

    print("NR minimum frequency: ", waveform_arguments['reference_frequency'])
    print("Injection and data initial frequency: ", waveform_arguments['minimum_frequency'])

    # setup the detectors with bilby
    detectors = bilby.gw.detector.InterferometerList(['H1', 'L1', 'V1'])
    fname_ifo = {'H1': 'aLIGO_ZERO_DET_high_P_psd.txt',
                 'L1': 'aLIGO_ZERO_DET_high_P_psd.txt',
                 'V1': 'AdV_psd.txt'}

    network_optimal_snr_squared = 0
    for det in detectors:

        det.duration = waveform_arguments['duration']
        det.minimum_frequency = waveform_arguments['minimum_frequency']
        det.sampling_frequency = waveform_arguments['sampling_frequency']
        det.start_time = waveform_arguments['geocent_time'] - det.duration / 2
        psd=bilby.gw.detector.PowerSpectralDensity.from_power_spectral_density_file(fname_ifo[det.name])
        det.psd = psd
 
        if gaussian_noise:
            # add gaussian noise to the detector
            det.set_strain_data_from_power_spectral_density(
                duration=det.duration, 
                sampling_frequency=det.sampling_frequency, 
                start_time=det.start_time, 
            )
        else:
            # set the strain data to zero
            det.set_strain_data_from_zero_noise(
                                          duration=det.duration, 
                                          sampling_frequency=det.sampling_frequency, 
                                          start_time=det.start_time)
        ts, wf_dict = get_time_domain_data(det, waveform_arguments=waveform_arguments)

        ts += det.strain_data.time_domain_strain
        # save to gwpy
        write_frame("nr_injection_" + det.name + ".gwf", f'{det.name}:INJ', ts)

        # Compute the optimal SNR for the detector
        # idk if this is the right way with gaussian noise tbh
        ts = gwpy_timeseries.from_pycbc(ts)
        det.set_strain_data_from_gwpy_timeseries(ts)
        network_optimal_snr_squared += det.optimal_snr_squared(det.frequency_domain_strain)
        print(f"\tSNR for {det.name}: {np.sqrt(det.optimal_snr_squared(det.frequency_domain_strain))}")

    print(f"Network SNR: {np.sqrt(network_optimal_snr_squared)}")


    # load the gwf files and plot them, to make sure
    import matplotlib.pyplot as plt
    data = gwpy_timeseries.read('nr_injection_H1.gwf', 'H1:INJ')
    data.plot()
    plt.show()


