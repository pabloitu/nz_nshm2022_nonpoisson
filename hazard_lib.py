import numpy as np
from pathlib import Path
import shutil
import subprocess
import os
import h5py

from openquake.hazardlib.source.point import PointSource
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.nodalplane import NodalPlane
from openquake.hazardlib.pmf import PMF
from openquake.hazardlib.scalerel.point import PointMSR
from openquake.hazardlib.mfd.evenly_discretized import EvenlyDiscretizedMFD
from openquake.hazardlib import sourcewriter
from openquake.hazardlib.tom import PoissonTOM, NegativeBinomialTOM


def run_model(folder_name):
    """
    Runs a folder where the openquake files are found: Job.ini, source model, source logic tree, grid and gmpe logic tree
    :param folder_name: Location of Openquake source files
    :return:
    """
    job_path = os.path.join(folder_name, 'job.ini')
    print('Running Openquake')

    oq_proc = subprocess.call(['oq', 'engine', '--run', job_path], stdout=subprocess.PIPE)
    if oq_proc == 0:
        # If success, gets calculation id
        process_2 = subprocess.run(['oq', 'engine', '--lhc'], stdout=subprocess.PIPE)
        out = str(process_2).split(r'\n')[-2]
        line = out.split()
        print(line)
        calc_id = int(line[0])
        src = os.path.join(Path.home(), 'oqdata', f'calc_{calc_id}.hdf5')
        # Copy calculation hdf5 database into model folder
        shutil.copy(src, folder_name)
        print(f'Done. hdf5 database copied to {folder_name}')
    else:
        #if oq throws error
        raise Exception('Run failed')


def parse_db(folder):
    """
    Reads a hdf5 database and parses the mean curves.  Intensity measure levels (imtl) are fixed.
    :param folder:
    :return:
    """
    files = list(os.walk(folder))[0][2]
    calcs = [i for i in files if 'calc' in i]
    ids = [int(i.split('_')[1].split('.')[0]) for i in calcs]
    calc_path = os.path.join(folder, calcs[int(np.argmax(ids))])


    db = h5py.File(os.path.join(folder, calc_path), 'r')

    grid = np.stack((db['sitecol']['lon'], db['sitecol']['lat'])).T

    hcurves = db['hcurves-stats'][:, 0, 0, :]

    imtl = np.logspace(-2, 0.2, 30)   # Same as all models' job.ini

    data = {'grid': grid,
            'hcurves': hcurves,
            'imtl': imtl}
    return data


class hazard_model:

    def __init__(self, name=None):
        self.name = name
        self.cells = np.array(0, np.dtype(float))
        self.rates = np.array(0, np.dtype(float))
        self.dispersion = np.array(0, np.dtype(float))
        self.rates_mbin = np.array(0, np.dtype(float))
        self.sources = []

        self.mmin = 5.
        self.mmax = 8.
        self.mbin = 0.1
        self.magnitudes = np.arange( self.mmin,  self.mmax + self.mbin,  self.mbin)

    @property
    def nsources(self):
        return self.cells.shape[0]

    def read_forecast(self, filename, eepas=False):
        """
        Parses forecast files
        :param forecast_name: path to forecast
        :param eepas: Flag to indicate if the forecast is just a spatial PDF (False) or, magnitude bins spatial PDFs (True)

        """
        if eepas is False:
            data = np.genfromtxt(filename, delimiter=',', skip_header=1)
            self.cells = np.array([[i[0] + i[1], i[2] + i[3]] for i in data[:, :4]])/2.
            self.rates = data[:, 8]
            self.dispersion = data[:, 9]

        else:
            data = np.genfromtxt(filename, delimiter=',', skip_header=1)
            self.cells = np.array([[i[0] + i[1], i[2] + i[3]] for i in data[:, :4]])/2.
            self.rates_mbin = data[:, 6:-1]
            self.rates = np.sum(self.rates_mbin, axis=1)

            self.dispersion = data[:, -1]

    def scale(self, scale):
        """
        Scales rates and magnitude rates
        :param scale:
        :return:
        """
        self.rates *= scale
        self.rates_mbin *= scale

    def set_trunc_gr(self, bval, eepas=False):
        """
        Projects spatial rates using a truncated Gutenberg-Richter magnitude-frequency-distribution into magnitude bins.
        If the model contains EEPAS, scales the magnitude bins pdf scaling each magnitude PDF by the weights of a truncated GR distribution
        :param bval: b-value of the distribution
        :return:
        """

        mags = self.magnitudes
        mbin = self.mbin
        mmin = self.mmin
        mmax = self.mmax
        rates = self.rates

        weights = (10 ** (- bval * (mags - mbin / 2.)) - 10 ** (- bval * (mags + mbin / 2.))) / (
                   10 ** (- bval * (mmin - mbin / 2.)) - 10 ** (- bval * (mmax + mbin / 2.)))

        if eepas is False:
            rates_mbin = np.zeros((rates.shape[0], mags.shape[0]))
            for i, mu in enumerate(rates):
                rates_mbin[i, :] = mu * weights
            self.rates_mbin = rates_mbin

        else:
            rates_mbin = self.rates_mbin

            for i, mu_m in enumerate(rates_mbin):
                self.rates_mbin[i, :] = weights * mu_m
                self.rates[i] = np.sum(self.rates_mbin[i, :])

    def get_point_srcs(self,
                       time_span=1,
                       tect_region='Active Shallow Crust',
                       rupt_mesh_spacing=5.0,
                       magnitude_scale_rel=PointMSR(),
                       rupt_aspect_ratio=1.0,
                       upper_seis_depth=0,
                       lower_seis_depth=30,
                       nodal_dist=[(1, NodalPlane(0, 90, 0.))],
                       hypo_dist=((0.5, 10.0), (0.5, 30.0))):
        """
        Creates the Openquake objects necessary to write a source nrml file. Assign a Negative Binomial temporal model
        if the dispersion is greater than 0, and a Poisson temporal model if 0.
        """
        npd = PMF(nodal_dist)
        hd = PMF(hypo_dist)

        for n in range(self.nsources):
            mfd = EvenlyDiscretizedMFD(self.mmin, self.mbin, self.rates_mbin[n, :])
            if self.dispersion[n] != 0:
                mu_tot = self.rates[n]
                alpha = self.dispersion[n]
                tom = NegativeBinomialTOM(time_span, mu_tot, alpha)
            else:
                tom = PoissonTOM(1)

            source = PointSource(source_id='%05i' % n,
                                 name='point%05i' % n,
                                 tectonic_region_type=tect_region,
                                 mfd=mfd,
                                 rupture_mesh_spacing=rupt_mesh_spacing,
                                 magnitude_scaling_relationship=magnitude_scale_rel,
                                 rupture_aspect_ratio=rupt_aspect_ratio,
                                 temporal_occurrence_model=tom,
                                 upper_seismogenic_depth=upper_seis_depth,
                                 lower_seismogenic_depth=lower_seis_depth,
                                 location=Point(self.cells[n, 0],
                                                self.cells[n, 1]),
                                 nodal_plane_distribution=npd,
                                 hypocenter_distribution=hd)
            self.sources.append(source)

    def write_source(self, source_path):

        sourcewriter.write_source_model(source_path,
                                        self.sources,
                                        name=self.name,
                                        investigation_time=1)


