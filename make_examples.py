from os.path import join, dirname, realpath
from hazard_lib import hazard_model, run_model, parse_db

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"xtick.bottom": True, "ytick.left": True,
                     "xtick.color": 'darkgrey', 'xtick.labelcolor': 'black',
                     "ytick.color": 'darkgrey', 'ytick.labelcolor': 'black'})


def multiplicative_model():
    """
    Parses the Multiplicative Hybrid forecast, writes the openquake file and
     runs the hazard calculation.
    :return:
    """

    N = 5.1
    bval = 0.929

    folder = dirname(realpath(__file__))

    model_folder = join(folder, 'examples/multiplicative')
    forecast_path = join(folder, 'forecasts/hybrids/src/m.csv')
    oqsource_path = join(model_folder, 'source.xml')
    grid_path = join(model_folder, 'grid.txt')

    m = hazard_model('Multiplicative')
    m.read_forecast(forecast_path)
    m.set_trunc_gr(bval)
    m.scale(N)
    m.get_point_srcs()
    m.write_source(oqsource_path)
    Auckland = np.array([[174.7, -36.8]])
    np.savetxt(grid_path, Auckland, fmt='%.2f')
    run_model(model_folder)


def poisson_floor():
    """
    Parses the Poisson Floor Ensemble forecast  for the multiplicative model,
     writes the openquake file and runs the hazard calculation.
    :return:
    """

    N = 5.1
    bval = 0.929
    folder = dirname(realpath(__file__))

    model_folder = join(folder, 'examples/poisson_floor_multiplicative')
    forecast_path = join(folder, 'forecasts/poisson_floors/src/fe_m.csv')
    oqsource_path = join(model_folder, 'source.xml')
    grid_path = join(model_folder, 'grid.txt')

    m = hazard_model('Poisson Floor')
    m.read_forecast(forecast_path)
    m.set_trunc_gr(bval)
    m.scale(N)
    m.get_point_srcs()
    m.write_source(oqsource_path)

    Auckland = np.array([[174.7, -36.8]])
    np.savetxt(grid_path, Auckland, fmt='%.2f')
    run_model(model_folder)


def negbinom_floor():
    """
    Parses the Negative Binomial Floor Ensemble forecast for the multiplicative
     model, writes the openquake file and runs the hazard calculation.
    :return:
    """

    N = 5.1
    bval = 0.929
    folder = dirname(realpath(__file__))

    model_folder = join(folder, 'examples/negbinom_floor_multiplicative')
    forecast_path = join(folder, 'forecasts/negbinom_floors/src/npfe_m.csv')
    oqsource_path = join(model_folder, 'source.xml')
    grid_path = join(model_folder, 'grid.txt')

    m = hazard_model('Negbinom Floor')
    m.read_forecast(forecast_path)
    m.set_trunc_gr(bval)
    m.scale(N)
    m.get_point_srcs()
    m.write_source(oqsource_path)

    Auckland = np.array([[174.7, -36.8]])
    np.savetxt(grid_path, Auckland, fmt='%.2f')
    run_model(model_folder)


def negbinom_floor_eepas():
    """
    Parses the Negative Binomial Floor Ensemble forecast for the multiplicative
     model, writes the openquake file and runs the hazard calculation.
    :return:
    """

    N = 5.1
    bval = 0.929
    folder = dirname(realpath(__file__))

    model_folder = join(folder, 'examples/negbinom_floor_eepas_multiplicative')
    forecast_path = join(folder,
                         'forecasts/negbinom_floors/src/npfe_m_eepas.csv')
    oqsource_path = join(model_folder, 'source.xml')
    grid_path = join(model_folder, 'grid.txt')

    m = hazard_model('Negbinom Floor')
    m.read_forecast(forecast_path, eepas=True)
    m.set_trunc_gr(bval, eepas=True)
    m.scale(N)
    m.get_point_srcs()
    m.write_source(oqsource_path)

    Auckland = np.array([[174.7, -36.8]])
    np.savetxt(grid_path, Auckland, fmt='%.2f')
    run_model(model_folder)


def plot_results():
    folder = dirname(realpath(__file__))
    multiplicative_folder = join(folder, 'examples/multiplicative')
    mult_data = parse_db(multiplicative_folder)

    poisson_floor_folder = join(folder,
                                'examples/poisson_floor_multiplicative')
    fe_m_data = parse_db(poisson_floor_folder)

    negbinom_floor_folder = join(folder,
                                 'examples/negbinom_floor_multiplicative')
    npfe_m_data = parse_db(negbinom_floor_folder)

    negbinom_floor_folder = join(
                 folder,
                 'examples/negbinom_floor_eepas_multiplicative')
    npfe_m_eepas_data = parse_db(negbinom_floor_folder)

    fig, ax = plt.subplots()
    ax.loglog(mult_data['imtl'], mult_data['hcurves'][0], color='blue',
              label='Multiplicative')
    ax.loglog(mult_data['imtl'], fe_m_data['hcurves'][0], color='red',
              label='Poisson Floor')
    ax.loglog(mult_data['imtl'], npfe_m_data['hcurves'][0], color='green',
              label='Negbinom Floor')
    ax.loglog(mult_data['imtl'], npfe_m_eepas_data['hcurves'][0],
              color='purple', label='Negbinom EEPAS Floor')

    xlims = [1e-2, 1.2]
    ylims = [1e-5, 0.1]
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    ax.axhline(0.002105, linestyle='--', linewidth=0.8, color='black')
    ax.text(min(xlims) * 1.1, 0.002205, '10% in ' + '%i yr.' % 50, fontsize=12)
    ax.axhline(0.000404, linestyle='--', linewidth=0.8, color='black')
    ax.text(min(xlims) * 1.1, 0.000454, '2% in ' + '%i yr.' % 50, fontsize=12)
    ax.set_xlabel(f'PGA $[g]$', fontsize=14)
    ax.set_ylabel('Probability of exceedance - %i years' % 50, fontsize=14)
    ax.legend(loc='upper right', fontsize=14)
    ax.grid(visible=True, which='minor', color='white', linestyle='-',
            linewidth=0.5)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    print('Running Multiplicative model')
    multiplicative_model()
    print('Running Poisson Floor model')
    poisson_floor()
    print('Running Negbinom Floor model')
    negbinom_floor()
    print('Running Negbinom Floor model with EEPAS')
    negbinom_floor_eepas()
    print('Visualizing results')
    plot_results()
    print('Ready')
