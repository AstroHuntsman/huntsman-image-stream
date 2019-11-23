from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from gunagala.imager import create_imagers
from gunagala.config import load_config
import gunagala.psf

from photutils.datasets import make_gaussian_sources_image

import numpy as np


def get_star_table(width,
                   height,
                   field_coordinates,
                   catalog=['I/340/ucac5'],
                   mag='Gmag',
                   mag_limit="<10",
                   row_limit=-1):

    v = Vizier(columns=['_RAJ2000', '_DEJ2000', mag],
               column_filters={mag: mag_limit},
               row_limit=row_limit)

    result = v.query_region(field_coordinates,
                            width=width,
                            height=height,
                            catalog=catalog)

    return(result[0])


def create_psf(field_coordinates):
    table = Table()
    table['amplitude'] = [10]
    table['x_mean'] = [50]
    table['y_mean'] = [50]
    table['x_stddev'] = [2.]
    table['y_stddev'] = [2.]
    table['theta'] = np.array([0.]) * np.pi / 180.

    psf_data = make_gaussian_sources_image((100, 100), table)

    psf = gunagala.psf.PixellatedPSF(psf_data,
                                     psf_sampling=2 * u.arcsec / u.pixel,
                                     oversampling=10)
    return(psf)


def make_noiseless_image(imager, star_table):

    imager
    pass


def main(exptime=0.005 * u.second):

    field_target_name = 'fornax cluster'
    coordinate_table = Simbad.query_object(field_target_name)
    field_coordinates = SkyCoord(coordinate_table['RA'][0],
                                 coordinate_table['DEC'][0],
                                 unit=(u.deg, u.deg))

    imager_filter = 'g'
    imagers = create_imagers(load_config(
        '/Users/lspitler/prog/GitHub/huntsman-ms/resources/performance_ms.yaml'))
    imager = imagers['one_zwo_canon_full_moon']

    # use digital PSF
    psf = create_psf(field_coordinates)
    imager = gunagala.imager.Imager(optic=imager.optic,
                                    camera=imager.camera,
                                    filters=imager.filters,
                                    psf=psf,
                                    sky=imager.sky)
    imager.set_WCS_centre(field_coordinates)

    gband_limit = imager.point_source_limit(total_exp_time=exptime,
                                            sub_exp_time=exptime,
                                            filter_name=imager_filter,
                                            snr_target=1.0)
    gband_limit_string = f'<{gband_limit.value}'
    ucac_filter_name = f'{imager_filter.upper()}mag'

    star_table = get_star_table(width=imager.field_of_view[0],
                                height=imager.field_of_view[1],
                                field_coordinates=field_coordinates,
                                mag=ucac_filter_name,
                                mag_limit=gband_limit_string)

    print(star_table)
    tuple_star_list = []
    for star in star_table:
        star_coord = SkyCoord(star['_RAJ2000'] * u.deg, star['_RAJ2000'] * u.deg)
        tuple_star_list.append((star_coord, star[ucac_filter_name]))

    # make_noiseless_image(imager, star_table)
    image_data = imager.make_noiseless_image(centre=field_coordinates,
                                             obs_time=exptime,
                                             filter_name=imager_filter,
                                             stars=tuple_star_list)
    image_data.write('out.fits')


if __name__ == '__main__':
    main()
