from astroquery.vizier import Vizier
from astroquery.simbad import Simbad

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.nddata import CCDData


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


def make_noiseless_image(imager,
                         imager_filter_name,
                         star_table):
    ucac_filter_name = f'{ imager_filter_name.upper() }mag'

    electrons = np.zeros((imager.wcs._naxis2,
                          imager.wcs._naxis1)) * u.electron / (u.second * u.pixel)

    # Calculate observed sky background
    sky_rate = imager.sky_rate[imager_filter_name]
    if hasattr(imager.sky, 'relative_brightness'):
        pixel_coords = imager.get_pixel_coords()
        relative_sky = imager.sky.relative_brightness(pixel_coords, obs_time)
        sky_rate = sky_rate * relative_sky
    electrons = electrons + sky_rate

    # compute stellar fluxes and locations
    star_rate = imager.ABmag_to_rate(
        star_table[ucac_filter_name].data * u.ABmag, imager_filter_name)
    pixel_coords = imager.wcs.all_world2pix(star_table['_RAJ2000'].data * u.deg,
                                            star_table['_DEJ2000'].data * u.deg, 0)

    table = Table()
    table['amplitude'] = star_rate.value
    table['x_mean'] = pixel_coords[0]
    table['y_mean'] = pixel_coords[1]
    table['x_stddev'] = np.ones(len(pixel_coords[0]))  # PSF width = 1 pixels
    table['y_stddev'] = np.ones(len(pixel_coords[0]))  # PSF width = 1 pixels
    table['theta'] = np.zeros(len(pixel_coords[0]))

    star_data = make_gaussian_sources_image((imager.wcs._naxis2,
                                             imager.wcs._naxis1), table) * star_rate.unit / u.pixel

    electrons = electrons + star_data
    noiseless = CCDData(electrons, wcs=imager.wcs)

    return(noiseless)


def create_noiseless_image(exptime=0.005 * u.second,
                           snr_limit=1.,
                           gunagala_imager_name='one_zwo_canon_full_moon',
                           gunagala_config_filename='/Users/lspitler/prog/GitHub/huntsman-ms/resources/performance_ms.yaml',
                           imager_filter_name='g',
                           field_target_name='fornax cluster',
                           output_fits_filename='out_noiseless.fits',
                           output_fits_file=False,
                           write_region_file=False):

    coordinate_table = Simbad.query_object(field_target_name)
    field_coordinates = SkyCoord(coordinate_table['RA'][0],
                                 coordinate_table['DEC'][0],
                                 unit=(u.deg, u.deg))

    imagers = create_imagers(load_config(gunagala_config_filename))
    imager = imagers[gunagala_imager_name]
    imager.set_WCS_centre(field_coordinates)

    band_limit = imager.point_source_limit(total_exp_time=exptime,
                                           sub_exp_time=exptime,
                                           filter_name=imager_filter_name,
                                           snr_target=snr_limit)
    band_limit_string = f'<{band_limit.value}'
    ucac_filter_name = f'{imager_filter_name.upper()}mag'

    # 5% buffer on edge
    star_table = get_star_table(width=imager.field_of_view[0] * .95,
                                height=imager.field_of_view[1] * .95,
                                field_coordinates=field_coordinates,
                                mag=ucac_filter_name,
                                mag_limit=band_limit_string)

    if write_region_file:
        star_table.write('region.reg',
                         format='ascii.fast_no_header',
                         include_names=['_RAJ2000', '_DEJ2000'])

    image_data = make_noiseless_image(imager,
                                      imager_filter_name,
                                      star_table)

    if output_fits_file:
        image_data.write(output_fits_filename, overwrite=True)

    real_data = imager.make_image_real(image_data, exptime)
    real_data.write('out_real.fits', overwrite=True)


if __name__ == '__main__':
    create_noiseless_image()
