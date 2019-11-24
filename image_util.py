"""Image utilities.
"""
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
    """Query Vizier catalog for stars of a given brightness.

    Args:
        width (astropy.Quantity): Field width in u.deg units.
        height (astropy.Quantity): Field width in u.deg units.
        field_coordinates (SkyCoord): Location of field centre.
        catalog (list, optional): List of catalog names to query
        mag (str, optional): Filter to cut on.
        mag_limit (str, optional): Filter magnitude limit criteria.
        row_limit (int, optional): Number of stars to use. -1 is unlimited.

    Returns:
        astropy.Table: Table of RA, Dec, Magnitude
    """
    v = Vizier(columns=['_RAJ2000', '_DEJ2000', mag],
               column_filters={mag: mag_limit},
               row_limit=row_limit)

    result = v.query_region(field_coordinates,
                            width=width,
                            height=height,
                            catalog=catalog)

    return(result[0])


def make_noiseless_data(imager,
                        imager_filter_name,
                        star_table):
    """Generate noiseless data of star field. Based upon
    gunagala.make_noiseless_data, which currently doesn't do PSFs
    very well.

    Args:
        imager (gunagala.Imager): Instance of a gunagala imaging system.
        imager_filter_name (str): Filter to use for noiseless data.
        star_table (astropy.Table): Table of star coordinates and magnitudes.

    Returns:
        CCDData: Noiseless imaging data.
    """
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


def generate_noiseless_image(exptime=0.005 * u.second,
                             snr_limit=1.,
                             gunagala_imager_name='one_zwo_canon_full_moon',
                             gunagala_config_filename='/Users/lspitler/prog/GitHub/huntsman-ms/resources/performance_ms.yaml',
                             imager_filter_name='g',
                             field_target_name='fornax cluster',
                             output_fits_filename='out_noiseless.fits',
                             output_fits_file=False,
                             write_region_file=False):
    """Summary

    Args:
        exptime (astropy.Quantity, optional): Exposure time in u.seconds
        snr_limit (float, optional): Signal-to-noise lower limit for selecting stars.
        gunagala_imager_name (str, optional): name of gunagala Imager to use.
        gunagala_config_filename (str, optional): Filepath to gunagala yaml config.
        imager_filter_name (str, optional): Name of filter to use.
        field_target_name (str, optional): Name of astronomical target.
        output_fits_filename (str, optional): Output filename of noiseless fits image.
        output_fits_file (bool, optional): Write out noiseless fits image.
        write_region_file (bool, optional): Write out simple RA,Dec text file for DS9 region overlay
    """

    exptime = exptime.to(u.second)
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

    image_data = make_noiseless_data(imager,
                                     imager_filter_name,
                                     star_table)

    if output_fits_file:
        image_data.write(output_fits_filename, overwrite=True)

    return(image_data, imager)


if __name__ == '__main__':
    generate_noiseless_image()
