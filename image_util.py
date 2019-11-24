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


def create_psf(field_coordinates, pixel_scale, psf_sampling=0.5):
    table = Table()
    table['amplitude'] = [10]
    table['x_mean'] = [50]
    table['y_mean'] = [50]
    table['x_stddev'] = [2.]
    table['y_stddev'] = [2.]
    table['theta'] = np.array([0.]) * np.pi / 180.

    print(pixel_scale)
    psf_data = make_gaussian_sources_image((100, 100), table)

    psf = gunagala.psf.PixellatedPSF(psf_data,
                                     psf_sampling=0.5 * u.arcsec / u.pixel,
                                     pixel_scale=pixel_scale,
                                     oversampling=10)
    return(psf)


def make_noiseless_image(imager,
                         imager_filter_name,
                         star_table):
    ucac_filter_name = f'{ imager_filter_name.upper() }mag'

    electrons = np.zeros((imager.wcs._naxis2,
                          imager.wcs._naxis1)) * u.electron / (u.second * u.pixel) * 1.

    # Calculate observed sky background
    sky_rate = imager.sky_rate[imager_filter_name]
    if hasattr(imager.sky, 'relative_brightness'):
        pixel_coords = imager.get_pixel_coords()
        relative_sky = imager.sky.relative_brightness(pixel_coords, obs_time)
        sky_rate = sky_rate * relative_sky
    electrons = electrons + sky_rate

    star_rate = imager.ABmag_to_rate(
        star_table[ucac_filter_name].data * u.ABmag, imager_filter_name)
    pixel_coords = imager.wcs.all_world2pix(star_table['_RAJ2000'].data * u.deg,
                                            star_table['_DEJ2000'].data * u.deg, 0)
    # pixel_coords[0] = np.array(pixel_coords[0])  # - imager.wcs.wcs.crpix[0]
    # pixel_coords[1] = np.array(pixel_coords[1])  # - imager.wcs.wcs.crpix[1]

    table = Table()
    if 0:  # testing
        table['amplitude'] = [10]
        table['x_mean'] = [50]
        table['y_mean'] = [50]
        table['x_stddev'] = [2.]
        table['y_stddev'] = [2.]
        table['theta'] = np.array([0.]) * np.pi / 180.
    else:
        table['amplitude'] = star_rate.value
        table['x_mean'] = pixel_coords[0]
        table['y_mean'] = pixel_coords[1]
        table['x_stddev'] = np.ones(len(pixel_coords[0]))  # PSF width = 1
        table['y_stddev'] = np.ones(len(pixel_coords[0]))  # PSF width = 1
        table['theta'] = np.zeros(len(pixel_coords[0]))

    star_data = make_gaussian_sources_image((imager.wcs._naxis2,
                                             imager.wcs._naxis1), table) * star_rate.unit / u.pixel

    electrons = electrons + star_data

    noiseless = CCDData(electrons, wcs=imager.wcs)

    return(noiseless)


def make_image_real(imager, exp_time, noiseless_data):
    noisy_data = apply_poisson_noise(noiseless_data)
    noisy_data


def main(exptime=0.005 * u.second, snr_limit=1.):

    field_target_name = 'fornax cluster'
    coordinate_table = Simbad.query_object(field_target_name)
    field_coordinates = SkyCoord(coordinate_table['RA'][0],
                                 coordinate_table['DEC'][0],
                                 unit=(u.deg, u.deg))

    imager_filter_name = 'g'
    imagers = create_imagers(load_config(
        '/Users/lspitler/prog/GitHub/huntsman-ms/resources/performance_ms.yaml'))
    imager = imagers['one_zwo_canon_full_moon']

    # use digital PSF
    # psf = create_psf(field_coordinates, imager.pixel_scale)
    # psf = gunagala.psf.PixellatedPSF(imager.psf.pixellated(),
    #                                 psf_sampling=imager.pixel_scale,
    #                                 pixel_scale=imager.pixel_scale)

    # imager = gunagala.imager.Imager(optic=imager.optic,
    #                                camera=imager.camera,
    #                                filters=imager.filters,
    #                                psf=psf,
    #                                sky=imager.sky)
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

    star_table.write('region.reg', format='ascii', include_names=['_RAJ2000', '_DEJ2000'])
    # print(star_table['_RAJ2000'], star_table['_DEJ2000'])
    #tuple_star_list = []
    # for star in star_table:
    #    star_coord = SkyCoord(star['_RAJ2000'] * u.deg, star['_DEJ2000'] * u.deg)
    #    tuple_star_list.append((star_coord, star[ucac_filter_name] * u.ABmag))

    # print(tuple_star_list)

    image_data = make_noiseless_image(imager,
                                      imager_filter_name,
                                      star_table)
    # image_data = imager.make_noiseless_image(centre=field_coordinates,
    #                                         obs_time=exptime,
    #                                         filter_name=imager_filter,
    #                                         stars=tuple_star_list)

    real_data = imager.make_image_real(image_data, exptime)
    image_data.write('out_noiseless.fits', clobber=True)
    real_data.write('out_real.fits', clobber=True)


if __name__ == '__main__':
    main()
