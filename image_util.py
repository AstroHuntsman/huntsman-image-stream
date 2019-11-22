from astroquery.vizier import Vizier

import astropy.units as u

from gunagala.imager import create_imagers
from gunagala.config import load_config


def get_star_table(width,
                   height,
                   location=,
                   catalog=['I/340/ucac5'],
                   mag='Gmag',
                   mag_limit="<10",
                   row_limit=-1):

    v = Vizier(columns=['_RAJ2000', '_DEJ2000', mag],
               column_filters={mag: mag_limit},
               row_limit=row_limit)

    result = v.query_region(location,
                            width=width,
                            height=height,
                            catalog=catalog)

    return(result[0])


def make_noiseless_image(imager, star_table):

    imager
    pass


def main(exptime=0.005 * u.second):

    field_target_name = 'fornax cluster'

    imagers = create_imagers(load_config(
        '/Users/lspitler/prog/GitHub/huntsman-ms/resources/performance_ms.yaml'))
    imager = imagers['one_zwo_canon_full_moon']

    gband_limit = imager.point_source_limit(total_exp_time=exptime,
                                            sub_exp_time=exptime,
                                            filter_name='g',
                                            snr_target=1.0)

    # print(gband_limit)

    gband_limit_string = f"<{gband_limit.value}"

    star_table = get_star_table(width=imager.field_of_view[0],
                                height=imager.field_of_view[1],
                                mag_limit=gband_limit_string)
    # print(star_table)
    make_noiseless_image(imager, star_table)


if __name__ == '__main__':
    main()
