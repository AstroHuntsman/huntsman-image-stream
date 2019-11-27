import astropy.units as u

from image_util import generate_noiseless_image

from astropy.table import Table

from photutils.datasets import make_gaussian_sources_image


class ImageStream:

    # From https://treyhunner.com/2018/06/how-to-make-an-iterator-in-python/

    def __init__(self,
                 num_frames=10,
                 gunagala_config_filename='/Users/lspitler/Downloads/performance_detailed.yaml',
                 imager_filter_name='g',
                 fraction_of_field=1,
                 exptime=0.005 * u.s,
                 fps=20,
                 n_cameras=1,
                 flash_list=None,
                 ):
        self.exptime = exptime
        self.curr_frame_num = 0
        self.imager_filter_name = imager_filter_name
        self.gunagala_config_filename = gunagala_config_filename

        # prep flash info
        self.flash_info_dict = {}
        if flash_list is not None:
            for flash_info in flash_list:
                # (x_location, y_location, signal_to_noise, frame_no)
                self.flash_info_dict[flash_info[-1]] = flash_info

        print('Preparing noiseless image - this will take a minute or so')
        self.noiseless_image, self.imager = generate_noiseless_image(
            gunagala_config_filename=self.gunagala_config_filename,
            exptime=self.exptime,
            imager_filter_name=self.imager_filter_name)
        print('Created noiseless image - ready to run iterator')

    def __iter__(self):
        return self

    def __next__(self,
                 write_real_fits_file=False):

        self.curr_frame_num += 1

        if self.curr_frame_num in self.flash_info_dict.keys():
            # inject flash
            flash_image = generate_flash_image(self.imager,
                                               self.flash_info_dict[self.curr_frame_num])
            real_data = self.imager.make_image_real(self.noiseless_image + flash_image,
                                                    self.exptime)
        else:
            real_data = self.imager.make_image_real(self.noiseless_image,
                                                    self.exptime)
        return(real_data)
        if write_real_fits_file:
            real_data.write('out_real.fits', overwrite=True)

    def generate_flash_image(self, imager, flash_info):
        x_location, y_location, signal_to_noise, _ = flash_info

        flash_magnitude = imager.point_source_limit(total_exp_time=self.exptime,
                                                    sub_exp_time=self.exptime,
                                                    filter_name=self.imager_filter_name,
                                                    snr_target=signal_to_noise)
        flash_rate = self.imager.ABmag_to_rate(
            flash_magnitude, self.imager_filter_name)

        table = Table()
        table['amplitude'] = flash_rate.value
        table['x_mean'] = x_location
        table['y_mean'] = y_location
        # PSF width = https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/
        table['x_stddev'] = self.imager.psf.fwhm / 2.355
        table['y_stddev'] = self.imager.psf.fwhm / 2.355
        table['theta'] = 0

        flash_image = make_gaussian_sources_image((self.imager.wcs._naxis2,
                                                   self.imager.wcs._naxis1), table) * flash_rate.unit / u.pixel
        return(flash_image)


if __name__ == '__main__':
    """Test the stream generator."""
    import numpy as np

    # (x_location, y_location, signal_to_noise, frame_no)
    flash_list = [(100, 100, 1000, 5)]

    for i in ImageStream(num_frames=10, flash_list=flash_list):
        print(np.mean(i))
