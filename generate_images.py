"""Python generator synthetic image stream.
"""
import astropy.units as u

from image_util import generate_noiseless_image

from astropy.table import Table
from astropy.nddata import CCDData

from photutils.datasets import make_gaussian_sources_image


class ImageStream:

    """Python generator synthetic image stream.

    Attributes:
        curr_frame_num (int): current frame number
        exptime (astropy.Quantity): Exposure time in u.seconds.
        flash_info_dict (dict): Dictionary of flashes, frame number as key.
        gunagala_config_filename (str): absolute path to gunagala configuration file.
        imager_filter_name (str): Imaging filter name.
        num_frames (int): Number of frames to be generated.
    """

    # From https://treyhunner.com/2018/06/how-to-make-an-iterator-in-python/

    def __init__(self,
                 num_frames=10,
                 gunagala_config_filename='/Users/lspitler/Downloads/performance_detailed.yaml',
                 imager_filter_name='g',
                 exptime=0.005 * u.s,
                 flash_list=None,
                 snr_limit=3.
                 ):
        """Init function to setup the generator.

        Args:
            num_frames (int, optional): Number of frames to be generated.
            gunagala_config_filename (str, optional): absolute path to gunagala configuration file.
            imager_filter_name (str, optional): Imaging filter name.
            exptime (astropy.Quantity, optional): Exposure time in u.seconds.
            flash_list (list, optional): List of flash info tuples (x, y, S/N, frame #).
            snr_limit (float, optional): Include stars with luminosities S/N > snr_limit.
        """
        # following are not implemented yet
        #         fraction_of_field=1,
        #         fps=20,
        #         n_cameras=1,

        self.exptime = exptime
        self.curr_frame_num = 0
        self.imager_filter_name = imager_filter_name
        self.gunagala_config_filename = gunagala_config_filename
        self.num_frames = num_frames

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
            imager_filter_name=self.imager_filter_name,
            snr_limit=snr_limit)
        print('Created noiseless image - ready to run iterator')

    def __iter__(self):
        """Access point to new generator instance.

        Returns:
            ImageStream: New instance of ImageStream.
        """
        return self

    def generate_flash_image(self, imager, flash_info):
        """Generate a single noiseless image containing only a single flash.

        Args:
            imager (gunagala.Imager): Instance of a gunagala imaging system.
            flash_info (tuple): Tuple with flash info: (x, y, S/N, frame #).

        Returns:
            CCDData: Noiseless image containing a flash.
        """
        x_location, y_location, signal_to_noise, _ = flash_info

        flash_magnitude = imager.point_source_limit(total_exp_time=self.exptime,
                                                    sub_exp_time=self.exptime,
                                                    filter_name=self.imager_filter_name,
                                                    snr_target=signal_to_noise)
        flash_rate = self.imager.ABmag_to_rate(
            flash_magnitude, self.imager_filter_name)

        table = Table()
        table['amplitude'] = [flash_rate.value]
        table['x_mean'] = [x_location]
        table['y_mean'] = [y_location]
        # PSF width = https://brainder.org/2011/08/20/gaussian-kernels-convert-fwhm-to-sigma/
        table['x_stddev'] = [self.imager.psf.fwhm / 2.355]
        table['y_stddev'] = [self.imager.psf.fwhm / 2.355]
        table['theta'] = [0]

        naxis1, naxis2 = self.imager.wcs.pixel_shape
        flash_image = make_gaussian_sources_image(
            (naxis2, naxis1), table) * flash_rate.unit / u.pixel

        return(CCDData(flash_image, wcs=self.imager.wcs))

    def __next__(self):
        """Main python generator routine. Returns
        a noisey image with a flash if appropriate.

        Returns:
            CCDData: Noisy imaging data.

        Raises:
            StopIteration: When current frame is greater than requested frame numbers, stop loop.
        """
        self.curr_frame_num += 1

        if self.curr_frame_num > self.num_frames:
            raise StopIteration

        if self.curr_frame_num in self.flash_info_dict.keys():
            # inject flash
            flash_image = self.generate_flash_image(self.imager,
                                                    self.flash_info_dict[self.curr_frame_num])
            flash_image.data = flash_image.data + self.noiseless_image.data
            real_data = self.imager.make_image_real(flash_image,
                                                    self.exptime)
        else:
            real_data = self.imager.make_image_real(self.noiseless_image,
                                                    self.exptime)
        return(real_data)


if __name__ == '__main__':
    """Test the stream generator. Produce 1 flash."""

    # (x_location, y_location, signal_to_noise, frame_no)
    flash_list = [(100, 100, 1000, 3), (1000, 1000, 1000, 2)]

    for i, image in enumerate(ImageStream(num_frames=5, flash_list=flash_list)):
        image.write(f'out_synth_{i+1}.fits', overwrite=True)
