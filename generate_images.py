import astropy.units as u

from image_util import generate_noiseless_image


class ImageStream:

    # From https://treyhunner.com/2018/06/how-to-make-an-iterator-in-python/

    def __init__(self,
                 gunagala_config_filename='/Users/lspitler/Downloads/performance_detailed.yaml',
                 fraction_of_field=1,
                 exptime=0.005 * u.s,
                 fps=20,
                 n_cameras=1
                 ):
        print('Preparing noiseless image - this will take a minute or so')
        self.noiseless_image, self.imager = generate_noiseless_image(
            gunagala_config_filename=gunagala_config_filename, exptime=exptime)
        self.exptime = exptime
        print('Created noiseless image - ready to run iterator')

    def __iter__(self):
        return self

    def __next__(self,
                 write_real_fits_file=False):

        real_data = self.imager.make_image_real(self.noiseless_image,
                                                self.exptime)

        return(real_data)
        if write_real_fits_file:
            real_data.write('out_real.fits', overwrite=True)
