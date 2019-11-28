# huntsman-image-stream
Produce a python iterator of synthetic images.


## ImageStream()

Is a python generator class that does the following:

 - produce a fake noiseless image for a specific setup: [Canon lens](https://store.canon.com.au/lenses/ef-400mm-f-2-8l-is-iii-usm.html) + [ZWO ASI183 Pro Camera](https://astronomy-imaging-camera.com/product/asi183mm-pro-mono).
 - produces a noisey image based upon given exposure time
 - inject bright sources of light

It will eventually:

 - start sending data no faster than the specified rate from: FPS, n_cameras, exptime
 - record processing speed (e.g. processing at X FPS)
 - randomly dim stars to mimic a occultation

## Requirements

A few python modules, including a custom astronomy sensitivity suite [`gunagala`](https://github.com/AstroHuntsman/gunagala).

To install, just run pip in source directory:

```
pip install -r requirements.txt
```

You can test it out by going to the source directory and running these:
```
python image_util.py absolute_path_to_gunagala_configuration_file
python generate_images.py absolute_path_to_gunagala_configuration_file
```

## How to run it

Two ways:

```python
from generate_images import ImageStream
i = ImageStream(gunagala_config_filename=absolute_path_to_gunagala_configuration_file)
next(i)
next(i)
```

or:

```python
from generate_images import ImageStream
for i in ImageStream(gunagala_config_filename=absolute_path_to_gunagala_configuration_file):
   print(i)
```

Also check out the `if __name__ == '__main__':` example in `generate_images.py`.


