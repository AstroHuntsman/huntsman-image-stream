# huntsman-image-stream
Produce a python iterator of synthetic images.


## ImageStream()

Is a python generator class that does the following:

 - produce a fake noiseless image using gunagala of a given size
 - produces a noisey image based upon given exposure time

It will eventually:

 - start sending data no faster than the specified rate from: FPS, n_cameras, exptime
 - randomly inject transient bright sources
 - record processing speed (e.g. processing at X FPS)
 - randomly dim stars to mimic a occultation

## How to run it

Two ways:

```python
from generate_images import ImageStream
i = ImageStream()
next(i)
next(i)
```

or:

```python
from generate_images import ImageStream
for i in ImageStream():
   print(i)
```
