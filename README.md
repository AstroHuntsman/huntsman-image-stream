# huntsman-ms
 Millisecond astronomy tools


## generate_image_stream()

This will do the following:

 - produce a fake noiseless image using gunagala of a given size
 - calculate system noise given an exposure time
 - open a socket to send data
 - start sending data no faster than the specified rate (= exposure time by default)

 