from PIL import Image
import numpy as np
import h5py



hf = h5py.File('image_test_data.h5', 'w')

img = Image.open('ladybug.jpg')
ar = np.array(img)

#ar[height,width,R G or B]

hf.create_dataset('rgb', data=ar)


hf.close()


