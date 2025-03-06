# BagelFit 🥯
Fitting Torus Structures to Nuclear Membrane Density Maps  

BagelFit is a computational tool designed to fit toroidal structures onto nuclear membrane density maps using efficient search algorithms. Built for precision and scalability, it leverages integrative modeling techniques to optimize torus parameters and analyze nuclear morphology.


Features:

✅ Torus fitting with binary & non-binary density maps

✅ Cross-correlation-based optimization

✅ Dice coefficient calculation for overlap analysis

✅ Visualization of voxel intensity distributions

✅ Support for nuclear membrane structural studies

Whether you’re studying nuclear organization or just love bagels, BagelFit brings toroidal fitting to the next level! 🥯✨


## Documentation https://neeleshsoni21.github.io/BagelFit/

## Example Usages

```python
#-----------------------------------------------------#
# Example Usages
#-----------------------------------------------------#

import bagelfit as bf

#-----------------------------------------------------#
# Write a torus map using Torus parameters
#-----------------------------------------------------#
fitter = bf.BagelFitter()

tor_R = 670
tor_r = 160
tor_th = 85
extension = 0.0

best_torus = fitter.generate_binary_torus(
    tor_R, tor_r, tor_th, 
    extension=0.0, 
    boundingbox_length=2240, 
    voxel_size=10.0, 
    outmap_fname="./yeast_membrane/torus_yeast_fitted.mrc"
)   

#-----------------------------------------------------#
# Score a torus map with other torus maps or experimental maps
#-----------------------------------------------------#
fitter = bf.BagelFitter()

mapfile1 = "./yeast_membrane/Yeast_C8_Double_MR_center.mrc"
mapfile2 = "./yeast_membrane/torus_yeast_fitted.mrc"
fitter.score_torus_maps(mapfile1, mapfile2)

#-----------------------------------------------------#
# Fit several torus onto nuclear membrane 
# and write the highest overlapping torus into the file
#-----------------------------------------------------#
fitter = bf.BagelFitter()

# Below map has 11 million voxels
fitter.load_exprimental_map("./yeast_membrane/Yeast_C8_Double_MR_center.mrc")

tor_R_range = (670, 680, 10)
tor_r_range = (160, 180, 20)
tor_th_range = (85, 95, 10)
extension = 0.0

best_torus = fitter.fit_binary_torus(tor_R_range, tor_r_range, tor_th_range, extension)
# best_torus = fitter.fit_nonbinary_torus(tor_R_range, tor_r_range, tor_th_range, extension)

fitter.write_torusmap_to_file("./yeast_membrane/torus_yeast_fitted.mrc")