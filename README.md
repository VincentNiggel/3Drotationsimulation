# 3Drotationsimulation
These codes are related to the publication: "3-D rotation tracking from 2-D images of spherical colloids with textured surface" by V. Niggel et al., published in Soft Matter in 2023.
makenewsimparticle_final.m allows to create simulated images of the rotated raspberry particle while track3D_rotationscanning.m is the code used to analyze the result in the paper. 
track3D_rotationscanning_fast.m is a new version which will only compare the pixels inside the mask resulting in a faster code but slightly different comparison result between two images.
track3D_rotationscanning_parallel.m is a paralyzed version of the code makenewsimparticle_final.m.

The different python code are translation attempts of the matlab codes.
