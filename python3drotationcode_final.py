# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 15:37:05 2022

@author: vniggel
"""

import os
import numpy as np
#import cv2
import matplotlib.pyplot as plt
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image, ImageEnhance
from scipy.ndimage import gaussian_filter
from scipy import ndimage

def imadjustgray(x,a,b,c,d,gamma=1):
    # Similar to imadjust in MATLAB.
    # Converts an image range from [a,b] to [c,d].
    # The Equation of a line can be used for this transformation:
    #   y=((d-c)/(b-a))*(x-a)+c
    # However, it is better to use a more generalized equation:
    #   y=((x-a)/(b-a))^gamma*(d-c)+c
    # If gamma is equal to 1, then the line equation is used.
    # When gamma is not equal to 1, then the transformation is not linear.

    y = (((x - a) / (b - a)) ** gamma) * (d - c) + c
    return y
#%%
path="C:\\Users\\vniggel\\Documents\\MyPhD\\3Drotationsimulation\\Paper3D\\Simulation\\Sim1"
os.chdir(path)
list1=os.listdir(path)

#%%
Image_number=len(list1)
image_size=50
radius=30
R_part=40/(radius+1)*image_size/2
mask_r=0.9

gfilt=0.1
contrastp=np.array([0.2,1])*255
contraste=np.array([0,1])*255
ang_lim=16
ang_step=2 
angle_scanx=np.arange(-ang_lim, ang_lim+ang_step, ang_step, dtype=int)
angle_scany=angle_scanx
angle_scanz=angle_scanx

A= np.array(Image.open(list1[0]))
centerA=(A.shape[0]+1)/2
crop_r=math.floor(radius+1)

B_jnew,B_inew=np.meshgrid(np.linspace(centerA-crop_r,centerA+crop_r,image_size),np.linspace(centerA-crop_r,centerA+crop_r,image_size))
A1=ndimage.map_coordinates(A.astype(np.float32), [B_inew.ravel()-1, B_jnew.ravel()-1], order=1).reshape(B_inew.shape)
A1=imadjustgray(gaussian_filter(A1, sigma=gfilt),contrastp[0],contrastp[1],contraste[0],contraste[1])


Radius_view=image_size/2-1
mask_radius=mask_r*Radius_view
centerB=(image_size-1)/2
sizeB=image_size
y_pos_init,x_pos_init=np.meshgrid(np.linspace(-centerB,sizeB-centerB-1,sizeB),np.linspace(-centerB,sizeB-centerB-1,sizeB))
VV=(np.square(x_pos_init)+np.square(y_pos_init))<mask_radius**2
f=plt.figure()
f.add_subplot(1, 3, 1)
plt.imshow(A, cmap='gray', vmin=0, vmax=255)
f.add_subplot(1, 3, 2)
plt.imshow(A1, cmap='gray', vmin=0, vmax=255)
f.add_subplot(1, 3, 3)
plt.imshow(A1*VV, cmap='gray', vmin=0, vmax=255)
plt.show()


#%% Initialyze all wished rotation
image_B=[]
for ind in range(Image_number):
      A= np.array(Image.open(list1[ind]))
      A1=ndimage.map_coordinates(A.astype(np.float32), [B_inew.ravel()-1, B_jnew.ravel()-1], order=1).reshape(B_inew.shape)
      A1=imadjustgray(gaussian_filter(A1, sigma=gfilt),contrastp[0],contrastp[1],contraste[0],contraste[1])
      image_B.append(A1)
      
#%%
def rotx(t):
    y=np.array([[1,0,0],[0, math.cos(t), -math.sin(t)],[0, math.sin(t), math.cos(t)]])
    return y
def roty(t):
    y=np.array([[math.cos(t),0,math.sin(t)],[0, 1, 0],[-math.sin(t), 0, math.cos(t)]])
    return y
def rotz(t):
    y=np.array([[math.cos(t), -math.sin(t), 0],[math.sin(t), math.cos(t), 0],[0, 0, 1]])
    return y

#%%
interesting_pixel_mask=x_pos_init**2+y_pos_init**2<mask_radius**2
interesting_pixel=(np.square(x_pos_init)+np.square(y_pos_init))<np.square(R_part)
z_pos_init=np.sqrt((R_part**2-(x_pos_init**2+y_pos_init**2))*interesting_pixel)
pos_init=np.zeros((3, x_pos_init.size))
for ind in range(x_pos_init.shape[1]):
    pos_init[[0],ind*x_pos_init.shape[1]:(ind+1)*x_pos_init.shape[1]]=np.transpose(x_pos_init[0:x_pos_init.shape[1],ind])
    pos_init[[1],ind*x_pos_init.shape[1]:(ind+1)*x_pos_init.shape[1]]=np.transpose(y_pos_init[0:x_pos_init.shape[1],ind])
    pos_init[[2],ind*x_pos_init.shape[1]:(ind+1)*x_pos_init.shape[1]]=np.transpose(z_pos_init[0:x_pos_init.shape[1],ind])

number_ang=angle_scanx.shape[0]*angle_scany.shape[0]*angle_scanz.shape[0]
ang_progress=np.zeros((3,number_ang));
mask=[]
x_rot2=np.zeros((sizeB,angle_scanx.shape[0]*angle_scany.shape[0]*angle_scanz.shape[0]*sizeB));
y_rot2=np.zeros((sizeB,angle_scanx.shape[0]*angle_scany.shape[0]*angle_scanz.shape[0]*sizeB));
ind_progress=0;
for anglex in angle_scanx:
    for angley in angle_scany:
        for anglez in angle_scanz:
            pos_rot=np.matmul(rotz(anglez/180*math.pi),np.matmul(roty(angley/180*math.pi),np.matmul(rotx(anglex/180*math.pi),pos_init)))
            x_rot=np.transpose(np.reshape(pos_rot[0,0:x_pos_init.size],(sizeB,sizeB)))
            y_rot=np.transpose(np.reshape(pos_rot[1,0:x_pos_init.size],(sizeB,sizeB)))
            mask.append(np.ones((sizeB,sizeB))*(x_rot**2+y_rot**2<mask_radius**2)*interesting_pixel_mask)
            x_rot2[0:sizeB,ind_progress*sizeB:(ind_progress+1)*sizeB]=x_rot
            y_rot2[0:sizeB,ind_progress*sizeB:(ind_progress+1)*sizeB]=y_rot
            ang_progress[0:3,[ind_progress]]=np.array([[anglex],[angley],[anglez]])
            ind_progress=ind_progress+1
            
#%%
imstep=1
ind_image=19;

f=plt.figure()
f.add_subplot(1, 2, 1)
plt.imshow(image_B[ind_image], cmap='gray', vmin=0, vmax=255)
f.add_subplot(1, 2, 2)
plt.imshow(image_B[ind_image+imstep], cmap='gray', vmin=0, vmax=255)
plt.show()    
#%%
def mean2(x):
    y = np.sum(x) / np.size(x);
    return y

def corr2(a,b):
    a = a - mean2(a)
    b = b - mean2(b)
    r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum());
    return r

correlation_res=np.zeros((number_ang,Image_number-imstep))
#for ind_image in range(Image_number-imstep):
for ind_image in range(5):
    print(ind_image+1)
    image_interest=image_B[ind_image+imstep]
    new_image=ndimage.map_coordinates(image_interest, [x_rot2+centerB, y_rot2+centerB], order=1)
    correlation=0
    image_initial=image_B[ind_image]
    ind_progress=0
    for anglex in angle_scanx:
        for angley in angle_scany:    
            for anglez in angle_scanz:
                correlation_res[ind_progress,ind_image]=corr2(image_initial*mask[ind_progress],new_image[0:sizeB,ind_progress*sizeB:(ind_progress+1)*sizeB]*mask[ind_progress])
                ind_progress=ind_progress+1
ind_max=np.argmax(correlation_res, axis=0)
res_ang=ang_progress[0:3,np.argmax(correlation_res, axis=0)]
#%%
plt.plot(np.arange(50)[np.newaxis, :],res_ang[[0],0:50],'ro',np.arange(50)[np.newaxis, :],res_ang[[1],0:50],'go',np.arange(50)[np.newaxis, :],res_ang[[2],0:50],'bo')
plt.show()    
UU=res_ang[[0],0:]