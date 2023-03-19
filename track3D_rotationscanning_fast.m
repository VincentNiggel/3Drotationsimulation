clear all
close all
clc
%% images folder's path 
path='C:\Users\vniggel\Documents\MyPhD\3Drotationsimulation\Paper3D\Simulation\Sim1';
cd(path)
list= dir(strcat(path,'\*'));
%% Definition of the parameters
Image_number=size(list,1)-2; %Number of images to analyze
image_size=50; % Final size of the images for computing the different angles
radius=30; %window radius applied to the image to crop
R_part=40/(radius+1)*image_size/2; %Radius of the particle in the image_size scale
mask_r=0.9; %proportion of the image with mask

gfilt=1; %size of the gaussianfilter from the matlab function imgaussfilt
contrastp=[0.2 0.9]; % contrast parameter from the matlab function imadjust

%definition of the different angle that will be scanned
ang_lim=16; %max scan angle
ang_step=2; %angle step 
angle_scanx=[-ang_lim:ang_step:ang_lim];
angle_scany=angle_scanx;
angle_scanz=angle_scanx;

%Apply to the first image to check the resized, cropped and adjusted image
    A= imread(list(1+2).name); %load first image of the folder
    [A_j,A_i]=meshgrid(1:size(A,2),1:size(A,1)); 
    centerA=(size(A,1)+1)/2; 
    crop_r=floor(radius+1);
    [B_jnew,B_inew]=meshgrid(linspace(centerA-crop_r,centerA+crop_r,image_size));
    
    image_B(1:image_size,1:image_size,1) = imadjust(imgaussfilt(interp2(A_j,A_i,im2double(A),B_jnew,B_inew),gfilt),contrastp);
    Radius_view=image_size/2-1;
    mask_radius=mask_r*Radius_view;
    centerB=(image_size+1)/2;  
    sizeB=image_size;
    [y_pos_init,x_pos_init]=meshgrid((1:sizeB)-centerB);

figure(1)
subplot(1,3,1)
imshow(A)
xlabel('Original image')
subplot(1,3,2)
imshow(image_B(:,:,1))
xlabel('Resized image with gaussian filter')
subplot(1,3,3)
imshow(image_B(:,:,1).*(x_pos_init.^2+y_pos_init.^2<mask_radius^2))
xlabel('Resized image with gaussian filter and mask')


%% Initialyze all wished rotation

 for ind=1:Image_number
      C= im2double(imread(list(ind+2).name));
      image_B(1:image_size,1:image_size,ind) = imadjust(imgaussfilt(interp2(A_j,A_i,C,B_jnew,B_inew),gfilt),contrastp);
 end



 % Create mask
%R_part=Radius_view/sin(angle_view);
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ; %rotation around axis x with t the angle in radian
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ; %rotation around axis y with t the angle in radian
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ; %rotion around axis z with t the angle in radian

centerB=(image_size+1)/2;  
sizeB=image_size;
[y_pos_init,x_pos_init]=meshgrid((1:sizeB)-centerB); %define the position of the pixels with the center of the image as origin

interesting_pixel_mask=(x_pos_init.^2+y_pos_init.^2)<(mask_radius)^2;%initial mask image
interesting_pixel=(x_pos_init.^2+y_pos_init.^2)<R_part^2;

z_pos_init=sqrt((R_part^2-(x_pos_init.^2+y_pos_init.^2)).*interesting_pixel);%definition of the position in z

pos_init(1,:)=x_pos_init(:).';
pos_init(2,:)=y_pos_init(:).';
pos_init(3,:)=z_pos_init(:).';

number_ang=size(angle_scanx,2)*size(angle_scany,2)*size(angle_scanz,2);% number of all angles set we will scan
%mask=zeros(sizeB,sizeB,number_ang);
x_rot1=zeros(sizeB,sizeB,number_ang);
y_rot1=zeros(sizeB,sizeB,number_ang);
ang_progress=zeros(3,number_ang);
ind_progress=1;
for anglex=angle_scanx
    for angley=angle_scany
        for anglez=angle_scanz
            pos_rot=rotz(anglez/180*pi)*roty(angley/180*pi)*rotx(anglex/180*pi)*pos_init; %position after rotation
            x_rot=reshape(pos_rot(1,:).',[sizeB,sizeB]);
            y_rot=reshape(pos_rot(2,:).',[sizeB,sizeB]);
            mask(1:sizeB,1:sizeB,ind_progress)=...
                (x_rot.^2+y_rot.^2 < mask_radius^2).*interesting_pixel_mask; %compute the mask associated to this case
            x_rot1(1:sizeB,1:sizeB,ind_progress)=x_rot;
            y_rot1(1:sizeB,1:sizeB,ind_progress)=y_rot;
            ang_progress(1:3,ind_progress)=[anglex angley anglez].'; %store the angles associtaed to this case
            ind_progress=ind_progress+1;
        end
    end
end
%%
imstep=1; %indices differences of the two images to compare

ind_image=20;
figure(2)
subplot(1,2,1)
imshow(image_B(:,:,ind_image));
subplot(1,2,2)
imshow(image_B(:,:,ind_image+imstep));

%% Normal

correlation_res=zeros(number_ang,Image_number-imstep);
    for ind_image=1:Image_number-imstep
        
        ind_image
        image_interest=image_B(1:sizeB,1:sizeB,ind_image+imstep);
        new_image=interp2(y_pos_init,x_pos_init,image_interest,y_rot1,x_rot1);%compute all theoritical rotated images
        new_image(isnan(new_image))=0;    
        correlation=0;
        image_initial=image_B(:,:,ind_image);
        ind_progress=1;
        for anglex=angle_scanx
            for angley=angle_scany
                for anglez=angle_scanz
                    new_image1=new_image(:,:,ind_progress);
                    correlation_res(ind_progress,ind_image)=corr2(image_initial(logical(mask(:,:,ind_progress))),new_image1(logical(mask(:,:,ind_progress))));
                    ind_progress=ind_progress+1;
                end
            end
        end
        
    end
  


[value_max ind_max]=max(correlation_res);
res_ang=ang_progress(1:3,ind_max);
%%
plot(res_ang(1,:))
hold on
plot(res_ang(2,:))
plot(res_ang(3,:))