clear all
close all
clc
%% parameters
imsize=100;
angle_view=85/180*pi;
increase_iso=2;
N_berries=150;
gaussfilt_berr=3;
Radius_view=0.4*imsize;
R_part=Radius_view/sin(angle_view);
R_berries=0.13*R_part;
std_berries=0.1*R_berries;
ang_rot_limit=10;
Nb_image=500;
Max_ring=72.15;
%% Create ring
[y_image,x_image]=meshgrid((1:imsize)-(imsize+1)/2);
fx1=@(t) Max_ring/255*exp(-((t-Radius_view)/1.603/Radius_view*15.05).^2);
fx2=@(t) 0.0878*exp(-((t-Radius_view)/1.885/Radius_view*13.78).^2);
pixel_distance=sqrt(x_image.*x_image+y_image.*y_image);
ring_intensity=fx1(pixel_distance)*255+...
            normrnd(zeros(imsize,imsize,Nb_image),ones(imsize,imsize,Nb_image).*fx2(pixel_distance)*255);
ring_intensity(ring_intensity<0)=0;
ring_intensity(ring_intensity>130)=130;
% Create mask for value
mask=(x_image.*x_image+y_image.*y_image<Radius_view*Radius_view);
mask_value=(sqrt(abs(Radius_view^2-pixel_distance.*pixel_distance))...
            /Radius_view*(1-Max_ring/255)+Max_ring/255).*mask;
%%
%taken from Markus Deserno paper: How to generate equidistributed points on
%the surface of a sphere, https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf

N=ceil(Radius_view^2*pi*4/sin(angle_view)^2*increase_iso);
r=1;
rot = @(t,u) r*[sin(t)*cos(u);sin(t)*sin(u);cos(t)] ;
N_count=0;
a=4*pi*r^2/N;
d=sqrt(a);
M_v=round(pi/d);
d_v=pi/M_v;
d_p=a/d_v;
for m=0:M_v-1
    v=pi*(m+0.5)/M_v;
    M_p=round(2*pi*sin(v)/d_p);
    for n=0:M_p-1
        p=2*pi*n/M_p;
        pos_iso_init(:,N_count+1)=rot(v,p);
        N_count=N_count+1;
    end
end
Surface_im=zeros(1,size(pos_iso_init,2));
pos_iso_init=pos_iso_init*R_part;
%%
rotx = @(t) [1 0 0; 0 cos(t) -sin(t) ; 0 sin(t) cos(t)] ;
roty = @(t) [cos(t) 0 sin(t) ; 0 1 0 ; -sin(t) 0  cos(t)] ;
rotz = @(t) [cos(t) -sin(t) 0 ; sin(t) cos(t) 0 ; 0 0 1] ;

N=5000;
r=1;
rot = @(t,u) r*[sin(t)*cos(u);sin(t)*sin(u);cos(t)] ;
N_count=0;
a=4*pi*r^2/N;
d=sqrt(a);
M_v=round(pi/d);
d_v=pi/M_v;
d_p=a/d_v;
for m=0:M_v-1
    v=pi*(m+0.5)/M_v;
    M_p=round(2*pi*sin(v)/d_p);
    for n=0:M_p-1
        p=2*pi*n/M_p;
        pos_iso_init2(:,N_count+1)=rot(v,p);
        N_count=N_count+1;
    end
end
Count_berr=zeros(1,size(pos_iso_init2,2));
count=0;
ang_interest=zeros(3,size(pos_iso_init2,2));
while count<size(pos_iso_init2,2)
    ang_rot_init=rand(1,3)*2*pi;
    pos_iso1=rotz(ang_rot_init(1,3))*roty(ang_rot_init(1,2))*rotx(ang_rot_init(1,1))*[0;0;1];
    dii=(pos_iso1(1,1)-pos_iso_init2(1,:)).^2+(pos_iso1(2,1)-pos_iso_init2(2,:)).^2+(pos_iso1(3,1)-pos_iso_init2(3,:)).^2;
    [ind2,ind3]=min(dii);
    if Count_berr(1,ind3)==0
        count=count+1
        Count_berr(1,ind3)=Count_berr(1,ind3)+1;
        ang_interest(1:3,count)=ang_rot_init;
    end
end
    
for ind=1:10
    ang_interest=ang_interest(:,randperm(count));
end

%%
%create an ellipse

noisex=normrnd(0,std_berries, [1,N_berries]); %associate each berry with a noise term for its shape
noisey=normrnd(0,std_berries, [1,N_berries]);
berr_ang=rand(1,N_berries)*pi;
fx3= @(t,u,sigmax,sigmay) exp(-((t/sigmax).^2/2+(u/sigmay).^2/2));

ang_rot_ind=ceil(rand(1,N_berries)*size(ang_interest,2));
for ind = 1:N_berries %generate each mask indiviudally
    ang_rot_init=ang_interest(1:3,ang_rot_ind(1,ind));
    pos_iso=rotz(ang_rot_init(3,1))*roty(ang_rot_init(2,1))*rotx(ang_rot_init(1,1))*pos_iso_init;
    radiusX = R_berries+ noisex(ind);
    radiusY = R_berries+ noisey(ind);
    Iso_berr_pos(1,:)=pos_iso(1,:)*cos(berr_ang(ind))+pos_iso(2,:)*sin(berr_ang(ind));
    Iso_berr_pos(2,:)=pos_iso(1,:)*sin(berr_ang(ind))-pos_iso(2,:)*cos(berr_ang(ind));
    Iso_berr_value=(Iso_berr_pos(1,:)).^2 ./ radiusY^2 ...
                 + (Iso_berr_pos(2,:)).^2 ./ radiusX^2 ;
    Iso_berr_ellipse=(Iso_berr_value<=1);
    Iso_berr_value=fx3(Iso_berr_pos(1,:),Iso_berr_pos(2,:),radiusY*0.7,radiusX*0.7)*200/255+rand(1,size(pos_iso_init,2))*(255-200)/255;
    Iso_berr_value=Iso_berr_value.*Iso_berr_ellipse;
    Iso_berr_value=Iso_berr_value.*(pos_iso(3,:)>0);
    Surface_im=max(Surface_im,Iso_berr_value);  
end

%%
x_image=x_image.*mask;
y_image=y_image.*mask;
pos_iso=pos_iso_init;
ang_rot_sim=(rand(3,Nb_image)-0.5)*2*ang_rot_limit/180*pi; %define the rotation angles  between each images
%% Creates the 2D images
im_rotation=zeros(imsize,imsize,Nb_image);
for ind_image=1:Nb_image
    ind_image
    interesting_berries=pos_iso.*(pos_iso(3,:)>R_part*cos(angle_view+5/180*pi));
    Interest_pos=interesting_berries;
    Interest_pos(:,~any(Interest_pos,1)) = [];
    Interest_value=Surface_im;
    Interest_value(:,~any(interesting_berries,1)) = [];
    im_rotation_inter=scatteredInterpolant(Interest_pos(2,:).',Interest_pos(1,:).',Interest_value.');
    im_rotation(1:imsize,1:imsize,ind_image)=imgaussfilt(max(im_rotation_inter(y_image,x_image).*...
                                             mask_value*255,ring_intensity(:,:,ind_image)),0.7);
    pos_iso=rotz(ang_rot_sim(3,ind_image))*roty(ang_rot_sim(2,ind_image))*rotx(ang_rot_sim(1,ind_image))*pos_iso;
end
%% Save the created images
for ind_image=1:Nb_image
 imwrite(uint8(im_rotation(1:imsize,1:imsize,ind_image)),strcat('Vincent_simulations_ellipse_',num2str(ind_image,'%04d'),'.tif'));
end
%% To verify how the image looks like
figure(1)
imshow(uint8(imgaussfilt(im_rotation(:,:,1),0.7)))

