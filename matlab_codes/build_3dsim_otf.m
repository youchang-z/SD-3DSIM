function [sim3d_otf] = build_3dsim_otf(paras)
% unit: nm 

x_pixelNumber = paras(1);
x_pixelSize = 40;
y_pixelNumber = paras(2);
y_pixelSize = 40;
z_pixelNumber = paras(3); 
z_pixelSize = paras(4);
exi_wavelength = paras(5);
emi_wavelength = paras(6);


%% illumination pattern simulation 
f0 = 1.1; % amplitude of the 0th-order beam at the objective lens back pupil
f1 = 1; % amplitudes of the four 1st-order beams at the objective lens back pupil

% set illumination angle
n = 1.4; % silicone oil: 1.4, water: 1.33, oil: 1.52, live cell: 1.38
theta = asin(exi_wavelength/(2*n*250)); % set illumination angle, 84um/166/2=250nm
phi = 0;

% set coordinates
x = x_pixelSize*ceil(-x_pixelNumber/2):x_pixelSize:x_pixelSize*ceil(x_pixelNumber/2-1);
y = y_pixelSize*ceil(-y_pixelNumber/2):y_pixelSize:y_pixelSize*ceil(y_pixelNumber/2-1);
Z = 0;
% building a 3D matrix is extremely computationally expensive
% so we set z=0 here to only simulate the illumination pattern of one z slice
% using projection-slice theorem, we can re-build the true 3D illumination pattern in the Fourier space
[X,Y] = meshgrid(x,y);

% define k vectors
% polarization effect is considered
k_amp = n*2*pi/exi_wavelength;
k_input = [0;0;1];
k_0 = k_amp*GJM(phi,0,k_input);
k_p10 = k_amp*GJM(phi,theta,k_input);
k_n10 = k_amp*GJM(phi+pi,theta,k_input);
k_0p1 = k_amp*GJM(phi+pi/2,theta,k_input);
k_0n1 = k_amp*GJM(phi+3*pi/2,theta,k_input);

% 0th beam polarization
e = 1; % ellipticity, e=0: s linear polarized at y, e=1: circularly polarized, e = infinity, p polarized
spin = -1; % clockwise: -1, counterclockwise: -1
px = spin*1i*1/sqrt(2)*sqrt(2*e^2/(1+e^2));
py = 1/sqrt(2) * sqrt(2/(1+e^2));
pz = 0;
p_input = [px;py;pz];
p_0 = GJM(phi,0,p_input);

% 4 1st beam polarization
e = 1; % ellipticity, e=0: s linear polarized at y, e=1: circularly polarized, e = infinity, p polarized
spin = -1; % clockwise: -1, counterclockwise: -1
px = spin*1i*1/sqrt(2)*sqrt(2*e^2/(1+e^2));
py = 1/sqrt(2) * sqrt(2/(1+e^2));
pz = 0;
p_input = [px;py;pz];
p_p10 = GJM(phi,theta,p_input);
p_n10 = GJM(phi+pi,theta,p_input);
p_0p1 = GJM(phi+pi/2,theta,p_input);
p_0n1 = GJM(phi+3*pi/2,theta,p_input);

% illumination pattern calculation
E_0_X = f0*exp(1i*k_0(1)*X + 1i*k_0(2)*Y + 1i*k_0(3)*Z)*p_0(1);
E_0_Y = f0*exp(1i*k_0(1)*X + 1i*k_0(2)*Y + 1i*k_0(3)*Z)*p_0(2);
E_0_Z = f0*exp(1i*k_0(1)*X + 1i*k_0(2)*Y + 1i*k_0(3)*Z)*p_0(3);
E_p10_X = f1*exp(1i*k_p10(1)*X + 1i*k_p10(2)*Y + 1i*k_p10(3)*Z)*p_p10(1);
E_p10_Y = f1*exp(1i*k_p10(1)*X + 1i*k_p10(2)*Y + 1i*k_p10(3)*Z)*p_p10(2);
E_p10_Z = f1*exp(1i*k_p10(1)*X + 1i*k_p10(2)*Y + 1i*k_p10(3)*Z)*p_p10(3);
E_n10_X = f1*exp(1i*k_n10(1)*X + 1i*k_n10(2)*Y + 1i*k_n10(3)*Z)*p_n10(1);
E_n10_Y = f1*exp(1i*k_n10(1)*X + 1i*k_n10(2)*Y + 1i*k_n10(3)*Z)*p_n10(2);
E_n10_Z = f1*exp(1i*k_n10(1)*X + 1i*k_n10(2)*Y + 1i*k_n10(3)*Z)*p_n10(3);
E_0p1_X = f1*exp(1i*k_0p1(1)*X + 1i*k_0p1(2)*Y + 1i*k_0p1(3)*Z)*p_0p1(1);
E_0p1_Y = f1*exp(1i*k_0p1(1)*X + 1i*k_0p1(2)*Y + 1i*k_0p1(3)*Z)*p_0p1(2);
E_0p1_Z = f1*exp(1i*k_0p1(1)*X + 1i*k_0p1(2)*Y + 1i*k_0p1(3)*Z)*p_0p1(3);
E_0n1_X = f1*exp(1i*k_0n1(1)*X + 1i*k_0n1(2)*Y + 1i*k_0n1(3)*Z)*p_0n1(1);
E_0n1_Y = f1*exp(1i*k_0n1(1)*X + 1i*k_0n1(2)*Y + 1i*k_0n1(3)*Z)*p_0n1(2);
E_0n1_Z = f1*exp(1i*k_0n1(1)*X + 1i*k_0n1(2)*Y + 1i*k_0n1(3)*Z)*p_0n1(3);
I = (E_0_X + E_p10_X + E_n10_X + E_0p1_X + E_0n1_X).*conj(E_0_X + E_p10_X + E_n10_X + E_0p1_X + E_0n1_X) + ...
    (E_0_Y + E_p10_Y + E_n10_Y + E_0p1_Y + E_0n1_Y).*conj(E_0_Y + E_p10_Y + E_n10_Y + E_0p1_Y + E_0n1_Y) + ...
    (E_0_Z + E_p10_Z + E_n10_Z + E_0p1_Z + E_0n1_Z).*conj(E_0_Z + E_p10_Z + E_n10_Z + E_0p1_Z + E_0n1_Z);

I = abs(I);
I = I./(max(I,[],'all'));
sim3d_ill_FD_singlePlane = fftshift(fft2(fftshift(I)));
sim3d_ill_FD_singlePlane = rescale(abs(sim3d_ill_FD_singlePlane));


%% simulate how 3D SIM shifts the wide-field OTF when the disk does not spin
max_shifted_freq_region = sim3d_ill_FD_singlePlane(...
    round(y_pixelNumber/2)-5:round(y_pixelNumber/2)+5,...
    ceil(x_pixelNumber/2+0.1)+20:end);
max_shifted_freq = max(max_shifted_freq_region,[],'all');
[row,col] = ind2sub(size(sim3d_ill_FD_singlePlane),...
                    find(sim3d_ill_FD_singlePlane == max_shifted_freq));
f_cut_xy = (max([row;col])-ceil(x_pixelNumber/2+0.1))*2;
f_cut_z = f_cut_xy/2*(1-cos(theta))/sin(theta);
ratio = (1/(z_pixelNumber*z_pixelSize))/(1/(x_pixelNumber*x_pixelSize));
f_cut_z = f_cut_z/ratio;

sim3d_otf_shift = zeros(x_pixelNumber,y_pixelNumber,z_pixelNumber);
step = round(f_cut_xy/4);
sim3d_otf_shift(:,:,ceil(size(sim3d_otf_shift,3)/2+0.1)) = sim3d_ill_FD_singlePlane;

sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)) = 0;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)) = 0;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)) = 0;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)) = 0;

sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)+round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+3*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step)/2;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)-round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+3*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step)/2;

sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)+round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+3*step)/2;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+3*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)-round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+3*step)/2;

sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)+round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step)/2;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)-round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,2)/2+0.1)+1*step)/2;

sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)+round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step)/2;
sim3d_otf_shift(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                ceil(size(sim3d_otf_shift,2)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step,...
                ceil(size(sim3d_otf_shift,3)/2+0.1)-round(f_cut_z)) = ...
                sim3d_ill_FD_singlePlane(ceil(size(sim3d_otf_shift,1)/2+0.1)-1*step:ceil(size(sim3d_otf_shift,1)/2+0.1)+1*step,...
                                         ceil(size(sim3d_otf_shift,2)/2+0.1)-3*step:ceil(size(sim3d_otf_shift,2)/2+0.1)-1*step)/2;

sim3d_otf_shift = sim3d_otf_shift.^2; % as the fluorescence pass through the disk one more time


%% 3d sim frequency shifts after rotating average (the disk spins)
page(1) = ceil(size(sim3d_otf_shift,3)/2+0.1);
page(2) = ceil(size(sim3d_otf_shift,3)/2+0.1) - round(f_cut_z);
page(3) = ceil(size(sim3d_otf_shift,3)/2+0.1) + round(f_cut_z);

distance_x = ceil(-x_pixelNumber/2):1:ceil(x_pixelNumber/2-1);
distance_y = ceil(-y_pixelNumber/2):1:ceil(y_pixelNumber/2-1);
[distance_X,distance_Y] = meshgrid(distance_x,distance_y);
distance_mtx = sqrt(distance_X.^2 + distance_Y.^2);

sim3d_otf_shift_rotavg = zeros(size(sim3d_otf_shift));
for k = 1:length(page)
    dict = zeros(1,round(max(distance_mtx,[],'all')+1));
    for i = 1:1:length(distance_x)
        for j = 1:1:length(distance_y)
            dict(round(distance_mtx(i,j))+1) = ...
            dict(round(distance_mtx(i,j))+1) + sim3d_otf_shift(i,j,page(k));
        end
    end
    count = 1:1:length(dict);
    dict = dict./(2*pi*count);
    for i = 1:1:length(distance_x)
        for j = 1:1:length(distance_y)
            sim3d_otf_shift_rotavg(i,j,page(k)) = dict(round(distance_mtx(i,j)+1));
        end
    end
end


%% build wide-field OTF (3D matrix)
wf_otf = build_wf_otf(paras);


%% build 3d sim OTF (3D matrix)
temp1 = ifftshift(ifftn(ifftshift(wf_otf)));
temp2 = ifftshift(ifftn(ifftshift(sim3d_otf_shift_rotavg)));
sim3d_otf = fftshift(fftn(fftshift(temp1.*conj(temp2))));
sim3d_otf = real(sim3d_otf);
sim3d_otf = rescale(sim3d_otf);
sim3d_otf(sim3d_otf<1e-5) = 0;
% imshow((log10(squeeze(sim3d_otf(:,ceil(y_pixelNumber/2+0.1),:))))',[])
% colormap('jet')


end