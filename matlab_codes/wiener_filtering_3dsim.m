clear all; close all; clc;

folder_path = '';
input_file_name = '';
zSpacing = 100; % 40nm for the beads, 100nm for the myocytes
exi = 488;
emi = 532; % 515 for beads, 532 for myocytes 


%% load image stack
info = imfinfo([folder_path input_file_name '.tif']);
num_images = numel(info);
for k = 1:num_images
    A(:,:,k) = im2double(imread([folder_path input_file_name '.tif'],k));
end
paras = [2048 2048 num_images zSpacing exi emi];
bk = 110; % camera thermal noise
A = A - bk/65535;
A(A<0) = 0;


%% calibrate intensity variation caused by photobleaching
int_tot = zeros(1,num_images);
for k = 1:num_images
    int_tot(k) = sum(A(:,:,k),"all");
end
x = (1:1:length(int_tot))';
y = int_tot';

fluo_decay = fittype('a + b1*exp(c1*(x+t)) + b2*exp(c2*(x+t))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a','b1','c1','b2','c2','t'});
options = fitoptions(fluo_decay);
options.Upper = [bk/65535*2048*2048 int_tot(1)*2 0 int_tot(1) 0 1000];
options.Lower = [0 0 -1 0 -1 0];
options.StartPoint = [0.5*bk/65535*2048*2048 int_tot(1) -0.01 int_tot(1) -0.01 10];
fit_result = fit(x,y,fluo_decay,options);
bk = fit_result.a/(2048*2048)*65535;

% remove bk
A = A - bk/65535;
A(A<0) = 0;

% correct bleaching 
bleaching_calib = 1./(fit_result(x)-fit_result.a);
for k = 1:num_images
    A(:,:,k) = A(:,:,k)*bleaching_calib(k);
end
A = A/max(A,[],'all');


%% avoid edge problem on Z
add_num = 10;
top = repmat(A(:,:,1),1,1,add_num);
bot = repmat(A(:,:,end),1,1,add_num);
A = cat(3,top,A,bot);
num_images = num_images + add_num*2;
paras(3) = paras(3) + add_num*2;


%% build frequency coordinates
x_pixelNumber = paras(1);
x_pixelSize = 40;
y_pixelNumber = paras(2);
y_pixelSize = 40;
z_pixelNumber = paras(3); 
z_pixelSize = paras(4);
exi_wavelength = paras(5);
emi_wavelength = paras(6);
n = 1.4;

fx_step = 1/(x_pixelSize*2*x_pixelNumber/2);
fy_step = 1/(y_pixelSize*2*y_pixelNumber/2);
fz_step = 1/(z_pixelSize*2*z_pixelNumber/2);
fx = fx_step*ceil(-x_pixelNumber/2):fx_step:fx_step*ceil(x_pixelNumber/2-1);
fy = fy_step*ceil(-y_pixelNumber/2):fy_step:fy_step*ceil(y_pixelNumber/2-1);
fz = fz_step*ceil(-z_pixelNumber/2):fz_step:fz_step*ceil(z_pixelNumber/2-1);
fcut_ampValue = 2*n/emi_wavelength; % the radius of the amplitude edwald sphere
fx = fx./fcut_ampValue;
fy = fy./fcut_ampValue;
fz = fz./fcut_ampValue;
[Fx,Fy,Fz] = meshgrid(fx,fy,fz);


%% compute FD of the input image stack
image_FD = fftshift(fftn(fftshift(A)));
powerSpec = abs(image_FD).^2;
figure()
imshow(squeeze(log10(powerSpec(:,ceil(size(powerSpec,2)/2+0.1),:)))',[]);


%% compute 3d sim OTF
[otf] = build_3dsim_otf(paras);
% figure()
% imshow(squeeze(log10(otf(:,ceil(size(otf,2)/2+0.1),:)))',[]);


%% suppress lateral low freq to aviod artifacts
strip = exp(-(Fx.^2+Fy.^2)/0.0002);
strip = 1- strip + 1e-1;
% figure()
% plot(squeeze((strip(:,ceil(size(strip,2)/2+0.1),ceil(size(strip,3)/2+0.1)))))


%% build signal and noise mask for the OTF
index = otf > 1e-10;
signal_mask = zeros(size(otf));
signal_mask(index) = 1;
signal_mask = logical(signal_mask);
noise_mask = ones(size(otf)) - signal_mask;
noise_mask = logical(noise_mask);


%% calculate power spectral density of signal (a scaler)
temp = signal_mask.*powerSpec;
temp = reshape(temp,1,size(temp,1)*size(temp,2)*size(temp,3));
index = temp == 0;
temp(index) = [];
sig_avg = mean(temp,'all');


%% calculate power spectral density of noise (a scaler)
temp = noise_mask.*powerSpec;
temp = reshape(temp,1,size(temp,1)*size(temp,2)*size(temp,3));
index = temp == 0;
temp(index) = [];
noise_avg = mean(temp,'all');


%% + apodization to avoid negative rings around PSF
apo = exp(-(Fx.^2+Fy.^2+8*Fz.^2)/4);
% figure()
% plot(squeeze((apo(:,ceil(size(strip,2)/2+0.1),ceil(size(strip,3)/2+0.1)))))


%% calculate wiener filter
w_mtx = (noise_avg./sig_avg)^1.45.*ones(size(image_FD)); 
% 1.7 for beads, 1.45 for myocytes
% this value should be around 1-2, but the exactly value should be
% adjusted based on SNR. Lower SNR -> smaller value.
wienerF = strip.*apo.*(otf)./(otf.^2 + w_mtx);
wienerF = rescale(wienerF,1,50);

% figure()
% subplot(211)
% imshow(squeeze((abs(wienerF(:,ceil(size(powerSpec,2)/2+0.1),:))))',[]);
% subplot(212)
% line = wienerF(:,ceil(paras(2)/2+0.1),ceil(paras(3)/2+0.1));
% plot(line);
% title('x axis')


%% restore the image by iterative inverse filtering
% the number of iteration is adjusted empirically
for i = 1:1:5
    if i == 1
        image_restored_FD = image_FD.*wienerF;
    else
        image_restored_FD = image_restored_FD.*wienerF;
    end

    image_restored = ifftshift(ifftn(ifftshift(image_restored_FD)));
    image_restored = real(image_restored);
    index = image_restored < 0;
    image_restored(index) = 0;
    image_restored = image_restored./max(image_restored,[],'all');
    image_restored_FD = fftshift(fftn(fftshift(image_restored)));

%     stepsize = 0.7^(i); % adjusted empirically
    image_restored_FD = (image_restored_FD + ...
        stepsize*(image_FD.*signal_mask - image_restored_FD.*otf));

end

image_restored = ifftshift(ifftn(ifftshift(image_restored_FD)));
image_restored = real(image_restored);
image_restored = 0.99*image_restored./max(image_restored,[],'all');

% figure()
% temp = (abs(fftshift(fftn(fftshift(image_restored)))));
% subplot(121)
% imshow(squeeze(log10(abs(temp(:,ceil(size(powerSpec,2)/2+0.1),:))))',[]);
% subplot(122)
% imshow(squeeze(log10(abs(temp(:,:,ceil(size(powerSpec,3)/2+0.1)))))',[]);


%% output image
output = image_restored./max(image_restored,[],'all');
output = im2uint16(output);
num_images = size(image_restored,3);
for k = 1+add_num:num_images-add_num
    imwrite(output(:,:,k),[folder_path input_file_name '_processed' '.tif'],...
        "tif",'WriteMode','append');
end




