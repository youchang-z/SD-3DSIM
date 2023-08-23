function wf_otf = build_wf_otf(paras)
% unit: nm 

x_pixelNumber = paras(1);
x_pixelSize = 40;
y_pixelNumber = paras(2);
y_pixelSize = 40;
z_pixelNumber = paras(3);
z_pixelSize = paras(4);
exi_wavelength = paras(5);
emi_wavelength = paras(6);

fx_step = 1/(x_pixelSize*2*x_pixelNumber/2);
fy_step = 1/(y_pixelSize*2*y_pixelNumber/2);
fz_step = 1/(z_pixelSize*2*z_pixelNumber/2);

fx = fx_step*ceil(-x_pixelNumber/2):fx_step:fx_step*ceil(x_pixelNumber/2-1);
fy = fy_step*ceil(-y_pixelNumber/2):fy_step:fy_step*ceil(y_pixelNumber/2-1);
fz = fz_step*ceil(-z_pixelNumber/2):fz_step:fz_step*ceil(z_pixelNumber/2-1);

n = 1.4; % silicone oil: 1.4, water: 1.33, oil: 1.52, live cell: 1.38
fcut_ampValue = n/emi_wavelength; % the radius of the amplitude edwald sphere
fx = fx./fcut_ampValue;
fy = fy./fcut_ampValue;
fz = fz./fcut_ampValue;

[Fx,Fy,Fz] = meshgrid(fx,fy,fz);
index = Fz == 0;
Fz(index) = 1e-10;

l = sqrt(Fx.^2 + Fy.^2);
K = sqrt(l.^2 + Fz.^2);
p = 2*abs(l).*sqrt(1-K.^2/4)./(K.*abs(Fz));
alpha = 72/180*pi;
otf = 4 * acos( (2*cos(alpha)./abs(Fz)+1)./p ) ./(pi*K);
otf = abs(otf);

index = (1-K.^2/4) < 0;
otf(index) = 0;
index = log10((2*cos(alpha)./abs(Fz)+1)./p) > 0;
otf(index) = 0;

wf_otf = rescale(otf);

wf_otf_xy_proj = sum(wf_otf,3);
ratio = wf_otf_xy_proj(ceil(size(otf,1)/2+0.1)+1,ceil(size(otf,2)/2+0.1)) /...
        wf_otf_xy_proj(ceil(size(otf,1)/2+0.1)+2,ceil(size(otf,2)/2+0.1));
center_value = wf_otf_xy_proj(ceil(size(otf,1)/2+0.1)+1,ceil(size(otf,2)/2+0.1)) * ratio;
wf_otf(ceil(size(otf,1)/2+0.1),ceil(size(otf,2)/2+0.1),ceil(size(otf,3)/2+0.1)) = center_value;

wf_otf = rescale(wf_otf);


end