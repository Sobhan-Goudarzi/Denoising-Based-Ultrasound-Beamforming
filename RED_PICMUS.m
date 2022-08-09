% Author: Sobhan Goudarzi
% IMage Processing And Characterization of Tissue (IMPACT) Group
% Concordia University
% email address: sobhan.goudarzi@concordia.ca
% August 2022
clc
clear all
close all
warning off
%% data type
acquisition = 'simulation';       % simulation || experiments || in_vivo
phantom = 'resolution_distorsion';      % resolution_distorsion || contrast_speckle || carotid_cross || carotid_long
suffix = 'simu';                   % simu || expe
path_dataset = ['./PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_dataset_rf','.hdf5'];
path_scan = ['./PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_scan','.hdf5'];
addpath('./minFunc_2012')% update the folder directory
%% loading dataset
%-- load scan and dataset
scan_f = linear_scan();
scan_f.read_file(path_scan);
dataset = us_dataset();
dataset.read_file(path_dataset);
%-- define scan based on time axis
time = (150:1400).'/dataset.sampling_frequency+dataset.initial_time;
z_axis= time*dataset.c0/2;
scan = linear_scan(scan_f.x_axis,z_axis);
%% Image reconstruction using DAS
%-- receive apodization
%-- dynamically expanding receive aperture with hanning apodization
rx_f_number = 0.5;
rx_aperture = scan.z/rx_f_number;
rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'hanning');
%-- angular apodization
angular_apodization = ones(scan.pixels,dataset.firings);
%-- beamforming loop
beamformed_data = zeros(scan.pixels,1);
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
w0 = 2*pi*dataset.modulation_frequency;
wb = waitbar(0,'DAS beamforming');
%
tgc1 = time_vector'./max(time_vector);
tgc1 = exp(4*tgc1);
%
for pw=38
    %-- transmit delay
    transmit_delay = scan.z*cos(dataset.angles(pw))+scan.x*sin(dataset.angles(pw));
    for nrx=1:dataset.channels
        %-- progress bar
        step=(nrx + (pw-1)*dataset.channels)/length(dataset.angles)/dataset.channels;
        waitbar(step,wb,sprintf('DAS-RF beamforming %0.0f%%',step*100));
        %-- receive delay
        receive_delay = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
        %-- total delay
        delay = (transmit_delay+receive_delay)/dataset.c0;
        %-- phase shift
        phase_shift = exp(1i.*w0*(delay-2*scan.z/dataset.c0));
        %-- beamformed data
        beamformed_data = beamformed_data + phase_shift.*angular_apodization(:,pw).*receive_apodization(:,nrx).*interp1(time_vector,tgc1.*dataset.data(:,nrx,pw),delay,'spline',0);
    end
end
close(wb);
beamformed_data(isnan(beamformed_data))=0;
DAS_RF = beamformed_data;
DAS_RF = DAS_RF./max(abs(DAS_RF));
DAS_RF = double(DAS_RF);
DAS_ENV = abs(hilbert(reshape(DAS_RF,[scan.Nz scan.Nx])));
%% Inverse beamforming step
load('Phi.mat'); %loading the weighting matrix
y = dataset.data(:,:,pw); 
y = tgc1.*y; 
y = double(y(:));
u0 = DAS_RF;
v0 = DAS_RF;
lambd0 = DAS_RF;
eps = 5e-2;
mu = 1e5;
beta = 1000;
maxiter = 100;
inner_it = 1;
[u, v, lambd, objective] = RED_v0(Phi, y, u0, v0, lambd0, eps, mu, beta, maxiter, inner_it);
L1_RF = u(:,3)./max(abs(u(:,3)));
L1_ENV = abs(hilbert(reshape(L1_RF,[scan.Nz scan.Nx])));
%% Visualization
dynamic_range = 50;
%-- setting axis limits (mm)
x_lim = [min(scan.x_axis) max(scan.x_axis)]*1e3;
z_lim = [min(scan.z_axis) max(scan.z_axis)]*1e3;
%-- compute dB values
DAS_Bmode = 20*log10(DAS_ENV./max(DAS_ENV(:))+0.001);
B_Bmode = 20*log10(L1_ENV./max(L1_ENV(:))+0.001);
vrange = [-dynamic_range 0];
%-- display image
figure
subplot(121)
imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,DAS_Bmode);
shading flat; colormap gray; caxis(vrange); colorbar; hold on;
axis equal manual;
xlabel('x [mm]');
ylabel('z [mm]');
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
axis([x_lim z_lim]);
title('DAS')
drawnow; hold off;
subplot(122)
imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,B_Bmode);
shading flat; colormap gray; caxis(vrange); colorbar; hold on;
axis equal manual;
xlabel('x [mm]');
ylabel('z [mm]');
set(gca,'YDir','reverse');
set(gca,'fontsize',16);
axis([x_lim z_lim]);
title('RED')
drawnow; hold off;