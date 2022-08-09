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
path_dataset = ['../PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_dataset_rf','.hdf5'];% update the folder directory
path_scan = ['../PICMUS/database/',acquisition,'/',phantom,'/',phantom,'_',suffix,'_scan','.hdf5'];% update the folder directory
%% loading dataset
%-- load scan and dataset
scan = linear_scan();
scan.read_file(path_scan);
dataset = us_dataset();
dataset.read_file(path_dataset);
%-- define scan based on time axis
time = (150:1400).'/dataset.sampling_frequency+dataset.initial_time;
z_axis= time*dataset.c0/2;
scan = linear_scan(scan.x_axis,z_axis);
%% Image reconstruction
%-- receive apodization
%-- dynamically expanding receive aperture with hanning apodization
rx_f_number = 0.5;
rx_aperture = scan.z/rx_f_number;
rx_aperture_distance = abs(scan.x*ones(1,dataset.channels)-ones(scan.pixels,1)*dataset.probe_geometry(:,1).');
receive_apodization = tools.apodization(rx_aperture_distance,rx_aperture*ones(1,dataset.channels),'hanning');
%-- angular apodization
angular_apodization = 1;
%-- matrix Phi construction
ins = 38; %select the desired insonification angle
time_vector = dataset.initial_time+(0:(dataset.samples-1))/dataset.sampling_frequency;
tau_m = time_vector;
tau_0 = dataset.initial_time;
tau_t = scan.z*cos(dataset.angles(ins))+scan.x*sin(dataset.angles(ins));
b_k = dataset.data(:,:,ins); b_k = b_k(:);
rr = length(b_k);
cc = length(scan.x);
count=0;
wb = waitbar(0,'matrix construction');
for nrx=1:dataset.channels
    %-- progress bar
    step = nrx/dataset.channels;
    waitbar(step,wb,sprintf('matrix construction %0.0f%%',step*100));
    %
    tau_r = sqrt((dataset.probe_geometry(nrx,1)-scan.x).^2+(dataset.probe_geometry(nrx,3)-scan.z).^2);
    delay = (tau_t+tau_r)/dataset.c0;
    A_k = sparse(rr,cc);
    for nrz=1:length(tau_m)
        count=count+1;
        diff1 = abs(delay-tau_m(nrz));
        diff2 = abs(delay-tau_m(nrz))<=(1/dataset.sampling_frequency);
        if (~isempty(nonzeros(diff2)))
            diff1 = diff1.*diff2;
            diff1 = diff1/max(diff1); diff1 = 1-diff1;
            aa = (diff1.*diff2).*receive_apodization(:,nrx);
            [i1,j1,s1] = find(transpose(aa));
            A_k = A_k + sparse(count*i1,j1,double(s1),rr,cc);
        end
    end
    save(['A_k_',num2str(nrx)],'A_k')
end
close(wb);