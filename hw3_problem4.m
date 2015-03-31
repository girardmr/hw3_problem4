clear all;
clc;

load('anecoicCystData.mat');

data = veraStrct.data(80:end,:,:);
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

for beam = 1:128
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5 + depth;
    time_to_point(ii) = d(ii)/speed;
end

delay_matrix(jj,:,beam) = time_to_point; %delays

end

for aa = 1:128
    delayed_channel(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
end


end


for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end
summed_channels = sum(delayed_channel,2);

%part a

beam_data = summed_channels(:,64);
beam_data_FT = abs(fftshift(fft(beam_data)));
beam_data_FT_dB = 20*log10(beam_data_FT);
figure;
plot(beam_data_FT_dB);
title('Magnitude of Fourier transform of beam 64');
xlabel('Frequency (radians)');
ylabel('dB');

%part b

for rf_beam = 1:128
    beam_data_rf = summed_channels(:,rf_beam);
    beam_data_FT_rf = abs(fftshift(fft(beam_data_rf)));
    beam_data_FT_dB_rf(:,rf_beam) = 20*log10(beam_data_FT_rf);
    averaged_FT_beam_rf(:,rf_beam) = mean(beam_data_FT_dB_rf(:,rf_beam));
end
figure;
plot(averaged_FT_beam_rf);
title('Average magnitude of FT of each RF line');
xlabel('Beam');
ylabel('Average magnitude');

%part c
normalized = beam_data_FT_dB-max(beam_data_FT_dB);
figure;
plot(normalized);
title('Normalized FT, beam 64');
%-6 dB frequencies ~600, 1800;
w1 = 570/2353;
w2 = 1790/2353;
coeff = fir1(40,[w1 w2]);
filtered = filtfilt(coeff,1,beam_data);
filtered_FT = abs(fftshift(fft(filtered)));
filtered_FT_dB = 20*log10(filtered_FT);
figure;
plot(filtered_FT_dB);
title('Magnitude of Fourier transform of beam 64, filtered');

%part d

for rf_line = 1:128
    filtered_beam(:,rf_line) = filtfilt(coeff,1,summed_channels(:,rf_line));
end
figure;
imagesc(20*log10(abs(hilbert(filtered_beam(:,:)))));
colormap('gray');
title('Anecoic cyst data filtered');

