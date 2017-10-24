%% Running the Test with 1 second Windows

clear
clc

addpath Necessary
addpath test

% Set the length of the window in seconds
win_len = 1;

%Set signal to noise ratios
SNR = 20:-5:0;
SNR_len = length(SNR);

cd test
% Find all MP3 files
listing = dir('*.mp3');
cd ..

tks = struct2cell(listing)';
tks(:,2:5) = [];

% Initialized Variables
clear_hashtable
add_tracks(tks);
num_tks = length(tks);
solution = {};
num_wins = [];
num_corr_elements = [];
query_t = {};
query_call = {};

[query, fs] = audioread(tks{1});

for track = 1:num_tks
    % Load a query waveform (recorded from playback on a laptop)
    % Loads and queries 1 second of every mp3 with varying SNR ratios, and
    % the 1 second slides by a half second across the entire file.
    [query,fs] = audioread(tks{track});
    query_len = length(query);
    num_wins = floor(query_len/(fs*win_len)*2)-1;
    temp_index = [];
    
    % Check if the video has stereo or mono audio, if stereo shrink to mono
    if size(query,2) == 2
        stereo_flag = 2;
        query = (query(:,1) + query(:,2))./2;
    else
        stereo_flag = 1;
    end
    
    %Buffer the track out into separated windows
    query_t = buffer(query, fs*win_len, fs*win_len/2);
    pwr = rms(query_t).^2;
    noise_pwr = repmat(pwr, [1,1,SNR_len]) .* repmat(reshape(10.^(-SNR./10),[1,1,SNR_len]),1,length(pwr));
    query_t = repmat(query_t(:,1:num_wins),1,1,SNR_len) + repmat(randn(fs*win_len,1),1,num_wins,SNR_len).*repmat(noise_pwr(:,1:num_wins,:),fs*win_len,1,1);
    
    query_call{track} = query_t;
    
end

disp(['Completed Moving stuff around']);

for track = 1:num_tks
    query = query_call{track};
    for sn = 1:SNR_len
        shift_max = size(query,2);
        for shift = 1:shift_max
            solution{track, shift, sn} = match_query(query(:,shift,sn),fs);
        end
    end
    disp(['Completed Track ', num2str(track), ' of ', num2str(num_tks),'.']);
end

        
    
% Stores the match_query information into a 3D cell, the first
% dimension is the track number, the second dimension is the location
% of the query in the track, and the third dimension is the signal to
% noise ratio.

    

% illustrate_match(dt,srt,tks);
% colormap(1-gray)

save(['Shazam_Data_', num2str(win_len),'_sec.mat'],'solution');

%% Checking for the number of correct elements

for i = 1:num_tks
    num_corr_elements(:,i) = reshape(sum(cellfun(@(x) x(1), solution(i,1:num_wins(i),:)) == i), [SNR_len 1 1]);
    num_corr_elements(:,i) = num_corr_elements(:,i)./num_wins(i);
end

plot(SNR, num_corr_elements) %Change to SNR
title('Percentage of Correctly Matched Tracks vs SNR for One Second Windows')
xlabel('Signal to Noise Ratio (dB)')
ylabel('Correct Match Percentage')
legend(tks{:})