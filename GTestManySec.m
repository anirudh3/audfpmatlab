%% Running the Test with 1 second Windows

clear
clc

addpath Necessary
addpath test

% Set the length of the window in seconds
win_len = 10;

% Set anchor number threshold
anch = 4;

%Set signal to noise ratios
SNR = 20:-5:-10;
SNR_len = length(SNR);

cd test
% Find all MP3 files
listing = dir('*.mp3');
cd ..

% set listing to the first 50 things
% listing = listing(1:50);

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

parfor track = 1:num_tks
    % Load a query waveform (recorded from playback on a laptop)
    % Loads and queries 1 second of every mp3 with varying SNR ratios, and
    % the 1 second slides by a half second across the entire file.
    [query,fs] = audioread(tks{track});
    query_len = length(query);
    num_wins = floor(query_len/(fs*win_len)*2)-1;
    shift_max(track) = num_wins;
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

query = [];

disp(['Completed moving stuff around']);

for track_2 = 1:num_tks
    query = query_call{track_2};
    for sn = 1:SNR_len
        for shift = 1:shift_max(track_2)
            Rin = match_query(query(:,shift,sn),fs);
            if size(Rin,1) == 0
                Rin = [0 0 0 0];
            end
            solution{track_2, shift, sn} = Rin;
        end
    end
    disp(['Completed Track ', num2str(track_2), ' of ', num2str(num_tks),'.']);
end

        
    
% Stores the match_query information into a 3D cell, the first
% dimension is the track number, the second dimension is the location
% of the query in the track, and the third dimension is the signal to
% noise ratio.

    

% illustrate_match(dt,srt,tks);
% colormap(1-gray)

save(['Shazam_Data_', num2str(win_len),'_sec.mat'],'solution');

num_total_elements = [];


%% Checking for the number of correct elements

for i = 1:num_tks
    ind = {};
    track_solution = reshape(solution(i,:,:),[max(shift_max),SNR_len]);
    for j = 1:SNR_len
        solmask(1:shift_max(i),j,i) = logical(cell2mat(cellfun(@(x) sum(sum(x(:,2)>anch)), track_solution(1:shift_max(i),j), 'UniformOutput', false)));
        ind(i,j) = {find(solmask(:,j,i) == 1)};
        num_corr_elements(i,j) = length({track_solution{ind{i,j}}});
    end
    num_total_elements(i,:) = sum(solmask(:,:,i));
end

%% Determining the Confusion Matrix
CM = zeros(num_tks,num_tks,SNR_len);

for   i = 1:num_tks
    ind = [];
    for j = 1:SNR_len
        ind = find(solmask(:,j,i)== 1);
        len = length(ind);
        for k = 1:len
            vals = solution{i,ind(k),j}(solution{i, ind(k), j}(:,2) > anch,1:2);
            len_2 = size(vals,1);
            for l = 1:len_2
                CM(i,vals(l,1),j) = CM(i,vals(l,1),j) + vals(l,2);
            end
        end
    end
end


plot(SNR, num_corr_elements./num_total_elements) %Change to SNR
title(['Percentage of Correctly Matched Tracks vs SNR for ', num2str(win_len), ' Second Windows'])
xlabel('Signal to Noise Ratio (dB)')
ylabel('Correct Match Percentage')
legend(tks{:})