%% Running the Test with 1 second Windows

clear
clc

addpath Necessary
addpath test

%Set signal to noise ratios
SNR = 20:-5:0;
SNR_len = length(SNR);

cd test
listing = dir('*.mp3');
cd ..

tks = struct2cell(listing)';
tks(:,2:5) = [];

clear_hashtable
add_tracks(tks);
num_tks = length(tks);
solution = {};
num_wins = [];
num_corr_elements = [];

for track = 1:num_tks
    % Load a query waveform (recorded from playback on a laptop)
    % Loads and queries 1 second of every mp3 with varying SNR ratios, and
    % the 1 second slides by a half second across the entire file.
    [query,fs] = audioread(tks{track});
    query_len = length(query);
    num_wins(track) = floor(query_len/fs/2.5) - 1;
    
    for j = 1:num_wins(track)
        query_t = query((j-1)*2.5*fs+1:(j-1)*2.5*fs+1+5*fs,:);
        pwr = rms(sum(query_t,2))^2;
        for k = 1:SNR_len
            noise_pwr = pwr * 10.^(-SNR./10);
            len = length(query_t);
            w_n = SNR(k)*rand(len,1);
            query_t = query_t + repmat(w_n,1,2);
            R = match_query(query_t,fs);
            if size(R) == [0,4]
                solution{track, j, k} = [0 0 0 0];
            else
                solution{track, j, k} = R(1,:);
            end
        end
    end
    % Stores the match_query information into a 3D cell, the first
    % dimension is the track number, the second dimension is the location
    % of the query in the track, and the third dimension is the signal to
    % noise ratio.
    
    disp(['Completed Track ', num2str(track), ' of ', num2str(num_tks),'.']);
end

% illustrate_match(dt,srt,tks);
% colormap(1-gray)

save('ShazamData5sec.mat','solution');

%% Checking for the number of correct elements

for i = 1:track
    num_corr_elements(:,i) = reshape(sum(cellfun(@(x) x(1), solution(i,1:num_wins(i),:)) == i,2), [SNR_len 1 1]);
    num_corr_elements(:,i) = num_corr_elements(:,i)./num_wins(i);
end

plot(SNR, num_corr_elements) %Change to SNR
title('Percentage of Correctly Matched Tracks vs SNR for Five Second Windows')
xlabel('Signal to Noise Ratio (dB)')
ylabel('Correct Match Percentage')
legend(tks{:})