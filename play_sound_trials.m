dir_code = '/Users/pinheirochagas/Pedro/Stanford/code/brain_music';
data_fn = [dir_code filesep 'data.mat'];
load(data_fn)

% Select
electrodes = [105, 15, 61, 81, 70];
elect_regions = {'lateral occipital gyrus (LOG)', 'posterior inferior temporal gyrus (pITG)','anterior intraparietal sulcus (aIPS)','superior parietal lobe (SPL)','premotor cortex (PM)'};
trials = data.trialinfo_all{1};
trials_no = find(strcmp(trials.condNames, 'math') & trials.spike_hfb==false & trials.RT < 5 & trials.RT > 1);
trials_final = trials(trials_no,:);
[~,idx] = sort(trials_final.RT);
trials_final = trials_final(idx,:);

data_plot = squeeze(data.wave(trials_no,3,:));
data_plot = data_plot(idx,:);

% Preapare sound
amp = 10;
fs = 20500;  % sampling frequency
duration = 0.05;
values = 0:1/fs:duration;

% sound_all = [];
% Transform data to sound

for ielec = 1:length(electrodes)
    data_plot = squeeze(data.wave(trials_no,ielec,:));
    data_plot = data_plot(idx,:);

    b_all = [];
    for i = 1:10%:size(data_plot,1) % take first 10 trials
        RT_idx = (round(trials_final.RT(i)*data.fsample) + abs(min(data.time)*data.fsample));
        times = 1:RT_idx;
        data_plot_tmp = data_plot(i,:);
        data_plot_tmp(times(end):(abs(data.time(1))+data.time(end))*data.fsample+1) = 0;
        
        data_plot_tmp(data_plot_tmp<=0) = 0;
        peak_vals_trans = round(data_plot_tmp*100);
        peak_vals_trans(peak_vals_trans==0) = 1;
        peak_vals_trans = downsample(peak_vals_trans,10);
        
        peak_vals_trans(peak_vals_trans<100) = 99999;  
        
        b = [];
        for iss = 1:length(peak_vals_trans)
            peak_vals_trans(iss);
            freq = peak_vals_trans(iss);
            a = amp*sin(2*pi* freq*values);
            if freq == 99999
                a(:) = 0;
            else
            end
            a_all(iss,:) = a;
            b = [b a]; %zeros(100,1)'
        end
        b_all = [b_all;b];
        
%         audiowrite(sprintf('%s/mmr_audio/mmr_audio_%s.wav',dir_code, num2str(i)),b,10000);
    end
    b_mean = mean(b_all);  
    b_sum = sum(b_all);  
    
    %% Plot spectrogram of all trials
%     subplot(length(electrodes),1,ielec)
%     pspectrum(b_sum,10000,'spectrogram','TimeResolution',0.1,'OverlapPercent',0,'MinThreshold',-20);
%     ylim([0 2])
%     xlim([0 11])
%     xlabel('Time')
%     set(gca,'xtick',[])
%     set(gca,'fontsize',16)
%     title(elect_regions{ielec}) 
%     
%     colormap(viridis)
    b_sum = bandpass(b_sum,[300 2000], 10000);
    audiowrite(sprintf('%s/mmr_audio/mmr_audio_elect_%s.wav', dir_code, num2str(electrodes(ielec))),b_sum,10000);
end
fname = [dir_code filesep 'spectrogram_elecs.png'];
savePNG(gcf,600,fname)


figure('units', 'normalized', 'outerposition', [0 0 0.4 0.5])
subplot(3,1,1)
plot(data.time, data_plot(i,:), 'LineWidth',2, 'Color', 'k')
xlim([-0.2 5])
ylabel('zscore HFB')
xlabel('Time (sec.)')
title('HFB of one epoch')
subplot(3,1,2)
plot(peak_vals_trans, 'LineWidth',2, 'Color', 'k')
ylabel('zscore HFB')
xlabel('Time (sec.)')
title('Transformed and downsampled signal of one epoch only until RT')
subplot(3,1,3)
plot(b(1,1:end/3), 'LineWidth',1, 'Color', 'k')
title('Example of sine waves of one given epoch')
fname = [dir_code filesep 'example_signal.png'];
savePNG(gcf,300,fname)



for ie = 1:length(electrodes)
    
    figure('units', 'normalized', 'outerposition', [0 0 0.2 0.6])
    data_plot = squeeze(data.wave(trials_no,ie,:));
    data_plot = data_plot(idx,:);
    
    n_plot = size(data_plot,1);
    cols = viridis(size(data_plot,1));
    for i = 1:size(data_plot,1)
        hold on
        RT_idx = (round(trials_final.RT(i)*data.fsample) + abs(min(data.time)*data.fsample));
        times = 1:RT_idx;
        data_plot_tmp = data_plot(i,:);
        data_plot_tmp(times(end)+data.fsample/5:(abs(data.time(1))+data.time(end))*data.fsample) = nan;
        [~,peaks] = findpeaks(data_plot_tmp+10*i);
        [peak_vals,~] = findpeaks(data_plot_tmp);
        peak_vals_avg(i,ie) = mean(peak_vals(peak_vals>0));
        %         plot(data.time(peaks), data_plot_tmp(peaks)+10*i, '.', 'Color', cols(i,:), 'MarkerSize', 6)
        plot(data.time, data_plot_tmp+10*i, 'k', 'LineWidth', 1.5, 'Color', cols(i,:))
        plot(data.time(RT_idx), data_plot_tmp(RT_idx)+10*i, '.', 'MarkerSize', 5, 'Color', 'k')
        
    end
    xlim([-0.250 5])
    xline(0, 'LineWidth', 2, 'Color', 'k')
    set(gcf,'color', 'w')
    set(gca,'color','w')
    axis on
    box on
    ylabel('Trials sorted by RT')
    xlabel('RT (s)')
    ylim([0 1500])
    set(gca,'fontsize',16)
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    fname = sprintf('%s/erp_plot_%s.png',dir_code, num2str(ie));
    savePNG(gcf,300,fname)
    close all
end









% Preapare sound

% Sine
amp = 10;
fs = 20500;  % sampling frequency
duration = 0.05;
values = 0:1/fs:duration;
freq = 400;
a = amp*sin(2*pi* freq*values);

%Gaussian
x = [-3:.1:3];
y = normpdf(x,0,1);
plot(x,y)


c = a*y


%% Directly to MIDI

notes = {'d', 'd#','r', 'r#', 'm', 'f', 'f#', 'so', 'so#', 'l', 'l#', 'si'}';
notes_all = cell(length(notes),1,1);
count = 1;
for io = 1:10
    start_oc = num2str(io);
    for i = 1:length(notes)
        notes_all{count} = [notes{i} start_oc];
        count = count+1;
    end
end
notes_all(1) =[]
notes_all_full = cell(127,1,1);
notes_all_full(1:size(notes_all,1)) = notes_all
notes_all_full(end-7:end) = {'d11', 'd#11','r11', 'r#11', 'm11', 'f11', 'f#11', 'so11',};

midi_mapping = table;
midi_mapping.values = [1:127]';
midi_mapping.notes = notes_all_full;

notes_elec = {'d6', 'm6', 'so6', 'si6', 'd7'};


data_elec = [];
for ielec = 1:length(electrodes)
    data_plot = squeeze(data.wave(trials_no,ielec,:));
    data_plot = data_plot(idx,:);

    data_plot_filt = [];
    for i = 1:size(data_plot,1) % take first 10 trials
        RT_idx = (round(trials_final.RT(i)*data.fsample) + abs(min(data.time)*data.fsample));
        times = 1:RT_idx;
        data_plot_tmp = data_plot(i,:);
        data_plot_tmp(times(end):(abs(data.time(1))+data.time(end))*data.fsample+1) = nan;
        data_plot_tmp(data_plot_tmp<=0) = 0;
        data_plot_tmp = downsample(data_plot_tmp,10);
        data_plot_filt(i,:) = data_plot_tmp;
    end
    data_elec_sum  = sum(data_plot_filt);
    data_elec_sum = round(data_elec_sum/max(data_elec_sum) * 127); % find better way to scale;
    data_elec_sum(isnan(data_elec_sum)) = 1;
    
    
    data_elec(:,ielec) = data_elec_sum;
    vols = round(data_elec_sum/max(data_elec_sum) * 100)';
%     
%     single_note = data_elec_sum;
%     single_note(single_note>1) = 100;
    % Make MIDI
    
    note = midi_mapping{find(strcmp(midi_mapping.notes, notes_elec{ielec})),1};
    
    
    % initialize matrix:
    N = size(data_elec,1);
    M = zeros(N,6);
    M(:,1) = 1;         % all in track 1
    M(:,2) = 1;         % all in channel 1
    M(:,3) = repmat(note,size(data_elec,1), 1);  % random note numbers data_elec(:,ielec)
    M(:,4) = vols;  % random volumes
    M(:,5) = linspace(1,size(data_elec,1), size(data_elec,1));
    M(:,6) = M(:,5) + .1;  % random duration .2 -> 1.2 seconds
    midi_new = matrix2midi(M);
    writemidi(midi_new, sprintf('%s/midi_%s.mid', dir_code, num2str(ielec)));
end





% initialize matrix:
N = 200;
M = zeros(N,6);

M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,3) = 30 + round(60*rand(N,1));  % random note numbers
M(:,4) = 60 + round(40*rand(N,1));  % random volumes
M(:,5) = 10 * rand(N,1);
M(:,6) = M(:,5) + .2 + rand(N,1);  % random duration .2 -> 1.2 seconds

midi_new = matrix2midi(M);
writemidi(midi_new, [dir_code '/testout5.mid']);


% initialize matrix:
N = 261;
M = zeros(N,6);

M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,3) = peak_vals_trans';  % random note numbers
M(:,4) = 60 + round(40*rand(N,1));  % random volumes
M(:,5) = 10 * rand(N,1);
M(:,6) = M(:,5) + .2;  % random duration .2 -> 1.2 seconds

midi_new = matrix2midi(M);
fname = sprintf('%s/erp_plot_%s.mid',dir_code, num2str(ie));
writemidi(midi_new, fname);


% initialize matrix:
N = 200;
M = zeros(N,6);

M(:,1) = 1;         % all in track 1
M(:,2) = 1;         % all in channel 1
M(:,3) = 30 + round(60*rand(N,1));  % random note numbers
M(:,4) = 60 + round(40*rand(N,1));  % random volumes
M(:,5) = 1.01:0.1:21;
M(:,6) = M(:,5)+0.1;  % random duration .2 -> 1.2 seconds

midi_new = matrix2midi(M);
writemidi(midi_new, [dir_code '/testout3.mid']);


M(:,3) = peak_vals_trans(1:200)

