% Alyssa Holman
% Last edited: 10/11/2023
% Calcium Analysis

%% Set paths and open files
% First, set the path below to open the images of interest
clc
clear
cd('/Users/Alyssa/Dropbox/PhD UCSD/Thesis/Ca2+ Imaging/2D_3D/Finalized Ca2+/D25 ex for plot/');

% import linescanning images
imagefiles = dir('*.tif')  % information about each .tif  
nfiles = length(imagefiles);    % number of files found

% get info from files
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   currentimage = imread(currentfilename);
   cropcurrent = currentimage;
   images{ii} = cropcurrent; % extract pixel intensity for each image
end

% set the conversion from frames to seconds
fps = 8; % frames per second
all_frames = 200; % number of frames total

all_pix = length(mean(images{1},2));
conv = (all_frames/fps)/all_pix; % gets conversion based on time / frames


%% Finding peaks and minima for intensity values
% Here, we want to find the maxima and minima so we can calculate the
% features within the ca2+ flux plots

for x=1:length(images)
    y = 1:size(images{x},1);
    m_z = mean(images{x},2);
    z = smooth(m_z, 0.05, 'loess');
    
    inv = -z;
    
    % preliminary peaks to get ranges
    [IntMaxs_tmp, TimeMaxs_tmp] = findpeaks(z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', 2);
    [IntMins_temp, TimeMins_tmp] = findpeaks(-z, y,'MinPeakDistance', 0.5, 'MinPeakProminence', 2);
    IntMins_tmp = -IntMins_temp;
    
    % get actual max and min values based on baseline
    prom_min = 0.2*(mean(IntMaxs_tmp) - mean(transpose(IntMins_tmp)));
    prom_max = 0.2*(mean(IntMaxs_tmp) - mean(transpose(IntMins_tmp)));
    
    [IntMaxs, TimeMaxs, width, prom] = findpeaks(z, y,'MinPeakDistance', 200, 'MinPeakProminence', prom_max);
    [IntMins_neg, TimeMins] = findpeaks(-z, y,'MinPeakDistance', 200, 'MinPeakProminence', prom_min);
    IntMins = -IntMins_neg;
    
    % save values
    IntMax{x} = IntMaxs; % save Intensity Maximum vals (i.e. intensity at peak)
    TimeMax{x} = TimeMaxs; % save Time Maximum vals (i.e. time of peak)
    IntMin{x} = IntMins; % save Intensity Minimum vals (i.e. intensity at minima)
    TimeMin{x} = TimeMins; % save Time Minimum vals (i.e. time of minima)
    Proms{x} = prom;
    width_fin = width*conv;
    Widths{x} = width_fin;
    
    % plot to make sure peaks and minima are correct
    sub = length(images)/12;
    subplot(ceil(sub),12,x); 
    hold on
        plot(y, z)
        title('Resting and Maximum Points');
        plot(TimeMins, IntMins, '*') % adds * where mins are
        plot(TimeMaxs, IntMaxs, '*') % adds * where maxs are
        title(imagefiles(x).name);
        xlabel('Time in Pixels');
        ylabel('Intensity');
    hold off
    
end

%% Finding Delta F and Fo values from intensity values
% First, calculate resting flux and use this to get F/Fo

for x=1:length(images)
    
    % get resting values
    restings = mean(IntMin{x}(2:end,1)); % skipped first val because sometimes not accurate
    resting{x} = restings; % save resting values
    
    % get deltaF values
    m_z = mean(images{x},2);
    z = smooth(m_z, 0.05, 'loess');
    deltF = [];
    
    for i = 1:size(z,1)
        v = (((z(i,1) - restings)/restings)+0.06); % take (intensity - resting) / resting
        deltF = [deltF, v]; % these are the delta(F/Fo) values
    end
    
    % save deltaF values
    deltaF{x} = deltF;
    
    % plot with deltaF values
    sub = length(images)/2;
    subplot(ceil(sub),2,x); 
    y = (1:size(images{x},1));
    plot(conv*y,deltF, 'LineWidth',2.0, 'color', 'k');
    set(gca,'YTick',[])
    title(imagefiles(x).name);
    xlabel('Time in Seconds');
    ylabel('(F-Fo)/Fo');  

    clear restings
    clear deltF
end


%% Get delta values for first peak

for x=1:length(images)

    % import values
    deltF = deltaF{x};
    TimeMins = TimeMin{x};
    TimeMaxs = TimeMax{x};
    IntMins = IntMin{x};
    IntMaxs = IntMax{x};
    deltF = deltaF{x};
    restings = resting{x};

    % skip file if a min/max value is skipped from findpeaks()
    if length(TimeMaxs) > length(TimeMins)+1 == 1 | length(TimeMaxs)+1 < length(TimeMins) == 1;
            continue;
    end

    % skip file if there aren't enough peaks
    num_peak=1:length(TimeMins)-1;
    if isempty(num_peak) == 1;
        continue;
    end
    
    % assuming enough peaks and correct determination of mins/maxs:
    % loop to get first peak or "n" peak
    for num_peak=1:length(TimeMins)-1;
  
    % first "if" statement is if 'findpeaks()' or 'islocalmin()' improperly found
    % 2 datapoints back to back
        if TimeMins(1,num_peak+1) < TimeMaxs(1,num_peak) 
            valo = TimeMins(1,num_peak+1); % first time minima (i.e. before peak)
            valt = TimeMins(1,num_peak+2); % second time minima (i.e. after peak)
            val_max = TimeMaxs(1,num_peak); % maxima time (i.e. peak)
            deltFo = deltF(1,TimeMins(1,num_peak+1)); % F/F0 before peak
            deltFmax = deltF(1,TimeMaxs(1,num_peak)); % F/F0 at peak

        % second "if" statement is if the first peak is cut off, so the second
        % peak will be selected
        elseif (IntMins(num_peak,1) < (restings - 75)) | (IntMins(num_peak,1) > (restings + 75))
            valo = TimeMins(1,num_peak);
            valt = TimeMins(1,num_peak+1);
            val_max = TimeMaxs(1,num_peak);
            deltFo = deltF(1,TimeMins(1,num_peak)); 
            deltFmax = deltF(1,TimeMaxs(1,num_peak));

        % otherwise, peak should be normal
        else
            valo = TimeMins(1,num_peak);
            valt = TimeMins(1,num_peak+1);
            val_max = TimeMaxs(1,num_peak);
            deltFo = deltF(1,TimeMins(1,num_peak)); 
            deltFmax = deltF(1,TimeMaxs(1,num_peak));
        end 
        
        % save values
        valone{num_peak} = valo;
        valtwo{num_peak} = valt;
        val__max{num_peak} = val_max;
        raw_tau_{num_peak} = val_max - valo; % time to peak
        delta_F1{num_peak} = deltFo;
        delta_Fmax{num_peak} = deltFmax;
        delta_F_diff{num_peak} = deltFmax - deltFo;
        
    end
    
    val1{x} = valone;
    val2{x} = valtwo;
    valmax{x} = val__max;
    raw_tau{x} = raw_tau_; % time to peak
    deltaF1{x} = delta_F1;
    deltaFmax{x} = delta_Fmax;
    deltaF_diff{x} = delta_F_diff;
   
    % plot with deltaF values
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    y = valone{1}:valtwo{1};
    plot(y*conv,deltF(valone{1}:valtwo{1})) % conversion 'conv' set before this loop
    title(imagefiles(x).name);
    xlabel('Time in Seconds');
    ylabel('(F-Fo)/Fo');
    ylim([-0.4 (max(deltF)*1.2)]);
    
    clear valo
    clear valt
    clear val_max
    clear deltFo
    clear deltFmax
    clear valone
    clear valtwo
    clear val__max
    clear raw_tau_
    clear delta_F1
    clear delta_Fmax
    clear delta_F_diff
    clear deltF
    clear TimeMins
    clear TimeMaxs
    clear IntMins
    clear IntMaxs
    clear restings
end


%% Calculate tau value
% tau is time at which F/Fo is 63% to peak
clear length

percent = 0.63;

for x=1:length(images)

    % import in
    deltF_diff = deltaF_diff{x};
    deltFo = deltaF1{x};
    deltFmax = deltaFmax{x};
    valo = val1{x};
    valt = val2{x};
    val_max = valmax{x};
    deltF = deltaF{x};

    % skip file if there aren't enough peaks
    num_peak=1:length(deltF_diff)-1;
    if isempty(num_peak) == 1;
        continue;
    end
    
    % get all tau
    for num_peak=1:length(deltF_diff)-1;
        
        if isempty(num_peak) | length(deltF_diff) < 2;
            f_tau = NaN;
            tau = NaN;
        	norm_tau = NaN;
        	lambda = NaN;
        	norm_lambda = NaN;
            print("empty")
            continue
        
        else    
            % F tau, flux intensity when tau occurs
            f_tau = (deltF_diff{num_peak}*percent) + deltFo{num_peak}; % F tau val

            % set range for tau
            range = 0.01;

            % time tau
            if valo{1} > val_max{1}
                for j=valo{num_peak}:val_max{num_peak+1}
                    if deltF(j) < f_tau + range && deltF(j) > f_tau - range % finding approximation
                        tau = j; % not normalized to time = 0
                        norm_tau = j - valo{num_peak}; % normalized to start of first peak at time = 0
                        lambda = 1/tau;
                        norm_lambda = 1/norm_tau;
                        break
                    end
                end
                for j=valo{num_peak}:val_max{num_peak+1}
                    if exist('tau', 'var') == 0
                        tau = (val_max{num_peak}-valo{num_peak})*percent + valo{num_peak};
                        norm_tau = tau - valo{num_peak}; % normalized to start of first peak at time = 0
                        lambda = 1/tau;
                        norm_lambda = 1/norm_tau;
                        break
                    end
                end
            elseif valo{1} < val_max{1}
                for j=valo{num_peak}:val_max{num_peak}
                    if deltF(j) < f_tau + range && deltF(j) > f_tau - range % finding approximation
                        tau = j; % not normalized to time = 0
                        norm_tau = j - valo{num_peak}; % normalized to start of first peak at time = 0
                        lambda = 1/tau;
                        norm_lambda = 1/norm_tau;
                    end
                end
                for j=valo{num_peak}:val_max{num_peak+1}
                    if exist('tau', 'var') == 0
                        tau = (val_max{num_peak}-valo{num_peak})*percent + valo{num_peak};
                        norm_tau = tau - valo{num_peak}; % normalized to start of first peak at time = 0
                        lambda = 1/tau;
                        norm_lambda = 1/norm_tau;
                    end
                end
            end
            
            f_taus{num_peak} = f_tau;
            taus{num_peak} = tau;
            taus_ms{num_peak} = tau*conv; % convert to ms
            taus_ms_norm{num_peak} = norm_tau*conv;
            lambdas{num_peak} = lambda;
            lambdas_ms{num_peak} = lambda*conv;
            lambdas_ms_norm{num_peak} = norm_lambda*conv;
        end
        
    end
    
    un = unique(cell2mat(taus_ms_norm));
    mod_taus = isoutlier(un, 'mean',  'ThresholdFactor', 1); % if calc incorrectly, remove
    
    % save values
    allftau{x} = f_taus;
    alltau{x} = taus;
    alltau_ms{x} = taus_ms; % convert to ms
    alltau_ms_norm{x} = un(~mod_taus);
    alllambda{x} = lambdas;
    alllambda_ms{x} = lambdas_ms;
    alllambda_ms_norm{x} = lambdas_ms_norm;

    clear taus_ms_norm
    clear un
    clear mod_taus
    clear taus_ms

    
    % val you choose to plot
    choice = 1;
    
    % plot with tau values to make sure they are correct
    sub = length(images)/4; % to make a grid of plots in 4 columns
    subplot(ceil(sub),4,x); % to specify number of rows
    
    hold on
    y = valo{choice}:valt{choice};
    plot(y,deltF(valo{choice}:valt{choice}))
    plot(alltau{x}{choice}, allftau{x}{choice}, '*') % add * where tau is
    
    title(imagefiles(x).name);
    xlabel('Time in Milliseconds');
    ylabel('(F-Fo)/Fo');
    ylim([-0.4 (max(deltF)*1.2)]);
    hold off
    
end

%% Find and plot time to peak using all peaks in data

% get peaks
for x=1:length(images)

    TimeMins = TimeMin{x};
    TimeMaxs = TimeMax{x};

    n = 1:(length(TimeMaxs)-2);
    if isempty(n) == 1;
        continue;
    end
    
    % calculating time to peak and cycle length
    for n = 1:(length(TimeMaxs)-2)
        if TimeMins(1,1) > TimeMaxs(1,1);
            TimetoPeakss = TimeMaxs(1,n+1) - TimeMins(1,n);
            CycleLength(n) = (TimeMins(1,n+1) - TimeMins(1,n))*conv;
            TimetoPeaks(n) = TimetoPeakss*conv;
        elseif TimeMins(1,1) < TimeMaxs(1,1)
            TimetoPeakss = TimeMaxs(1,n) - TimeMins(1,n);
            CycleLength(n) = (TimeMins(1,n+1) - TimeMins(1,n))*conv;
            TimetoPeaks(n) = TimetoPeakss*conv;
        end
    
        % remove outliers
        Mod_TimetoPeaks = isoutlier(TimetoPeaks, 'mean',  'ThresholdFactor', 1.75);
        TimetoPeak{x} = TimetoPeaks(~Mod_TimetoPeaks);
        
        % remove outliers
        new_cycle = isoutlier(CycleLength, 'mean', 'ThresholdFactor', 1);
        
        Cycle_Length{x} = CycleLength(~new_cycle);
        val_std = std(TimetoPeaks(~Mod_TimetoPeaks));
        val_std_peak{x} = val_std;
        val_error = val_std/length(TimetoPeaks(~Mod_TimetoPeaks));
        val_error_peak{x} = val_error;
        
    end
    
    AvgTimetoPeak = mean(TimetoPeaks(~Mod_TimetoPeaks));
    AvgTimetoPeaks{x} = AvgTimetoPeak;
    AllTimetoPeaks{x} = TimetoPeaks(~Mod_TimetoPeaks);
    AvgCycleLength = mean(CycleLength(~new_cycle));
    AvgCycleLengths{x} = AvgCycleLength;
    AllCycleLength{x} = CycleLength(~new_cycle);
    
    clear TimeMins
    clear TimeMaxs
    clear new_cycle
    clear Mod_TimetoPeaks
    clear TimetoPeakss
    clear CycleLength
    clear TimetoPeaks
end

%% make table with variables and export
% rename files so Matlab accepts the format (letter then number)
for x=1:length(images)
    filename = imagefiles(x).name;
    filename_split = strsplit(filename,'.');
    filenames{x} = strcat(filename_split{1},filename_split{2});
end

%% Gets all of the values into an excel sheet
A = nan(length(images)*6, 30); % if more values, increase from 30
variable_labels = [];
for i=1:length(alltau_ms_norm)
    horiz = horzcat(alltau_ms_norm{i});
    A(i, 1:length(horiz)) = horiz;
    name = 'tau(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(AllTimetoPeaks)
    horiz = horzcat(AllTimetoPeaks{i});
    new = i+length(alltau_ms_norm);
    A(new, 1:length(horiz)) = horiz;
    name = 'timetopeak(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(deltaF_diff)
    horiz = horzcat(cell2mat(deltaF_diff{i}));
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks);
    A(new, 1:length(horiz)) = horiz;
    name = 'deltaF/F0';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Cycle_Length)
    horiz = horzcat(Cycle_Length{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff);
    A(new, 1:length(horiz)) = horiz;
    name = 'cyclelength(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Widths)
    horiz = horzcat(Widths{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff)+length(Cycle_Length);
    A(new, 1:length(horiz)) = horiz;
    name = 'widths(ms)';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end
for i=1:length(Proms)
    horiz = horzcat(Proms{i});
    new = i+length(alltau_ms_norm)+length(AllTimetoPeaks)+length(deltaF_diff)+length(Cycle_Length)+length(Widths);
    A(new, 1:length(horiz)) = horiz;
    name = 'prominance';
    variable_labels = [variable_labels, sprintf("%1$s_%2$s", filenames{i}, name)];
    clear horiz
end

final_tab = array2table(A);
new_table = [ table(transpose(variable_labels), 'VariableNames',{'Labels'})  final_tab];
writetable(new_table, 'ca2+_data.csv'); % writes to folder

