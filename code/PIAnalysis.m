% To read the data

% LTM Sessions will be used as 50 and 51
% Total number of test sessions (prior to LTM tests) for Group 1 is 29 and
% for Group 2 is 31

 clear all
 Timing.Subject(1).Session(1).FileName=[];

% tmp = matlab.desktop.editor.getActive;
% cd(fileparts(tmp.Filename));

filePath = mfilename('fullpath');
cd(fileparts(filePath));

cd('..');

addpath('./code');

a(1).b='./data/peak_data';
%%
geno = readtable('data/genotypes.csv');
%%

ses=30;
cd(a(1).b);
files=[];
files=dir;
sub_count=0;
subjects=[];

for f=4:length(files)   % may need to change initial value
    inVect=[];
    clear Data
    filename=files(f).name;
    inVect=textread(filename,'%s','delimiter','\n');
    
    Fullsub=str2num(filename(1:5));
    if  ismember(Fullsub, subjects)
        sub=sub_count;
    else
        subjects(end+1)=Fullsub;
        sub_count=sub_count+1;
        sub=sub_count;
    end
    
    %ses=str2num(filename(9:10))-24; % this needs change accordingly if you are using only last couple of sessions
    sesStart = strlength(filename) - 5;
    sesEnd = strlength(filename) - 4;
    ses=str2num(filename(sesStart:sesEnd));
    %    if isempty(StartRow)==1; UseStartRow=4; else UseStartRow=StartRow; end %changed to 5
    
    idg = ismember(geno.mouse, Fullsub);
    Geno = geno{idg,2};
    Old = geno{idg,6};
    
    inVect(1) = [] ; %inVect = inVect(2:length(inVect));
    
    for i=1:length(inVect); % If it gives error, the initial value might need to change.
        Data(i,:)=strsplit(inVect{i,:},',');
    end
    
    Timing.Subject(sub).Session(ses).FileName=filename;
    Timing.Subject(sub).Session(ses).Data=Data;
    Timing.Subject(sub).SubjectID=Fullsub;
    Timing.Subject(sub).Geno=Geno;
    Timing.Subject(sub).Old=Old;
    
    
    
end
%%

a(1).b='..';
cd(a(1).b);
%% 
% FirstBatch Subjects [28 29 30]; SecondBatch Subjects [35 36 37 45 49 52 53]; 
% third batch [39 40 41 42 43 44]; fourth batch [49 55 57 59]

subVector = vertcat(Timing.Subject.SubjectID);
subList = [];
for i = 1:length(subVector);
    subList(i) = mod(subVector(i), 1000); %extractBetween((subVector(i)), 2, 4);
end
%%

for sub = 1:sub_count; %[28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59] %[28 29 30] %[28 29 30 35 36 37 45 49 52 53] %[35 36 37 45 49 52 53]%[28 29 30]%% [28 29 30]%[35 36 37 45 49] %[21 22 23 28 29 30]%[1:7 9:18]
    
    % gets max session number for subj
    tmpSes = length(Timing.Subject(sub).Session);
    
    for ses = 1:tmpSes; %51:52
        
        PIIndex=[];
        FIIndex=[];
        PIEnd = [];
        TrayNPIn=[];
        TrayNPOut=[];
        startrow = [];
        UseStart = [];
        
        Data = Timing.Subject(sub).Session(ses).Data;
        
        %'BLOCK' is the start of the experiment
        for i=1:length(Data)
            if contains([Data{i,2}], 'BLOCK') && contains([Data{i,3}], 'START')==1 
                startrow(end+1) = i; UseStart = startrow(1,1) ;
            end
        end
        
        FeedingNum=0;
        for i=1:length(Data)
            
            if contains([Data{i,3}], 'FeederTimerEvent')==1;
                FeedingNum=FeedingNum+1;
            end
        end
        
        for i=UseStart:length(Data); %this will find the lines of PIStart and PIEnd
            % 'PIEvent' is the start of PI trial
            if contains([Data{i,3}], 'PIEvent')==1 % This is when PI trial starts
                PIIndex(end+1) = i;
            end
            % 'PITimerEvent' is the end of PI trial
            if contains([Data{i,3}], 'PITimerEvent')==1 % This is when PI trial end
                PIEnd(end+1) = i;
            end
            
            if contains([Data{i,3}], 'FIEvent')==1 % This is when PI trial end
                FIIndex(end+1) = i;
            end
            
            %now find out the NP in out responses in the data, realize that
            %this are not restricted to the trial limits, we will do it
            %later
            if isempty(strfind([Data{i,3}], 'Hole5NP') & strfind([Data{i,4}], '-1'))==0 % This is NP in
                TrayNPIn(end+1) = i;
            end
            
            if isempty(strfind([Data{i,3}], 'Hole5NP') & strfind([Data{i,4}], '0'))==0 % This is NP out
                TrayNPOut(end+1) = i;
            end      
        end
        
        %quick fix if experiment ends before the PI trial is over
        if length(PIEnd)<length(PIIndex)
            PIIndex(end)= [];
        end        
         
        
        TrialDelimits=[];
        TrialDelimits(1:length(PIIndex),1)=PIIndex;
        TrialDelimits(1:length(PIIndex),2)=1;
         
        Timing.Subject(sub).Session(ses).TrialDelimits=TrialDelimits;
        Timing.Subject(sub).Session(ses).PIIndex=PIIndex;
        Timing.Subject(sub).Session(ses).FIIndex=FIIndex;
        Timing.Subject(sub).Session(ses).PIEnd=PIEnd;
        Timing.Subject(sub).Session(ses).TrayNPIn=TrayNPIn;
        Timing.Subject(sub).Session(ses).TrayNPOut=TrayNPOut;
        Timing.Subject(sub).Session(ses).NumberofReinforcements=FeedingNum;
        
%         Timing.Subject(sub).Session(ses).Data=[];
    end
end

%TrialDelimits(end+1,:)=[length(Data) 1]; % something like this should be
%done so that the last data in the last trialwould be included-I fixed it
%in the following (2nd) section

%% 
% 
%% Range of sessions

%START FROM HERE:

% clear all
% load PeakASD Timing
%%
for sub= 1:sub_count %[28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59]% [28 29 30]%  %[35 36 37 45 49]% [21 22 23 28 29 30]%[1:7 9:18]
    % We are actually excluding the first three subjects from the first
    % batch because of the tube priming issue and the resultant low number
    % of responses in the first animals. Tube was not fully primed in some
    % sessions until the middle of the session for the third animal in
    % Exp 1. Subject 1 of Experiment 2 is excluded for the same reason but
    % since it still responded enough, the tube was primed during its
    % session.
    % 21 was tested first in the first batch
    % 34 was tested first in the second batch
    % 39 was tested first in the third batch
    
    % gets max session number
    tmpSes = length(Timing.Subject(sub).Session);
    
    % change startses and stopses to look at different range of sessions
    startses= 1 ; %tmpSes - 5;%51;
    stopses= tmpSes; %tmpSes;%52;
    figure;
    RespRate=[];
    RespRateCons=[];
    emptyTrials = [];
    nTrials = [];
    rasterloc=0;
    Block = [];
    tri=0;
    for ses= startses:stopses
        tri=tri+1;
        
        TrialDelimits=Timing.Subject(sub).Session(ses).TrialDelimits;
        %I add a final element to trial delimits that will be excluded when calculating PI responses
        %this way we will be able to include the responses of final PI trial that is
        %omitted in the following lines due to the way the code is organized-(it gets the start of the lst PI but does not define an ending line for this)
        try 
            TrialDelimits(end+1,1:2) = [Timing.Subject(sub).Session(ses).PIEnd(end)  TrialDelimits(end,2)]; % %i fixed that now it will check until the end of the last PI end instead of the start of the last PI
        catch 
            TrialDelimits(end+1,1:2) = [NaN NaN];
        end
        
        Data=Timing.Subject(sub).Session(ses).Data;
        
        for FITrial=2:size(TrialDelimits,1)
            % figure;
            TempIn=[];
            TempOut=[];
            TempInTime=[];
            TempOutTime=[];
            rasterloc=rasterloc+1;
            
            TrayNPIn=Timing.Subject(sub).Session(ses).TrayNPIn;
            TrayNPOut=Timing.Subject(sub).Session(ses).TrayNPOut;
            
            TempIn=TrayNPIn(find(TrayNPIn>TrialDelimits(FITrial-1,1) & TrayNPIn<TrialDelimits(FITrial,1))); % finding np ins which are in between pi trial n and pi trial n+1 
            TempOut=TrayNPOut(find(TrayNPOut>TrialDelimits(FITrial-1,1) & TrayNPOut<TrialDelimits(FITrial,1))); % finding np outs which are in between pi trial n and pi trial n+1 
            
            
            if isempty(TempIn)==0
                if TempOut(1)<TempIn(1) % first in should be before first out-fixed
                    TempOut(1)=[];
                end
                
                if isempty(TempOut)==0
                    if TempOut(end)<TempIn(end) % last out cannot be smaller than lst in-fixed
                        TempIn(end)=[];
                    end
                end
            end
            
            TempInTime=[];
            TempOutTime=[];
            
            for ii=1:length(TempIn)
                TempInTime(ii)=str2num(Data{TempIn(ii),1})-str2num(Data{TrialDelimits(FITrial-1),1}); % this line normalizes the response times within the trial-second element is the start of the trial
            end
            TempInTime= TempInTime(TempInTime<45000); %previously code was taking all nose pokes from one pi trial to another which might include the response from a fi trial in between so I fixed this problem by adding a interval limit-45000 is the end of the trial for PI
            
            for ii=1:length(TempOut)
                TempOutTime(ii)=str2num(Data{TempOut(ii),1})-str2num(Data{TrialDelimits(FITrial-1),1}); % this line normalizes the response times within the trial 
            end
            TempOutTime= TempOutTime(TempOutTime<45000); %fixed same as in tempIn
            
            if isempty(TempOutTime)==0                  %EG
                if TempOutTime(end)<TempInTime(end) %checking in outs, after 45000 fix once more
                    TempInTime(end)=[];
                end
            end
            
            TempArray=[];
            TempArray=[floor(TempInTime*.001)+1; floor(TempOutTime*.001)+1]';  %converts time to seconds
            if isempty(TempArray)==0 
                TempInTime = []; TempInTime = TempArray(:,1);
                TempOutTime = []; TempOutTime = TempArray(:,2);
            end
            
            % subplot(1,3,1);
            % 
            % for ti=1:length(TempOutTime)
            %     hold on;
            %     plot([TempInTime(ti) TempOutTime(ti)],[rasterloc rasterloc],'b');
            % end
            % 
            % if length(TempInTime)>length(TempInTime)
            %     plot([TempInTime(ti) 45],[rasterloc rasterloc],'b');
            % end
            % xlabel('Trial Time','FontSize',12);
            % ylabel('Peak Interval Trial Number','FontSize',12);

            

            
            Timing.Subject(sub).Session(ses).PITrial(FITrial-1).Responses=[TempInTime, TempOutTime];
            
            if isempty(TempArray)==0 % so empty trials are not count
                
                Block(end + 1) = fix((ses - 1) / 5) + 1;
                RespRate(end+1,1:45)=0; 
                
                for Ri=1:length(TempArray(:,1)) % This way of doing it takes into account each response.
                    RespRate(end,TempArray(Ri,1):TempArray(Ri,2))=RespRate(end,TempArray(Ri,1):TempArray(Ri,2))+1;
                end
                
                RespRateCons(end+1,1:45)=0; 
                
                for Ri=1:length(TempArray(:,1)) % This way of doing it takes into account one resp in each bin.
                    RespRateCons(end,TempArray(Ri,1):TempArray(Ri,2))=1;
                end
            else
                emptyTrials(end + 1) = ses;
            end
        end
        nTrials(ses) = size(RespRate, 1);
    end
    %         xlim([0 45]);
    % 
    % 
    % subplot(1,3,2);
    % plot([15 15], [0 5], 'k:');
    % hold on;
    % w = plot(smooth(mean(RespRate),5));xlim([0 45]);
    % hold on;
    % x = plot(smooth(mean(RespRate),5)/max(smooth(mean(RespRate),5)),'r');xlim([0 45]);
    % xlabel('Trial Time','FontSize',12);
    % ylabel('Response Rate (all responses)','FontSize',12);
    % 
    % subplot(1,3,3);
    % plot([15 15], [0 1], 'k:');
    % hold on;
    % y = plot(smooth(mean(RespRateCons),5));xlim([0 45]);
    % hold on;
    % z = plot(smooth(mean(RespRateCons),5)/max(smooth(mean(RespRateCons),5)),'r');xlim([0 45]);
    % xlabel('Trial Time','FontSize',12);
    % ylabel('Response Rate (binary responses)','FontSize',12);
    % 
    % figureName = sprintf('Sub%d', sub);
    % 
    % print (figureName, '-dpng');
    % 
    % close gcf
    
%     %just to check data:
%     SteadyState(sub,1:45)=mean(RespRate);
%     SteadyStateCons(sub,1:45)=mean(RespRateCons);
    
    Timing.Subject(sub).RespRate=RespRate;
    Timing.Subject(sub).RespRateCons=RespRateCons;
    Timing.Subject(sub).emptyTrials = emptyTrials;
    Timing.Subject(sub).nTrials = nTrials;
    Timing.Subject(sub).Block = Block;
    
end
%%
% WTs = [28 29 30 36 37 40 42];
% TGs = [35 45 49 52 53 39 41 43 44];
% 
% for i=1:length(WTs);
%     WTPeaks(i,:)=mean(Timing.Subject(WTs(i)).RespRate);%27to29);
%     WTNormPeaks(i,:)=mean(Timing.Subject(WTs(i)).RespRateCons); %27to29);
% end
% for i=1:length(TGs);
%     TGPeaks(i,:)=mean(Timing.Subject(TGs(i)).RespRate); %27to29);
%     TGNormPeaks(i,:)=mean(Timing.Subject(TGs(i)).RespRateCons); %27to29);
% end
% 
% figure; plot(mean(WTPeaks)); hold on; plot(mean(TGPeaks),'r'); legend('WT', 'TG');
% figure; plot(mean(WTNormPeaks)); hold on; plot(mean(TGNormPeaks),'r'); legend('WT', 'TG');
%%
%  save PeakASD Timing
%%
% clear all
% load PeakASD Timing
%%
%For All Response:

    for sub= 1:sub_count; % [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59]%[21 22 23 28 29 30]
        dat=[];
        dat=Timing.Subject(sub).RespRate; %20to25; 
        for i=1:length(dat(:,1)) %for each trial
            data=[];
           data=dat(i,1:45);
           figure;
            try
                [StartStop] = getStartStop(data); % column 1 is start, 2 is stop, 3 is second start.
            catch 
                StartStop = [NaN NaN];
            end
            Timing.Subject(sub).Output(i,:)=StartStop; %added to save output
            close all;
            % Output refers to Ses20to25
        end
    end
    
%%
% save PeakASD Timing
%%
% MeanVals=[];
%  for sub= [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44]
%      MeanVals(sub,1:3)=mean(Timing.Subject(sub).Output); %20to25);
%      
%  end
% WTs = [28 29 30 36 37 40 42];
% TGs = [35 45 49 52 53 39 41 43 44];
% 
% ii=0;
% for sub=[28 29 30 36 37 40 42]
%     ii=ii+1;
%    WTVals(ii,:)= MeanVals(sub,1:3)
% end
% ii=0;
% for sub=[35 45 49 52 53 39 41 43 44]
%     ii=ii+1;
%    TGVals(ii,:)= MeanVals(sub,1:3)
% end
%%
%for Binary response
% clear all
% load PeakAD2 Timing

for sub= 1:sub_count; % [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59] %[1:7 9:18]
    dat=[];
    dat=Timing.Subject(sub).RespRateCons; %RespRateConsLTM;
    for i=1:length(dat(:,1)) %for each trial
        data=[];
        data=dat(i,1:45);
        figure;
        try
            [StartStop] = getStartStop(data); % column 1 is start, 2 is stop, 3 is second start.
        catch 
            StartStop = [NaN NaN];
        end
        Timing.Subject(sub).OutputCons(i,:)=StartStop;  %added to save output
        close all;
    end
end

%%
% save PeakASDnoTail Timing '-v7.3'
%%

% save PeakAD2 Timing
%%
% d=Timing.Subject(52).Output
% 
% figure;
% hold on;
% for i=1:length(d)
%     plot([d(i,1) d(i,2)],[i i],'r');
%     
% end
% plot([15 15],[1 i],'k--')
%%
% drop .Session so there is a compact version for import to R
% Timing2 = Timing;
% for i = 1:length(Timing2.Subject)
%     Timing2.Subject(i).Session = [];
% end
% save PeakASDShortNoTail Timing2

%%

for sub = 1:sub_count   
    block6 = Timing.Subject(sub).Block == 6;
    
    first_block6 = sum(block6 == 0) + 1;
    last_block6 = first_block6 + sum(block6 == 1) - 1;
    %RespRate = Timing.Subject(sub).RespRate(block6,:);
    RespRate = Timing.Subject(sub).RespRate(first_block6:last_block6,:);
    
    meanRespRate = mean(RespRate, 1); % averaged resp curve
    
    SmoothTest1(sub,:) = smooth(meanRespRate,5); % smoothed averaged resp curve
    
    normTest1(sub,:) = SmoothTest1(sub,:) / max(SmoothTest1(sub,:)); % normalized smoothed averaged resp curve
    
    pks = []; locs = []; w = []; p = []; % peak value % location % width % prominance
    
    data = normTest1(sub,:);
    
    figure
    
    %subplot(2,1,1)
    
    % plot(data,'black'); % plots data
    % 
    % hold on;
    % 
    % plot(SmoothTest1(sub,:), "blue");
    % 
    % plot(1:45, meanRespRate, "green");
    
    findpeaks(data); %marks the peak(s) on the plot
    
    [pks,locs, w, p] = findpeaks(data);
    
    PeakSmooth(sub).ParamsH1 = [pks;locs; w; p];
    
    PeakLoc = find(PeakSmooth(sub).ParamsH1(4,:) == max(PeakSmooth(sub).ParamsH1(4,:))) %get the one with highest prominance as the peak
    
    if length(PeakLoc) > 1
        peak_x = mean(PeakSmooth(sub).ParamsH1(2, PeakLoc));
    else
        peak_x = PeakSmooth(sub).ParamsH1(2, PeakLoc);
    end

    % xline(peak_x, 'red');
    % 
    % xlabel('Trial Time','FontSize',12);
    % ylabel('Response Rate','FontSize',12);
    % 
    % figureName = sprintf('Peaks_Sub%d', Timing.Subject(sub).SubjectID);
    % 
    % print (figureName, '-dpng')
    % 
    % hold off;
    
    Timing.Subject(sub).PeakLoc = PeakSmooth(sub).ParamsH1(2, PeakLoc);
end
%%
% make long format data
longTiming = [];
longTiming.subject = [];
longTiming.session = [];
longTiming.trial = [];
longTiming.second = [];
longTiming.responseRate = [];
longTiming.responseRateCon = [];
longTiming.outputStart = [];
longTiming.outputStop = [];
% longTiming.outputSecondStart = [];
longTiming.outputConStart = [];
longTiming.outputConStop = [];
% longTiming.outputConSecondStart = [];
longTiming.peaklocation1 = [];
longTiming.peaklocation2 = [];

for subj = 1:sub_count % [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59] %39
    subjID = Timing.Subject(subj).SubjectID;
    totalSessions = length(Timing.Subject(subj).Session);
    totalTrials = length(Timing.Subject(subj).RespRate);
    trialTracker = 1;
    for ses = 1:totalSessions
        if length(Timing.Subject(subj).Session(ses).Data) > 1
            if ses == 1
                sessionTrials = Timing.Subject(subj).nTrials(ses);
            else
                sessionTrials = Timing.Subject(subj).nTrials(ses) - Timing.Subject(subj).nTrials(ses - 1);
            end
            missingTrials = sum(Timing.Subject(subj).emptyTrials == ses);
            for trial = 1:(sessionTrials - missingTrials)
                for sec = 1:45
                    longTiming(end + 1).subject = subjID;
                    longTiming(end).session = ses;
                    longTiming(end).trial = trialTracker;
                    longTiming(end).second = sec;
                    longTiming(end).responseRate = Timing.Subject(subj).RespRate(trialTracker, sec);
                    longTiming(end).responseRateCon = Timing.Subject(subj).RespRateCons(trialTracker, sec);
                    longTiming(end).outputStart = Timing.Subject(subj).Output(trialTracker, 1);
                    longTiming(end).outputStop = Timing.Subject(subj).Output(trialTracker, 2);
%                     longTiming(end).outputSecondStart = Timing.Subject(subj).Output(trialTracker, 3);
                    longTiming(end).outputConStart = Timing.Subject(subj).OutputCons(trialTracker, 1);
                    longTiming(end).outputConStop = Timing.Subject(subj).OutputCons(trialTracker, 2);
%                     longTiming(end).outputConSecondStart = Timing.Subject(subj).OutputCons(trialTracker, 3);
                    if length(Timing.Subject(subj).PeakLoc) > 1
                        longTiming(end).peaklocation1 = Timing.Subject(subj).PeakLoc(1);
                        longTiming(end).peaklocation2 = Timing.Subject(subj).PeakLoc(2);
                    else
                        Timing.Subject(subj).PeakLoc
                        longTiming(end).peaklocation1 = Timing.Subject(subj).PeakLoc(1);
                        longTiming(end).peaklocation2 = 'NA';
                    end
                    longTiming(end).peak_m = mean(Timing.Subject(subj).PeakLoc);
                end
                trialTracker = trialTracker + 1;
            end
        end
    end
end
%%
delimits = 0;
subj = 39;
for i = 1:length(Timing.Subject(subj).Session)
    if length(Timing.Subject(subj).Session(i).Data) > 1
        delimits = delimits + length(Timing.Subject(subj).Session(i).PIEnd);
        if Timing.Subject(subj).Session(i).PIEnd(end) > 20*60
            Timing.Subject(subj).Session(i).PIEnd(end);
        end
    end
end
delimits
nTrials = length(Timing.Subject(subj).RespRate)
delimits - nTrials

count = 0
for i = 1:30
    count = count + length(Timing.Subject(subj).Session(i).PIIndex);
end
count

for i = 1:30
    if Timing.Subject(subj).Session(i).TrialDelimits(end,1) >20*60
        i
    end
    
end
%%

% save ASDlongTiming longTiming
%%
% write csv for use in R
struct2csv(longTiming, 'longDataNoTail.csv')    
%%
% change point analysis

% [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59]
% a(1).b='/Users/kyleroddick/University/Comps/Autism and Timing/ChangePointData'
% 
% cd(a(1).b);
% 
% for i = subList; % [28 29 30 35 36 37 45 49 52 53 39 40 41 42 43 44 49 55 57 59]
%     tmpData = Timing.Subject(i).Output(:,2);
%     for ii = 1:length(tmpData)
%         if tmpData(ii) < (15 * 2.5)
%             tmpData(ii) = 1;
%         else
%             tmpData(ii) = 0;
%         end
%     end
%     tmpID = num2str(Timing.Subject(i).SubjectID);
%     tmpPre = 'cpData';
%     tmpPost = '.csv';
%     tmpFileName = strcat(tmpPre, tmpID, tmpPost)
%     out = cpRL(tmpData, 5)
%     tmpOut = [];
%     tmpOut.subject = [];
%     tmpOut.a = [];
%     tmpOut.b = [];
%     tmpOut.c = [];
%     tmpOut.d = [];
%     for iii = 1:size(out,1)
%         tmpOut(end + 1).subject = tmpID;
%         tmpOut(end).a = out(iii, 1);
%         tmpOut(end).b = out(iii, 2);
%         tmpOut(end).c = out(iii, 3);
%         tmpOut(end).d = out(iii, 4);
%     end
% 
%     struct2csv(tmpOut, tmpFileName);
% end