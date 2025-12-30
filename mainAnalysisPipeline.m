%% 1. Load data

%load shockStabilityACC.mat
%load shockStabilityRSC.mat
%load tfcACC.mat
load tfcRSC.mat
%% 2. Extract data for each animal

poolMat = shockExp.calcium;
nCellsLong = shockExp.nCells;
nAnimals = length(shockExp.tensors);
freezeMat = shockExp.freezing;
dt = 20;
dtF = 15;
nSesh = 5;
for a=1:nAnimals
    tensors{a} = shockExp.tensors{a};
end

%IF YOU ARE DOING THE SHOCK STABILITY EXPERIMENT UNCOMMENT BELOW
% nTrials = 8.*ones(1,5);
% sessionNames = {'0hr','4hr','24hr','1wk','2wk'};
% sessions = {[1:8],[9:16],[17:24],[25:32],[33:40]};

%ELSE IF YOU ARE DOING TFC
nTrials = [8 8 6 6 6];
sessionNames = {'Hab','Acq','Ext1','Ext2','Ext3'};
sessions = {[1:8],[9:16],[17:22],[23:28],[29:34]};



%% 3a. Single cell dynamics over trials FOR SHOCK STABILITY
bl = 10; % zscore baseline
cap = 100; %ceiling for cell firing to compress variance
capForVisual = 50; %same but just to visualize
shockTime = 55; %time of shock

figure
for m=1:nSesh
for t=1:nTrials(m)
    orderCells=[];
    subplot(nSesh,nTrials(m),nTrials(m)*(m-1)+t)
    data = [];
    for a=1:nAnimals %concatenate tensors across animals
        data = [data tensors{a}(:,:,t+nTrials(m)*(m-1))];
    end
    data = (data-mean(data(1:bl*dt,:)))./std(data(1:bl*dt,:));
    data(data>cap)=cap;
    dat=data;
    dat(dat>capForVisual)=capForVisual;
    dat(dat<-capForVisual)=-capForVisual;
    shadedAreaPlot(dat,'r');
    xline(shockTime*dt)
    ylim([-1 6])
    yline(0)
    title(strcat('Trial-',num2str(t)))
    xlim([0 trialLength*dt])
    xticks([0 bl*dt shockTime*dt])
    xticklabels({'0','Shock','90'})
    xtickangle(45)
    xlabel('Seconds')
    ylabel('z-score')
    
    %get single cell dynamics around shock
    for c=1:size(data,2) 
        %take slice of data around shock
        time = [(shockTime-15)*dt:trialLength*dt];
        %smooth data
        signal_smooth=smooth(data(time,c),10);
        %take new baseline
        baseline=mean(signal_smooth(1:5*dt,:)); 
        %find the max values and index of data
        [peak_val, peak_idx] = max(signal_smooth); 
        t_peak = time(peak_idx); 
        %get time to peak
        time_to_peak = t_peak - time(1); 
        %fit 1/e decay
        decay_target = baseline + (peak_val - baseline) * exp(-.5);  
        %find time to half decay
        decay_idx = find(signal_smooth(peak_idx:end) <= decay_target, 1, 'first'); 
        if ~isempty(decay_idx)
            t_decay = time(peak_idx + decay_idx - 1);
            decay_time_constant = t_decay - t_peak;
        else
            decay_time_constant = NaN;  % decay not reached
        end
        %rise time (10% to 90% of peak)
        rise_10 = baseline + 0.1 * (peak_val - baseline);
        rise_90 = baseline + 0.9 * (peak_val - baseline);

        %find indices where signal crosses 10% and 90% of peak
        rise_10_idx = find(signal_smooth >= rise_10, 1, 'first');
        rise_90_idx = find(signal_smooth >= rise_90, 1, 'first');
        if ~isempty(rise_10_idx) && ~isempty(rise_90_idx)
            t_rise_10 = time(rise_10_idx);
            t_rise_90 = time(rise_90_idx);
            rise_time = t_rise_90 - t_rise_10;
        else
            rise_time = NaN;
        end
        
        %save dynamics for each cell, trial, and session
        amplitudes(c,t,m) = peak_val;
        peaktimes(c,t,m) = time_to_peak./20+40;
        decays(c,t,m) = decay_time_constant;
        rises(c,t,m) = rise_time;

    end
end

end

%prepare dynamic stats for export into Prism
for t=1:8
    amps(m,3*(t-1)+1:3*t) = [mean(amplitudes(:,t+32)) std(amplitudes(:,t+32)) length(amplitudes(:,t))];
    peaks(m,3*(t-1)+1:3*t) = [mean(peaktimes(:,t+32)) std(peaktimes(:,t+32)) length(peaktimes(:,t))];
    decs(m,3*(t-1)+1:3*t) = [nanmean(decays(:,t+32)) nanstd(decays(:,t+32)) length(decays(:,t))];
    rise(m,3*(t-1)+1:3*t) = [nanmean(rises(:,t+32)) nanstd(rises(:,t+32)) length(rises(:,t))];
end


%replace nans with medians
rises(isnan(rises))=nanmedian(rises,'all');
decays(isnan(decays))=nanmedian(decays,'all');
peaktimes(isnan(peaktimes))=nanmedian(peaktimes,'all');
amplitudes(isnan(amplitudes))=nanmedian(amplitudes,'all');

%identify and store neurons peaking before, during, and after shock
for m=1:nSesh
    for t=1:nTrials(m)
        preshockNeurons{t+nTrials(m)*(m-1)} = find(peaktimes(:,t,m)>shockTime-5 & peaktimes(:,t,m)<shockTime); 
        shockNeurons{t+nTrials(m)*(m-1)} = find(peaktimes(:,t,m)>=shockTime & peaktimes(:,t,m)<shockTime+3);
        postshockNeurons{t+nTrials(m)*(m-1)} = find(peaktimes(:,t,m)>=shockTime+3 & peaktimes(:,t,m)<shockTime+9);
        %calculate pooled fraction of these cells for each trial
        fracPool(t+nTrials(m)*(m-1),1) = length(preshockNeurons{t+8*(m-1)})./size(poolMat,2);
        fracPool(t+nTrials(m)*(m-1),2) = length(shockNeurons{t+8*(m-1)})./size(poolMat,2);
        fracPool(t+nTrials(m)*(m-1),3) = length(postshockNeurons{t+8*(m-1)})./size(poolMat,2);
    end
end
%plot pooled fractions
figure
plot(fracPool)
legend({'Pre-shock','Shock','Post-shock'})
title('Pooled fractions of cells')

%get fraction of these celltypes by animal and trial
x=0;
for a=1:nAnimals
    x=x+1;
    s1=cumsum(nCellsLong(1:a));
    s2=cumsum(nCellsLong(1:a+1));
    for t=1:nTrials(m)*nSesh
        fracAn(t,a,1) = mean(ismember([s1(end)+1:s2(end)],preshockNeurons{t}));
        fracAn(t,a,2) = mean(ismember([s1(end)+1:s2(end)],shockNeurons{t}));
        fracAn(t,a,3) = mean(ismember([s1(end)+1:s2(end)],postshockNeurons{t}));
    end
end

%plot mean ± SEM fraction cells
cols = {'c','m','g'};
figure
for p=1:3
shadedAreaPlot(fracAn(:,:,p),cols{p})
end

legend({'Pre-shock','Shock','Post-shock'})
title('Fractions of cells by animals')

%make cumulative density plots of each measured cell dynamic for the first
%session
cols = {'k','r',[255, 165, 0]./255,'y','g','c','b',[148,0,211]./255};
figure
subplot(1,4,1)
for t=1:nTrials(m)
h=cdfplot(amplitudes(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Amplitudes')
subplot(1,4,2)
for t=1:nTrials(m)
h=cdfplot(peaktimes(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Peak times')
subplot(1,4,3)
for t=1:nTrials(m)
h=cdfplot(rises(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Rise tau')
subplot(1,4,4)
for t=1:nTrials(m)
h=cdfplot(decays(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
xlim([0 80])
title('Decay tau')

%plot average trial activity of cell populations defined in last trial of 
%first session in each session
pops = {preshockNeurons{8}, shockNeurons{8}, postshockNeurons{8}}; %cell types from last trial of first session
pCols = {[255,154,0]./255,'r',[255,88,15]./255};

titles = {'0hrs', '4hrs', '24hrs', '1wk', '2wks'};
figure
for m=1:nSesh
    subplot(1,nSesh,m)
    data=[];
    for a=1:nAnimals
        data=[data tensors{a}(:,:,sessions{m})];
    end
    data=(data-mean(data(1:bl*dt,:,:)))./std(data(1:bl*dt,:,:));
    data(data>capForVisual)=capForVisual;
    data(data<-capForVisual)=-capForVisual;
    for p=1:length(pops)
    shadedAreaPlot(mean(data(:,pops{p},:),3),pCols{p});
    end
    yline(0)
    xline(shockTime*dt)
    ylim([-1 6])
    title(titles{m})
    ylabel('z-score')
    xlim([0 trialSamples])
    xticks([0 shockTime*dt trialSamples])
    xticklabels({0, 'Shock', 90})
end
% heatmaps of peaktimes in trial 1 vs. trial 8
trials = [1 8];
tis = {'Trial 1','Trial 8'};
figure
for t=1:2
    
    subplot(1,length(trials),t)
    mat=zeros(sum(nCellsLong(:,1)),trialLength);
    %collect peaktimes for each longitudinally registered cell
    for c=1:sum(nCellsLong(:,1))
    mat(c,round(peaktimes(c,trials(t)))) = 1;
    end
    %sort earliest to latest and take time slice used to calculate dynamics
    [B,I] = sort(peaktimes(:,trials(t)));
    data=mat(I,40:90);

    h=heatmap(data,'colorlimits',[0 .5]);
    grid off
    h.YDisplayLabels = [1 repmat("", 1, size(h.YData, 1)-2) size(data,1)];
    h.XDisplayLabels = [-15:35]'; %time clabels
    title(tis{t})
    ylabel('Cells')
    xlabel('Time wrt shock (sec)')
end




%calculate and plot overlaps
popNames = {'Pre-shock','Shock','Post-shock'};
for t=1:nTrials(m)*nSesh
    for r=1:nTrials(m)*nSesh
        overlaps(t,r,1) = mean(ismember(preshockNeurons{t},preshockNeurons{r}));
        overlaps(t,r,2) = mean(ismember(shockNeurons{t},shockNeurons{r}));
        overlaps(t,r,3) = mean(ismember(postshockNeurons{t},postshockNeurons{r}));
    end
end

figure
for p=1:3
subplot(1,3,p)
heatmap(overlaps(:,:,p),'colorlimits',[0 .6])
title(popNames{p})
end
colormap(palette('scheme',6))
grid off

%% 3b. Single cell dynamics over trials FOR TFC
bl = 10; % zscore baseline
cap = 100; %ceiling for cell firing to compress variance
capForVisual = 50; %same but just to visualize
shockTime = 55; %time of shock
m=2;
figure
for t=1:nTrials(m)
    orderCells=[];
    subplot(1,nTrials(m),t)
    data = [];
    for a=1:nAnimals %concatenate tensors across animals
        data = [data tensors{a}(:,:,t+nTrials(m)*(m-1))];
    end
    data = (data-mean(data(1:bl*dt,:)))./std(data(1:bl*dt,:));
    data(data>cap)=cap;
    dat=data;
    dat(dat>capForVisual)=capForVisual;
    dat(dat<-capForVisual)=-capForVisual;
    shadedAreaPlot(dat,'r');
    xline(shockTime*dt)
    ylim([-1 8])
    yline(0)
    title(strcat('Trial-',num2str(t)))
    xlim([0 trialLength*dt])
    xticks([0 bl*dt shockTime*dt])
    xticklabels({'0','Shock','90'})
    xtickangle(45)
    xlabel('Seconds')
    ylabel('z-score')
    
    %get single cell dynamics around shock
    for c=1:size(data,2) 
        %take slice of data around shock
        time = [(shockTime-15)*dt:trialLength*dt];
        %smooth data
        signal_smooth=smooth(data(time,c),10);
        %take new baseline
        baseline=mean(signal_smooth(1:5*dt,:)); 
        %find the max values and index of data
        [peak_val, peak_idx] = max(signal_smooth); 
        t_peak = time(peak_idx); 
        %get time to peak
        time_to_peak = t_peak - time(1); 
        %fit 1/e decay
        decay_target = baseline + (peak_val - baseline) * exp(-.5);  
        %find time to half decay
        decay_idx = find(signal_smooth(peak_idx:end) <= decay_target, 1, 'first'); 
        if ~isempty(decay_idx)
            t_decay = time(peak_idx + decay_idx - 1);
            decay_time_constant = t_decay - t_peak;
        else
            decay_time_constant = NaN;  % decay not reached
        end
        %rise time (10% to 90% of peak)
        rise_10 = baseline + 0.1 * (peak_val - baseline);
        rise_90 = baseline + 0.9 * (peak_val - baseline);

        %find indices where signal crosses 10% and 90% of peak
        rise_10_idx = find(signal_smooth >= rise_10, 1, 'first');
        rise_90_idx = find(signal_smooth >= rise_90, 1, 'first');
        if ~isempty(rise_10_idx) && ~isempty(rise_90_idx)
            t_rise_10 = time(rise_10_idx);
            t_rise_90 = time(rise_90_idx);
            rise_time = t_rise_90 - t_rise_10;
        else
            rise_time = NaN;
        end
        
        %save dynamics for each cell, trial, and session
        amplitudes(c,t) = peak_val;
        peaktimes(c,t) = time_to_peak./20+40;
        decays(c,t) = decay_time_constant;
        rises(c,t) = rise_time;

    end
end


%prepare dynamic stats for export into Prism
for t=1:8
    amps(m,3*(t-1)+1:3*t) = [mean(amplitudes(:,t)) std(amplitudes(:,t)) length(amplitudes(:,t))];
    peaks(m,3*(t-1)+1:3*t) = [mean(peaktimes(:,t)) std(peaktimes(:,t)) length(peaktimes(:,t))];
    decs(m,3*(t-1)+1:3*t) = [nanmean(decays(:,t)) nanstd(decays(:,t)) length(decays(:,t))];
    rise(m,3*(t-1)+1:3*t) = [nanmean(rises(:,t)) nanstd(rises(:,t)) length(rises(:,t))];
end


%replace nans with medians
rises(isnan(rises))=nanmedian(rises,'all');
decays(isnan(decays))=nanmedian(decays,'all');
peaktimes(isnan(peaktimes))=nanmedian(peaktimes,'all');
amplitudes(isnan(amplitudes))=nanmedian(amplitudes,'all');

m=2;
%identify and store neurons peaking before, during, and after shock
shockTime=55;
for t=1:nTrials(m)
    preshockNeurons{t} = find(peaktimes(:,t)>shockTime-5 & peaktimes(:,t)<shockTime); 
    shockNeurons{t} = find(peaktimes(:,t)>=shockTime & peaktimes(:,t)<shockTime+3);
    postshockNeurons{t} = find(peaktimes(:,t)>=shockTime+3 & peaktimes(:,t)<shockTime+9);
    %calculate pooled fraction of these cells for each trial
    fracPool(t,1) = length(preshockNeurons{t})./size(poolMat,2);
    fracPool(t,2) = length(shockNeurons{t})./size(poolMat,2);
    fracPool(t,3) = length(postshockNeurons{t})./size(poolMat,2);
end

%plot pooled fractions
figure
plot(fracPool)
legend({'Pre-shock','Shock','Post-shock'})
title('Pooled fractions of cells')
m=2;
%get fraction of these celltypes by animal and trial

x=0;
for a=1:nAnimals
    x=x+1;
    s1=cumsum(nCellsLong(1:a));
    s2=cumsum(nCellsLong(1:a+1));
    for t=1:nTrials(m)
        fracAn(t,a,1) = mean(ismember([s1(end)+1:s2(end)],preshockNeurons{t}));
        fracAn(t,a,2) = mean(ismember([s1(end)+1:s2(end)],shockNeurons{t}));
        fracAn(t,a,3) = mean(ismember([s1(end)+1:s2(end)],postshockNeurons{t}));
    end
end

%plot mean ± SEM fraction cells
cols = {'c','m','g'};
figure
for p=1:3
shadedAreaPlot(fracAn(:,:,p),cols{p})
end
legend({'Pre-shock','Shock','Post-shock'})
title('Fractions of cells by animals')

%make cumulative density plots of each measured cell dynamic for the first
%session
cols = {'k','r',[255, 165, 0]./255,'y','g','c','b',[148,0,211]./255};
figure
subplot(1,4,1)
for t=1:nTrials(m)
h=cdfplot(amplitudes(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Amplitudes')
subplot(1,4,2)
for t=1:nTrials(m)
h=cdfplot(peaktimes(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Peak times')
subplot(1,4,3)
for t=1:nTrials(m)
h=cdfplot(rises(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
title('Rise tau')
subplot(1,4,4)
for t=1:nTrials(m)
h=cdfplot(decays(:,t));
set(h, 'Color', cols{t}, 'LineWidth', 2);
hold on
end
xlim([0 80])
title('Decay tau')

%plot average trial activity of cell populations defined in last trial of 
%first session in each session

pops = {preshockNeurons{8}, shockNeurons{8}, postshockNeurons{8}}; %cell types from last trial of first session
pCols = {[255,154,0]./255,'r',[255,88,15]./255};

titles = {'0hrs', '4hrs', '24hrs', '1wk', '2wks'};
figure
for m=1:nSesh
    subplot(1,nSesh,m)
    data=[];
    for a=1:nAnimals
        data=[data tensors{a}(:,:,sessions{m})];
    end
    data=(data-mean(data(1:bl*dt,:,:)))./std(data(1:bl*dt,:,:));
    data(data>capForVisual)=capForVisual;
    data(data<-capForVisual)=-capForVisual;
    for p=1:length(pops)
    shadedAreaPlot(mean(data(:,pops{p},:),3),pCols{p});
    end
    yline(0)
    xline(shockTime*dt)
    ylim([-1 6])
    title(titles{m})
    ylabel('z-score')
    xlim([0 trialSamples])
    xticks([0 shockTime*dt trialSamples])
    xticklabels({0, 'Shock', 90})
end

% heatmaps of peaktimes in trial 1 vs. trial 8
trials = [1 8];
tis = {'Trial 1','Trial 8'};
figure
for t=1:2
    
    subplot(1,length(trials),t)
    mat=zeros(sum(nCellsLong(:,1)),trialLength);
    %collect peaktimes for each longitudinally registered cell
    for c=1:sum(nCellsLong(:,1))
    mat(c,round(peaktimes(c,trials(t)))) = 1;
    end
    %sort earliest to latest and take time slice used to calculate dynamics
    [B,I] = sort(peaktimes(:,trials(t)));
    data=mat(I,40:90);

    h=heatmap(data,'colorlimits',[0 .5]);
    grid off
    h.YDisplayLabels = [1 repmat("", 1, size(h.YData, 1)-2) size(data,1)];
    h.XDisplayLabels = [-15:35]'; %time clabels
    title(tis{t})
    ylabel('Cells')
    xlabel('Time wrt shock (sec)')
end



m=2
%calculate and plot overlaps
popNames = {'Pre-shock','Shock','Post-shock'};
for t=1:nTrials(m)
    for r=1:nTrials(m)
        overlaps(t,r,1) = mean(ismember(preshockNeurons{t},preshockNeurons{r}));
        overlaps(t,r,2) = mean(ismember(shockNeurons{t},shockNeurons{r}));
        overlaps(t,r,3) = mean(ismember(postshockNeurons{t},postshockNeurons{r}));
    end
end

figure
for p=1:3
subplot(1,3,p)
heatmap(overlaps(:,:,p),'colorlimits',[0 .6])
title(popNames{p})
end
colormap(palette('scheme',6))
grid off

%% 4. Calculate and plot neural synchrony
%arrange cell types
pops = {preshockNeurons{8} shockNeurons{8} postshockNeurons{8}};

%concatenate tensors of activity over animals
d=[];
for a=1:nAnimals
    d = [d tensors{a}];
end

%plot and calculate change in synchrony over 0.5 sec intervals
bin=10;
figure
for m=1:nSesh
    subplot(1,nSesh,m)
    clear corrs
    for p=1:3
        for l=1:nTrials(m)
            data=movmean(d(:,pops{p},sessions{m}(l)),dt);
            
            % generate correlation matrices between cells over 0.5 second
            % bins and only keep values under the diagonal. a neuron's
            % synchrony is its mean corr with all other neurons in its
            % group
            for t=1:trialSamples/bin-1
                r = corr(data(bin*(t-1)+1:bin*(t-1)+dt,:)');
                rs=r(tril(true(size(r)), -1));
                corrs(t,l)=mean(rs);
            end
        end
        corrcell{m,p} = corrs;
        toPlot=((corrs(:,:)-mean(corrs(1:bin,:)))./std(corrs(1:bin,:)));
        shadedAreaPlot2(movmean(toPlot,5),length(pops{p}),pCols{p})
        hold on
        
        xline(shockTime*2)
        yline(0,'--')

    end
end

%prepare stats for export
for p=1:3
    for m=1:5
        dat=corrcell{m,p};
        dat=((dat(:,:)-mean(dat(1:bin,:)))./std(dat(1:bin,:)));
        %shock synchrony
        stats(m,p,1,1:3) = [mean(dat(55*2+1:57*2,:),'all') std(mean(dat(55*2+1:57*2,:))) length(pops{p})];
        %pre-shock synchrony
        stats(m,p,2,1:3) = [mean(dat(50*2+1:55*2,:),'all') std(mean(dat(50*2+1:60*2,:))) length(pops{p})];
        %post-shock synchrony
        stats(m,p,3,1:3) = [mean(dat(58*2+1:63*2,:),'all') std(mean(dat(58*2+1:63*2,:))) length(pops{p})];
        %whole window synchrony
        stats(m,p,4,1:3) = [mean(dat(50*2+1:63*2,:),'all') std(mean(dat(50*2+1:63*2,:))) length(pops{p})];
        
        stats2(m,p,1,:) = [mean(dat(56*2+1:58*2,:))];
        stats2(m,p,2,:) = [mean(dat(50*2+1:55*2,:))];
        stats2(m,p,3,:) = [mean(dat(58*2+1:63*2,:))];
        stats2(m,p,4,:) = [mean(dat(50*2+1:63*2,:))];
    end
end

%% 5. Calculate and visualize PCA
pop=[1:size(poolMat,2)];

%uncomment these lines to knock out cell type of interest in calculations
%pop(preshockNeurons{8})=[];
%pop(shockNeurons{8})=[];
%pop(postshockNeurons{8})=[];

pc1=1;
pc2=2;

%calculate and plot 
figure
for m=1:nSesh
    subplot(1,nSesh,m)
    data = zscore(poolMat((sessions{m}(1)-1)*trialLength*dt+1:sessions{m}(end)*trialLength*dt,pop));
    
    data = zscore(poolMat((sessions{m}(1)-1)*trialLength*dt+1:sessions{m}(end)*trialLength*dt,pop));
    [co,sc,~,~,ex] = pca(data);
    scs{m}=sc;
    cos{m}=co;
    c = linspace(0,1,length(sessions{m}));

    x=[];
    y=[];

    for t=1:length(sessions{m})
        x=[x mean(sc(trialLength*dt*(t-1)+1:trialLength*dt*t,pc1))];
        y=[y mean(sc(trialLength*dt*(t-1)+1:trialLength*dt*t,pc2))];

        hold on
        bl = 1:10*dt;
        sti = 10*dt+1:35*dt;
        scatter(mean(sc(90*(t-1)*dt+sti,pc1)),mean(sc(90*(t-1)*dt+sti,pc2)),50,[.5 .5 .5],'filled')
        hold on
        tr = 35*dt+1:50*dt;
        pre = 50*dt+1:55*dt;
        shock = 55*dt+1:57*dt;
        scatter(mean(sc(90*(t-1)*dt+shock,pc1)),mean(sc(90*(t-1)*dt+shock,pc2)),50,'r','filled')
        hold on
        post = 57*dt+1:90*dt;
        triTimes = {bl, sti, tr, pre, shock,  post};
    end
    hold on
    scatter(x,y,150,c,'filled')
    hold on
    plot(x,y)
    xlabel('PC1')
    ylabel('PC2')
    ylim([-15 15])
    xlim([-15 15])
    zlim([-15 15])
    title(sessionNames{m})
end



%calcualte and plot distance traveled across trials for first two PCs;
%switch between them below
dim=1:2;

figure
for m=1:nSesh
    distTraveled=zeros(length(sessions{m}));
    subplot(1,5,m)
    for t=1:length(sessions{m})
        d1 = mean(scs{m}(trialLength*dt*(t-1)+1:trialLength*dt*(t),dim));
        for r=1:length(sessions{m})
            d2 = mean(scs{m}(trialLength*dt*(r-1)+1:trialLength*dt*(r),dim));
            distTraveled(t,r)=pdist2(d1, d2)./sqrt(2);
        end
    end
    distancesTraveled{m} = distTraveled;
    heatmap(distTraveled,'colorlimits',[0 10])
    colormap(palette('scheme',1))
    title(sessionNames{m})
    xlabel('Trials')
    ylabel('Trials')
end


% calculate and plot distance traveled within vs. between trials
 sessCols={'k','r','y','g','b'}
for m=1:nSesh
    distTraveled=zeros(length(sessions{m}),1);
    subplot(1,5,m)
    for t=1:length(sessions{m})
        d=scs{m}(trialLength*dt*(t-1)+1:trialLength*dt*(t),dim);
        d1 = mean(d(1:trialLength/2*dt,:));
        d2 = mean(d(trialLength/2*dt+1:trialLength*dt,:));
        distTraveled(t,1)=pdist2(d1, d2);
        neighborDistances{m} = distTraveled./sqrt(2);
    end
end

figure
for m=1:nSesh
    d=distancesTraveled{m}(tril(true(size(distancesTraveled{m})), -1));
    scatter((2*(m-1)+1).*ones(length(d),1)+1,d,30,sessCols{m},'filled')
    hold on
    d1=neighborDistances{m};
    scatter(2*m.*ones(length(d1),1)+1,d1,30,sessCols{m})  
    dmats{m,1} = d;
    dmats{m,2} = d1;
end
xlim([0 13])
xticks([2 4 6 8 10])
xticklabels(sessionNames)
ylabel('Distance traveled')
title('Distance traveled within (empty) and between (filled) trials')

for m=1:nSesh
    dmats{m,1}(end+1:length(dmats{1,1})) = nan;
    dmats{m,2}(end+1:length(dmats{1,2})) = nan;
end
betweenDistanceTable = table(dmats{1,1}, dmats{2,1}, dmats{3,1}, dmats{4,1}, dmats{5,1}, 'VariableNames', sessionNames);
withinDistanceTable = table(dmats{1,2}, dmats{2,2}, dmats{3,2}, dmats{4,2}, dmats{5,2}, 'VariableNames', sessionNames);
writetable(betweenDistanceTable,'betdist.csv')
writetable(withinDistanceTable,'witdist.csv')


%calculate and plot distance as a function of temporal lag
timeLagDistances=cell(nTrials(m)-1);
for m=1:nSesh
    for t = 1:length(sessions{m})-1
        timeLagDistances{t}=[timeLagDistances{t}; diag(distancesTraveled{m},-t)];
    end
end

m=2;
figure
for t=1:nTrials(m)-1
bar(t,mean(timeLagDistances{t}))
hold on
errorbar(t,mean(timeLagDistances{t}),std(timeLagDistances{t})./sqrt(length(timeLagDistances{t})-1))
hold on
end
xticks([1:7])
xlabel('Time lag between trials')
ylabel('Distances')
title('Distance increases as a function of temporal distance')


dim=1:2;
for m=1:nSesh
    sp=[];
    for s=1:length(dim)
        v= diff(scs{m}(:,dim(s)));
        sp=[sp v.^2];
    end
    speed{m} = sqrt(sum(sp,2))./sqrt(size(sp,2));
    
    accel{m} = diff(speed{m});
    
end

%calculate average pooled trajectory speed and acceleration in each trial period
%for each trial and animal
for m=1:nSesh
    for t=1:nTrials(m)
        dat=speed{m}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
        speedMat{m}(:,t) = dat;
        %dat=movmean(dat,20); 
        speedTrialPeriods{m}(1,t) = mean(dat(1*dt:5*dt-2,:));
        speedTrialPeriods{m}(2,t) = mean(dat(8*dt:13*dt-2,:));
        speedTrialPeriods{m}(3,t) = mean(dat(33*dt:38*dt-2,:));
        speedTrialPeriods{m}(4,t) = mean(dat(50*dt:55*dt-2,:));
        speedTrialPeriods{m}(5,t) = mean(dat(55*dt:57*dt-2,:));
        speedTrialPeriods{m}(6,t) = mean(dat(57*dt:62*dt-2,:));
        
        dat=accel{m}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
        accelMat{m}(:,t) = dat;
        %dat=movmean(dat,20); 
        accelTrialPeriods{m}(1,t) = mean(dat(1*dt:5*dt-2,:));
        accelTrialPeriods{m}(2,t) = mean(dat(8*dt:13*dt-2,:));
        accelTrialPeriods{m}(3,t) = mean(dat(33*dt:38*dt-2,:));
        accelTrialPeriods{m}(4,t) = mean(dat(50*dt:55*dt-2,:));
        accelTrialPeriods{m}(5,t) = mean(dat(55*dt:57*dt-2,:));
        accelTrialPeriods{m}(6,t) = mean(dat(57*dt:62*dt-2,:));
        
    end
end

for m=1:nSesh
    sps{m}=speedTrialPeriods{m}.*10;
    ac{m}=accelTrialPeriods{m}.*1000;
    sps{m}(:,end+1:8) =nan;
    ac{m}(:,end+1:8) =nan;
end

sp2 = [];
ac2 = [];
for m=1:nSesh
    s = [];
    a = [];
    for t=1:6
        s = [s sps{m}(t,:)];
        a = [a ac{m}(t,:)];
    end
    sp2=[sp2; s];
    ac2 = [ac2; a];
end
speedTable = table([sp2]);
accTable = table([ac2]);

writetable(speedTable,'speedtable.csv')
writetable(accTable,'acctable.csv')
%% 6. plot average speed / acceleration over trials. choose your parameters
for m=1:nSesh
    for t=1:nTrials(m)
        dat=speed{m}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
        speedMat{m}(:,t) = dat;
        dat=movmean(dat,20); 
        speedTrialPeriodsSmooth{m}(1,t) = mean(dat(1*dt:5*dt-2,:));
        speedTrialPeriodsSmooth{m}(2,t) = mean(dat(8*dt:13*dt-2,:));
        speedTrialPeriodsSmooth{m}(3,t) = mean(dat(33*dt:38*dt-2,:));
        speedTrialPeriodsSmooth{m}(4,t) = mean(dat(50*dt:55*dt-2,:));
        speedTrialPeriodsSmooth{m}(5,t) = mean(dat(55*dt:57*dt-2,:));
        speedTrialPeriodsSmooth{m}(6,t) = mean(dat(57*dt:62*dt-2,:));
        
        dat=accel{m}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
        accelMat{m}(:,t) = dat;
        dat=movmean(dat,20); 
        accelTrialPeriodsSmooth{m}(1,t) = mean(dat(1*dt:5*dt-2,:));
        accelTrialPeriodsSmooth{m}(2,t) = mean(dat(8*dt:13*dt-2,:));
        accelTrialPeriodsSmooth{m}(3,t) = mean(dat(33*dt:38*dt-2,:));
        accelTrialPeriodsSmooth{m}(4,t) = mean(dat(50*dt:55*dt-2,:));
        accelTrialPeriodsSmooth{m}(5,t) = mean(dat(55*dt:57*dt-2,:));
        accelTrialPeriodsSmooth{m}(6,t) = mean(dat(57*dt:62*dt-2,:));
        
    end
end

toPlot = "s";

if toPlot=="s"
    datas = speedTrialPeriodsSmooth;
    %if stab
    %ys = [0 .8];
    %if TFC
    ys = [0 .3];
else
    datas = accelTrialPeriodsSmooth;
    %ys=[-0.02 0.02];
    %if TFC
    ys=[-0.005 0.005];
end

figure
for m=1:5
    subplot(5,1,m)
    data=datas{m};
    
    dat=reshape(data,size(data,1)*size(data,2),size(data,3));
    
    shadedAreaPlot(dat,[.5 .5 .5])
    
    for t=1:length(sessions{m})
        hold on
        scatter(6*(t-1)+1,mean(dat(6*(t-1)+1,:),2),30,'k','filled')
        hold on
        scatter(6*(t-1)+2,mean(dat(6*(t-1)+2,:)),30,'b','filled')
        hold on
        scatter(6*(t-1)+3,mean(dat(6*(t-1)+3,:)),30,'c','filled')
        hold on
        scatter(6*(t-1)+4,mean(dat(6*(t-1)+4,:)),30,cols{3},'filled')
        hold on
        scatter(6*(t-1)+5,mean(dat(6*(t-1)+5,:)),30,'r','filled')
        hold on
        scatter(6*(t-1)+6,mean(dat(6*(t-1)+6,:)),30,cols{3},'filled')
        hold on
        xline(6.5+6*(t-1))
    end
%ylim([0 .6])
if m<3
xticks([1 7 13 19 25 31 37 43])
xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Trial 6','Trial 7','Trial 8'})
else
    xticks([1 7 13 19 25 31])
xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Trial 6'})
end
xlim([0 49])
ylabel('Speed')
title(strcat('Average trajectory speed around events in',sessionNames{m}))
if m==5
    legend({'Mean trajectory speed','SEM','Baseline','Tone','Trace','Pre-shock','Shock','Post-shock'})
end
ylim(ys)
end
%% 7. Calculate trajectories for individual animals; same process as for pooled

for m=1:nSesh
    for a=1:nAnimals
        s1=cumsum(nCellsLong(1:a,1));
        s2=cumsum(nCellsLong(1:a+1,1));
        range=s1(end)+1:s2(end);
        
        data = zscore(poolMat((sessions{m}(1)-1)*trialLength*dt+1:sessions{m}(end)*trialLength*dt,range));
        [co,sc,~,~,ex] = pca(data);
        scsA{m,a}=sc;
    end
end

% trajectory speed
for dim=1:2
    if dim==1
        dimname='pc1';
    else
        dimname='pc2';
    end
    for a=1:nAnimals
        for m=1:5
            
            
            if size(scsA{m,a},2)<dim
                speedA{m,a}=nans(size(speedA{m,1}));
                accelA{m,a}=nans(size(accelA{m,1}));
            else
                v= diff(scsA{m,a}(:,1:2));
                sp=v.^2;
                
                
                speedA{m,a} = sqrt(sum(sp,2))./sqrt(size(sp,2)); %sqrt(vx.^2 + vy.^2);
                
                accelA{m,a} = diff(speedA{m,a});
            end
            
        end
    end
    
    for a=1:nAnimals
        for m=1:nSesh
            
            for t=1:length(sessions{m})
                dat=speedA{m,a}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
                speedMatA{m}(:,t,a) = dat;
                dat=movmean(dat,20);
                speedTrialPeriodsA{m}(1,t,a) = mean(dat(1*dt:5*dt-2,:));
                speedTrialPeriodsA{m}(2,t,a) = mean(dat(8*dt:13*dt-2,:));
                speedTrialPeriodsA{m}(3,t,a) = mean(dat(33*dt:38*dt-2,:));
                speedTrialPeriodsA{m}(4,t,a) = mean(dat(50*dt:55*dt-2,:));
                speedTrialPeriodsA{m}(5,t,a) = mean(dat(55*dt:57*dt-2,:));
                speedTrialPeriodsA{m}(6,t,a) = mean(dat(57*dt:62*dt-2,:));
                
                
                dat=accelA{m,a}(trialLength*dt*(t-1)+1:trialLength*dt*t-2,:);
                accelMatA{m}(:,t,a) = dat;
                
                
                accelTrialPeriodsA{m}(1,t,a) = mean(dat(1*dt:5*dt-2,:));
                accelTrialPeriodsA{m}(2,t,a) = mean(dat(8*dt:13*dt-2,:));
                accelTrialPeriodsA{m}(3,t,a) = mean(dat(33*dt:38*dt-2,:));
                accelTrialPeriodsA{m}(4,t,a) = mean(dat(50*dt:55*dt-2,:));
                accelTrialPeriodsA{m}(5,t,a) = mean(dat(55*dt:57*dt-2,:));
                accelTrialPeriodsA{m}(6,t,a) = mean(dat(57*dt:62*dt-2,:));
                if m==4
                    dat=movmean(dat,20); %smooth for corr
                    smoothAccel(t,a) = mean(dat(55*dt:57*dt-2,:));
                end
            end
        end
    end
    
    
    
    ss=[];
    aa=[];
    for m=1:nSesh
        s=[];
        a=[];
        for t=1:6
            
            s=[s; squeeze(nanmean(speedTrialPeriodsA{m}(t,:,:),2)).*10];
            a=[a; squeeze(nanmean(accelTrialPeriodsA{m}(t,:,:),2)).*1000'];
            
        end
        ss=[ss; s'];
        aa=[aa; a'];
    end
    sss{dim}=ss;
    aaa{dim}=aa;
    for t=1:6
        accelIndivs(1+nAnimals*(dim-1):nAnimals*dim,:,t) = aa(:,nAnimals*(t-1)+1:nAnimals*t)';
    end
    csvwrite(strcat(dimname,'Accel.csv'),aa)
    aas(:,:,dim) = aa';
end

figure
for m=1:nSesh
subplot(1,nSesh,m)
shadedAreaPlot(squeeze(accelIndivs(1:10,m,:))','r')
hold on
shadedAreaPlot(squeeze(accelIndivs(10+1:10*2,m,:))','b')
title(sessionNames{m})
if m==5
legend({'PC1','','PC2',''})
end
yline(0)
ylim([-1.5 2])
xlim([0 7])
xticks([1:6])
xticklabels({'Baseline','Tone','Trace','Pre-shock','Shock','Post-shock'})
xtickangle(45)
end

%% 8. For TFC, analyze & visualize behavior
x=0;
dtF = 15;
%Organize freezing in same way as trajectory accelerations
for a=1:nAnimals
    data = freezeMat(:,a);
    for m=1:nSesh
        dat=data(trialLength*dtF*(sessions{m}(1)-1)+1:trialLength*dtF*(sessions{m}(end)),1);
        
        for t=1:length(sessions{m})
            d = dat(trialLength*dtF*(t-1)+1:trialLength*dtF*t,:);
            freMat{m}(:,t,a) = d;
            freDat{m}(t,a,1) = mean(d(5*dtF+1:10*dtF,:));
            freDat{m}(t,a,2) = mean(d(10*dtF+1:35*dtF,:));
            freDat{m}(t,a,3) = mean(d(35*dtF+1:50*dtF,:));
            freDat{m}(t,a,4) = mean(d(50*dtF+1:55*dtF,:));
            freDat{m}(t,a,5) = mean(d(55*dtF+1:57*dtF,:));
            freDat{m}(t,a,6) = mean(d(57*dtF+1:90*dtF,:));
        end
        meanFreeze(a,m) = mean(dat);
    end
end
%calculate fear recall and extinction
fearRecall = mean(mean(freDat{3}(1:2,:,2:3),3),1)'-mean(mean(freDat{1}(end-1:end,:,2:3),3),1)'; %unsmoothed both
fearExtinction = mean(mean(freDat{5}(1:2,:,2:3),3),1)'-mean(mean(freDat{3}(end-1:end,:,2:3),3),1)'; %smoothed both

%visualize
figure
plot([ones(nAnimals,1) 2.*ones(nAnimals,1)]',[fearRecall fearExtinction]','k')
hold on
scatter(ones(nAnimals,1),fearRecall,50,'r','filled')
hold on
scatter(2.*ones(nAnimals,1),fearExtinction,50,'c','filled')
xlim([0 3])
xticks([1 2])
yline(0)
xticklabels({'early Ext1-late Hab', 'early Ext3-late Ext1'})
title('∆ Freezing')

%relate to trajectory accelerations
figure
subplot(121)
scatter(squeeze(mean(accelTrialPeriodsA{2}(4,:,:),2)).*1000,fearRecall,50,'r','filled')
lsline
xlabel('Mean acceleration during Acq pre-shock')
ylabel('Fear recall')
subplot(122)
scatter(mean(smoothAccel)'.*1000,fearExtinction,50,'c','filled')
lsline
xlabel('Mean acceleration during Ext2 shock omission')
ylabel('Fear Extinction')



