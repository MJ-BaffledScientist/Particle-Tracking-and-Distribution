%calcFPP - orignal by SS Rogers, modified for parallel calculation
%
% usage [S,FPP]=calcFPP_parallel(tracks,framerate,pixscale,Lscale,interactive)
%
% Calculates First Passage Probability (FPP) as function of slowness S, as
% described in:
% "The first passage probability of intracellular particle trafficking"
% Salman S. Rogers, Neftali Flores-Rodriguez, Victoria J. Allan,
% Philip G. Woodman and Thomas A. Waigh
% Phys. Chem. Chem. Phys., 2010, 12, 3753-3761
%
% input: tracks - cell array of tracks
%        framerate - frames per second
%        pixscale - microns per pixel
%        Lscale - array of length scales at which to calculate FPP
%        interactive - set to 1/0 to plot output or not
%   N.B. to run in parallel, the user is required to start matlabpool
%   before executing the scipt.

function [S,FPP,counts]=calcFPP_parallel(tr_raw,framerate,pixscale,Lscale,interactive)

tracktoggle=true(size(tr_raw));

colours=jet(length(Lscale));
NEvents=0;
end2end=zeros(size(tr_raw));

if interactive
    figure
    hold on
end

for n=find(tracktoggle)
    NEvents=NEvents+length(tr_raw{n}(:,1));
    end2end(n)=sqrt(sum(diff(tr_raw{n}([1 end],1:2)).^2));
end
S=cell(size(Lscale));
FPP=cell(size(Lscale));

counts = zeros(1,length(Lscale));

parfor L=1:length(Lscale)
    timeWithin=cell(length(tr_raw),1);
    Lscale2=(Lscale(L)/pixscale)^2;
    ttoggleLhit=false(length(tr_raw),1);
    for n=find(tracktoggle)
        timeWithin{n}=NaN(size(tr_raw{n}(:,1)));
        for t=1:size(tr_raw{n},1)
            dr2end=(tr_raw{n}(end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(end,2)-tr_raw{n}(t,2)).^2;
            if dr2end>Lscale2 %if particle has made passage

                dr2=(tr_raw{n}(t:end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(t:end,2)-tr_raw{n}(t,2)).^2;
                indices=cumsum(~(dr2<Lscale2));
                
                if indices(end)
                    timeWithin{n}(t)=sum(~indices);
                    ttoggleLhit(n)=true;
                end
            end
        end
    end
    timeWithin=cell2mat(timeWithin);
    timeWithin=timeWithin(~isnan(timeWithin)&(timeWithin>0));
    %     timeWithin=timeWithin/framerate; % first passage time
%     timeWithin=timeWithin/Lscale(L)/framerate; % first passage inv vel
        timeWithin=1./(timeWithin/Lscale(L)/framerate); % first passage vel
    
    if ~isempty(timeWithin)
        switch 'a'
            case 'a'
                nbars=28;
                points=log(min(timeWithin)/2)+(0:1/nbars:1)*(log(max(timeWithin)*2)-log(min(timeWithin)/2));
                edges=points(1:2:end);
                S{L}=exp(points(2:2:end));
                N=histc(log(timeWithin),edges);
                N=N(:);
                %%%%%%%%%%%%%%%%%%%%
                FPP{L}=N(1:end-1)'/NEvents; %normalise with total number of time points
                disp(['L=' num2str(Lscale(L)) ', number of passages =' num2str(length(timeWithin))])
                counts(L) = length(timeWithin);
                %%%%%%%%%%%%%%%%%%%%
                FPP{L}=FPP{L}./diff(exp(edges)); %rescale in terms of prob density
            case 'b'
                [FPP{L},S{L}]=hist(timeWithin,20); disp('test')
                if interactive
                    plot(S{L},FPP{L},'displayname',['L = ' num2str(Lscale(L)) ' \mum'],'color',colours(L,:))
                end
        end
    else
        S{L}=[];
        FPP{L}=[];
    end
    disp(L)
end

if interactive
    for i = 1:length(FPP)
        plot(S{i},FPP{i},'DisplayName',['L = ' num2str(Lscale(i)) ' \mum'],'color',colours(i,:),'linewidth',2)
    end
%     xlabel('Inverse speed S  (s\mum^{-1})','fontsize',14)
 xlabel('Velocity (\mums^{-1})','fontsize',14)
    ylabel('Probability density (s{\mu}m^{-1})','fontsize',14)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    legend('location','best')
end