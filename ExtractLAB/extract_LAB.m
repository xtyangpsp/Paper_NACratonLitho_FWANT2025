close all;
ite_nm = 'ite_0.05deg_07';
velocitytag='S'; % 'P' for P velocities, 'PR' for Poisson's ratios.
% idx = 0; % 0, absolute velocity; 1, velocity perturbation
% savefigtag=1;
% figlabel={'A. ','B. ','C. ','D. ','E. ','F. ','G. ','H. ','I. ','J. '};
figlabel={'a. ','b. ','c. ','e. ','e. ','f. ','g. ','h. ','i. ','j. '};
%read previous model
fnm_conf=['./SeisFD3D.conf_' ite_nm];
dir_coord=['./updated_input_' ite_nm];
dir_media=['./updated_input_' ite_nm];
disp(['Read model... ' dir_media]);

id = 0; subs=[1,1,1];subc=[-1,-1,-1];subt=[1,1,1];
[snapinfo]=locate_snap(fnm_conf,id,'start',subs,'count',subc,'stride',subt);
[XSIM,YSIM,ZSIM]=gather_coord(snapinfo,'coorddir',dir_coord);
% convert from radian to degrees
XSIM=90-XSIM*180/pi; %latitude
YSIM=YSIM*180/pi;

%define the area of plot (exclude pmls)
npml=12; %number of pml layers
minlat=XSIM(end,1,end);maxlat=XSIM(1,1,end);
minlon=YSIM(1,1,end);maxlon=YSIM(1,end,end);

mrh=gather_media(snapinfo,'rho','mediadir',dir_media);
mmu=gather_media(snapinfo,'mu','mediadir',dir_media);
mla=gather_media(snapinfo,'lambda','mediadir',dir_media);
mvp=((mla+2*mmu)./mrh).^0.5;
mvs0=(mmu./mrh).^0.5/1000;
mvs=mvs0(npml:end-npml,npml:end-npml,:);
x=squeeze(YSIM(1,npml:end-npml,1))-360;
y=squeeze(XSIM(npml:end-npml,1,1));
% x=YSIM(1,:,1)-360;
% y=XSIM(:,1,1);
z=6371-squeeze(ZSIM(1,1,:))/1000;
load ../00masterdatafiles/us_states.mat;
canada_provs=shaperead('PROVINCE.SHP','UseGeoCoords',true);
%%
amask=nan(length(y),length(x));
% set a box to nan
load('RayCoverOutline_ite_0.05deg_07_60-100s_cutoff5.mat');
for i=1:size(amask,1)
    clear id00;
    id00=inpolygon(x,y(i)*ones(size(amask,2),1),raycover.data(:,1),...
            raycover.data(:,2));
    amask(i,id00)=1;
end
amask3d=nan(size(mvs));
for k=1:length(z)
    amask3d(:,:,k)=amask;
end
vplot=mvs.*amask3d;

%% LAB depth
%Detection and defination of low-velocity layer:
% Step-1: search for the first velocity drop below a certain top depth (e.g., 40 km). 
% The velocity drop is defined as negative gradient measuring from the top
% downward.
% Step-2: set acceptance threshold. Set a threshold for this
% gradient or velocity value.

disp('Detecting LAB ...')
myzlimlab=[40,250];
dz=0.5;
vgrad_min=0.01;
vdrop_min=0.3;
vmin=4.6;
zinterp=min(z):dz:max(z);

zlab=nan(size(squeeze(vplot(:,:,1))));
% vgrad_all=nan(size(squeeze(vplot(:,:,1))));
depthidxvlab=find(zinterp>=myzlimlab(1) & zinterp<=myzlimlab(2));
zsub=zinterp(depthidxvlab);
plottest=0;
if plottest
    figure('Position',[400 400 600 800]); hold on;
    grid on;
    axis on;
    box on;
    xlabel('Vs (km/s)')
    ylabel('depth (km)')
    xlim([3.8,5.4])
    ylim(myzlimlab);
    set(gca,'fontsize',14,'ydir','reverse')
end
subsample_step=1;
for i=1:subsample_step:size(vplot,1)
    for j =1:subsample_step:size(vplot,2)
        vz=squeeze(vplot(i,j,:));
        if isnan(sum(vz))
            continue;
        end
        vzinterp=interp1(z,vz,zinterp);
        vsub=vzinterp(depthidxvlab);
        [vdrop,firstneg,firstpos,idxtarget,vtarget]=find_firstdrop(vsub);
        if isnan(vdrop)
%             idxtarget_final=length(vsub);
            disp('no 1st drop')
            continue;
        else
            if firstpos < length(vsub)-1 &&  vdrop < vdrop_min && max(vsub(firstpos + 1:end)) >= vmin
                %search second drop
                [vdrop2,firstneg2,firstpos2,idxtarget2,vtarget2]=find_firstdrop(vsub(firstpos + 1:end));
                if isnan(vdrop2)
                    disp('no second drop')
                    idxtarget_final=idxtarget;
                else
                    idxtarget_final=firstpos + idxtarget2;
                end
            else
                idxtarget_final=idxtarget;
            end

        end
        vlab=vsub(idxtarget_final);
        zlabtemp=zsub(idxtarget_final);
        if plottest
            plot(vsub,zsub,'k');
            plot(vlab,zlabtemp,'r.','markersize',15);
            drawnow;
        end
%         pause;
        zlab(i,j)=zlabtemp;
    end
end

if plottest
    hold off;
end
%% merge high velocity patches
load('Vs_anomalypatches.mat');
vpatchplot=vpatch(1:2);
nlayers=length(vpatchplot);
vpatchx=vpatchplot{1}.x;
vpatchy=vpatchplot{1}.y;
vpatchgrid_hv=zeros(size(vpatchplot{1}.hv));
vpatchgrid_lv=zeros(size(vpatchplot{1}.lv));
for i=1:length(vpatchplot)
    vpatchgrid=vpatchplot{i}.hv;
    vpatchgrid_hv = vpatchgrid_hv +vpatchgrid;
    vpatchgrid=vpatchplot{i}.lv;
    vpatchgrid_lv = vpatchgrid_lv +vpatchgrid;
end   
%
vpatchgrid_hv(vpatchgrid_hv<1) = nan;
vpatchgrid_hv(vpatchgrid_hv>0) = 1;
vpatchgrid_lv(vpatchgrid_lv<1) = nan;
vpatchgrid_lv(vpatchgrid_lv>0) = 1;
vpatchgrid_regional=amask;
%%
vpatchgrid_hv_temp=vpatchgrid_hv;
vpatchgrid_hv_temp(isnan(vpatchgrid_hv))=0;
patchedge=1.0*edge(vpatchgrid_hv_temp);
patchedge(patchedge<1)=nan;
%% Plot extracted LAB depth
plotoutline=1;
crange=[60 180];
cinterval=5;


smoothx=5;
smoothy=smoothx;
smooth_surface=1;
zlab_smoothed=smooth2a(zlab,smoothx,smoothy);

figure('position',[200,400,400,400]);

hold on;

if smooth_surface
    h2=image(x,y,zlab_smoothed); 
    amask2=amask;
else
    h2=image(x,y,zlab); 
    amask2=amask;
    
end
h2.CDataMapping='scaled';
amask2(isnan(zlab))=0;
h2.AlphaData=amask2;

colormap(flip(parula(round(range(crange/cinterval)))))
hcbar=colorbar('eastoutside');
set(gca,'CLim',crange);
set(hcbar,'TickDirection','out','Ydir','reverse');
hcbar.Label.String='Extracted LAB depth (km)';

axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
daspect([1 cosd(mean(maparea.lat)) 1]);

for sb=1:length(state)
    plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[.5 .5 .5],'LineWidth',.75);
end
for sb=1:length(canada_provs)
    plot(canada_provs(sb).Lon, canada_provs(sb).Lat,'color',0.4*[1 1 1],'LineWidth',.25);
end
if plotoutline
    flist=dir('outline_*.txt');
    for k =1:length(flist)
        infile=flist(k).name;
        outline=load(infile);
        if contains(infile,'basin') || contains(infile,'rift') || ...
                contains(infile,'mcr') || contains(infile,'rr') || contains(infile,'dome')
            %close polygon
            outline(end+1,:)=outline(1,:);
        end
        disp(infile)
        if contains(infile,'ibasin') || contains(infile,'mbasin') %|| contains(infile,'dome')
            plot(outline(:,1),outline(:,2),'k','linewidth',1.75);%,'color',0.*[1 1 1]);
%         elseif contains(infile,'fcbasin') || contains(infile,'afbasin') %|| contains(infile,'dome')
%             plot(outline(:,1),outline(:,2),'k','linewidth',.75);%,'color',0.*[1 1 1]);
        elseif contains(infile,'dome')
            plot(outline(:,1),outline(:,2),'k--','linewidth',.75);%,'color',0.*[1 1 1]);
        elseif contains(infile,'rift') || contains(infile,'mcr') || contains(infile,'rr') 
            plot(outline(:,1),outline(:,2),'k:','linewidth',1.5);%,'color',0.*[1 1 1]);
%         elseif contains(infile,'border') 
%             plot(outline(:,1),outline(:,2),'k-','linewidth',1.5,'color',"#EDB120");
%         elseif contains(infile,'ndline')
%             plot(outline(:,1),outline(:,2),'k--','linewidth',1.5,'color',"#EDB120");
        end
    end
end
title('LAB depth')

hold off;
axis on;
box on;
set(gca,'TickDir','out','fontsize',14);
set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
eval(['print -depsc -r300 -painters extracted_LAB_depth.eps']);

%%
figure('position',[200,400,400,400]);

hold on;

%plot HV patch edges (merged)
amask3=amask;
p1=pcolor(vpatchx,vpatchy,patchedge);shading flat;
p1.AlphaData = amask3;
colormap(flag)
hcbar=colorbar('eastoutside');
set(gca,'CLim',crange);
set(hcbar,'TickDirection','out','Ydir','reverse');
hcbar.Label.String='Extracted LAB depth (km)';
axis([maparea.lon(1) maparea.lon(2) maparea.lat(1) maparea.lat(2)]);
daspect([1 cosd(mean(maparea.lat)) 1]);

for sb=1:length(state)
    plot(state(sb).polygon(:,1), state(sb).polygon(:,2),'color',[.5 .5 .5],'LineWidth',.75);
end
for sb=1:length(canada_provs)
    plot(canada_provs(sb).Lon, canada_provs(sb).Lat,'color',0.4*[1 1 1],'LineWidth',.25);
end
if plotoutline
    flist=dir('outline_*.txt');
    for k =1:length(flist)
        infile=flist(k).name;
        outline=load(infile);
        if contains(infile,'basin') || contains(infile,'rift') || ...
                contains(infile,'mcr') || contains(infile,'rr') || contains(infile,'dome')
            %close polygon
            outline(end+1,:)=outline(1,:);
        end
        disp(infile)
        if contains(infile,'ibasin') || contains(infile,'mbasin') %|| contains(infile,'dome')
            plot(outline(:,1),outline(:,2),'k','linewidth',1.75);%,'color',0.*[1 1 1]);
%         elseif contains(infile,'fcbasin') || contains(infile,'afbasin') %|| contains(infile,'dome')
%             plot(outline(:,1),outline(:,2),'k','linewidth',.75);%,'color',0.*[1 1 1]);
        elseif contains(infile,'dome')
            plot(outline(:,1),outline(:,2),'k--','linewidth',.75);%,'color',0.*[1 1 1]);
        elseif contains(infile,'rift') || contains(infile,'mcr') || contains(infile,'rr') 
            plot(outline(:,1),outline(:,2),'k:','linewidth',1.5);%,'color',0.*[1 1 1]);
%         elseif contains(infile,'border') 
%             plot(outline(:,1),outline(:,2),'k-','linewidth',1.5,'color',"#EDB120");
%         elseif contains(infile,'ndline')
%             plot(outline(:,1),outline(:,2),'k--','linewidth',1.5,'color',"#EDB120");
        end
    end
end
title('HV merged outlines')

hold off;
axis on;
box on;
set(gca,'TickDir','out','fontsize',14);
set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
eval(['print -depsc -r300 -painters extracted_LAB_depth_HVedges.eps']);
%%
points=[83, 156;
        83, 81;
        104,82;
        73,170];
        %72,186];
colorlist={'k','b','c','r'};
figure('Position',[400 400 1000 350]); 
figlabel={'a. ','b. ','c. ','d. ','e. ','f. ','g. ','h. ','i. ','j. '};

for k=1:size(points,1)
    i=points(k,1);
    j=points(k,2);
    vz=squeeze(vplot(i,j,:));
    if isnan(sum(vz))
        continue;
    end
    vzinterp=interp1(z,vz,zinterp);
    vsub=vzinterp(depthidxvlab);
    [vdrop,firstneg,firstpos,idxtarget,vtarget]=find_firstdrop(vsub);
    if isnan(vdrop)
        
        disp('no 1st drop')
        continue;
    else
        if firstpos < length(vsub)-1 &&  vdrop < vdrop_min && max(vsub(firstpos + 1:end)) >= vmin
            %search second drop
            [vdrop2,firstneg2,firstpos2,idxtarget2,vtarget2]=find_firstdrop(vsub(firstpos + 1:end));
            if isnan(vdrop2)
                disp('no second drop')
                idxtarget_final=idxtarget;
                idxtop=firstneg;
                idxbot=idxtarget_final;
            else      
                idxtarget_final=firstpos + idxtarget2;
                idxtop=firstpos + firstneg2;
                idxbot=firstpos + firstpos2;
            end
        else
            idxtarget_final=idxtarget;
            idxtop=firstneg;
            idxbot=firstpos;
        end
    end

    vtop=vsub(idxtop);
    vbot=vsub(idxbot);
    ztop=zsub(idxtop);
    zbot=zsub(idxbot);
    vlab=vsub(idxtarget_final);
    zlabtemp=zsub(idxtarget_final);

    subplot(1,4,k)
    hold on;
    
    
    plot(vsub,zsub,'DisplayName',['Case: ',num2str(k)],'linewidth',1,'color','b');
    h=plot(vlab,zlabtemp,'.','markersize',25,'color','b','DisplayName',['LAB']);
    yline(ztop,'k','linewidth',2);
    yline(zbot,'k--','linewidth',2);
    if k==1
        xbarloc=4.8;
        yloc=mean([ztop,zbot]);
    else
        xbarloc=4.2;
        yloc=ztop+range([ztop,zbot])*0.2;
    end
    plot([xbarloc,xbarloc],[ztop,zbot],'r:','linewidth',1.5)

    text(xbarloc+0.08,yloc,[num2str(vtop-vbot,2),' km/s'],'fontsize',12,'color','r');


    title([figlabel{k},'Case ',num2str(k)],'fontsize',14)
    
    hold off;
    grid on;
    axis on;
    box on;
    xlabel('Vs (km/s)')
    ylabel('Depth (km)')
    xlim([4,5.4])
    ylim(myzlimlab);

    legend(h,'location','southeast');
    set(gca,'fontsize',12,'ydir','reverse','Tickdir','out','TitleHorizontalAlignment','left')
end
orient landscape;
set(gcf,'renderer','Painters','PaperPositionMode','auto'); 
eval(['print -depsc -r300 -painters extracted_LAB_exampleprofiles.eps']);
