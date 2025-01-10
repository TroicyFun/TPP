clear
sdpath = 'D:\CG\TreeRingLab\other_data\';
mstr = 'P7C10';

%% replace by importing your data
% dimensions for the imported data "obs_aqc": year, month, stations, variables
[ chr,yr_chr,~ ] = read_chronology( '青藏高原_重建降水数据.xlsx' );
inn = find(isnan(chr));
chr(inn) = (chr(inn-1)+chr(inn+1))/2;

st_cal = tool_get_st_inp( '大柴旦 德令哈 刚察 格尔木 诺木洪 都兰 茶卡 兴海 曲麻莱 玛多' );
vr_cal = tool_get_st_inp( 'tmean pr rhmean' );
[ cal_num1,cal_name1,vr_cal1,ise1 ] = tool_get_num_name( st_cal,vr_cal,sdpath,0 );
[ cal_var,ise1 ] = tool_extract_station_data( cal_num1,vr_cal1,ise1,sdpath );
[LY,ND,ns,nv] = size(cal_var);
year_span = 1951:2020;
obs_aqc = zeros(LY,12,ns,nv)*nan;
for i = 1:nv
    for j = 1:ns
        obs_aqc(:,:,j,i) = quality_control( cal_var(:,:,j,i),year_span,'M',vr_cal1{i} );
    end
end

%% calculate and plot
mstr0 = 'P1C12';
obs_adj0 = zeros(LY,12*2,ns,nv)*nan;
obs_adj0(2:end,1:12,:,:) = obs_aqc(1:end-1,:,:,:);
obs_adj0(:,13:24,:,:) = obs_aqc;
m_flg0 = repmat({'J','F','M','A','M','J','J','A','S','O','N','D'},1,2);
flg_mth0 = repmat({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'},1,2);
flg_yr0 = reshape(repmat({'P';'C'},1,12)',1,24);
flg_im0 = repmat(1:12,1,2);
[ mth0 ] = get_mth( mstr0 );
[ mth1 ] = get_mth( mstr );
[~,IMTH,~] = intersect(mth0,mth1);
m_flg1 = m_flg0(IMTH);
flg_im1 = flg_im0(IMTH);
obs_adj1 = obs_adj0(:,IMTH,:,:);
flg_mth1 = flg_mth0(IMTH);
flg_yr1 = flg_yr0(IMTH);
LMI = length(IMTH);

[yr_int,ia,ib] = intersect(yr_chr,year_span);
LY = length(yr_int);

% Corr.
maxj = 12;
R = zeros(nv,LMI+1,ns)*nan;
P = R;
N = R;
for I = 1:ns
    chr_c = -chr(ia);
    for V = 1:nv
        obs_c = obs_adj1(ib,:,I,V);
        ik = 0;
        for i = 1:LMI
            o_to = zeros(LY,1);
            for j = i:min(LMI,i+maxj-1)
                ik = ik+1;
                o_to = o_to+obs_c(:,j);
                o_tom = o_to./(j-i+1);
                inn = ~isnan(o_to);
                o_to1 = o_tom(inn);
                chr_c1 = chr_c(inn);
                [ r0,p0 ] = corrcoef( [chr_c1,o_to1] );
                r_et(ik,:,I,V) = [r0(2) p0(2) sum(inn)];
                obs_et(:,ik,I,V) = o_tom;
                loc_et(ik,:,I,V) = [i j];
            end
        end
        IRD = find(loc_et(:,1,I,V)==loc_et(:,2,I,V));
        R(V,1:LMI,I) = r_et(IRD,1,I,V);
        P(V,1:LMI,I) = r_et(IRD,2,I,V);
        N(V,1:LMI,I) = r_et(IRD,3,I,V);
    end
end
for I = 1:ns
    [rm1,IRM1] = min(squeeze(r_et(:,2,I,:)));
    [rm2,IRM2] = min(rm1);
    IRM = IRM1(IRM2);
    R(:,LMI+1,I) = squeeze(r_et(IRM,1,I,:));
    P(:,LMI+1,I) = squeeze(r_et(IRM,2,I,:));
    N(:,LMI+1,I) = squeeze(r_et(IRM,3,I,:));
    nmax(I) = squeeze(r_et(IRM,3,I,IRM2));
    lc1 = loc_et(IRM,1,I,IRM2);
    lc2 = loc_et(IRM,2,I,IRM2);
    loc_rmax{I,1} = [flg_yr1{lc1},'_{\fontsize{5.5}',flg_mth1{lc1},'}-',flg_yr1{lc2},'_{\fontsize{5.5}',flg_mth1{lc2},'}'];
end


% Plots
stnm_e = {'DaChaiDan','DeLingHa','GangCha','GeErMu','NuoMuHong','DuLan',...
          'ChaKa','XingHai','QuMaLai','MaDuo'};
ftn = 'Arial';
fts = 7.5;
NR = 4;
NC = 3;
width = 220;
left = 50+(0:NC-1)*(width+15);
ttw = left(end)+width+17;
height = 105;
bottom = 52+(NR-1:-1:0)*(height+50);
tth = bottom(1)+height+22;

fignm = 'CORR_TPP_STT';
figure('Position',[50 50 ttw tth],'Name',fignm,'Color','w');
for I = 1:ns
subplot('Position',[left(mod(I-1,NC)+1)/ttw bottom(ceil(I/NC))/tth width/ttw height/tth]);
xt = [mth1 mth1(end)+1];
xtl = m_flg1;
xtl(end+1) = {''};
xl = [min(xt)-.8 max(xt)+.8];
yl = [-.8 .4];
xsg1 = .5;
xsg2 = mth1(end)+.5;
r_sg = rinv(.05,nmax(I)-2);
% for k = 1:nv
%     bscolor(k,:) = [(k-1)/V (k-1)/V 1-(k-1)/V];
% end
% bscolor(1,:) = [0 0 0];
bscolor = [255 131 65; 2 136 255; 132 132 132]./255;
hb = bar(xt,R(:,:,I)','DisplayName','mth','barwidth',.68);hold on
for k=1:V
    hb(k).FaceColor = bscolor(k,:);
    hb(k).EdgeColor = bscolor(k,:);
end
hbs = plot(xl,[0 0],'k-');hold on
hb(4) = plot(xl,[r_sg r_sg],'r--');hold on
hb(4) = plot(xl,-[r_sg r_sg],'r--');hold on
hsg = plot([xsg2 xsg2],yl,'color','k');hold on

axis([xl yl]);
set(gca,'xcolor','k','ycolor','k','xtick',xt,'xticklabel',xtl,'xticklabelrotation',0,'fontname',ftn,'fontsize',fts,'ytick',-1:.2:1);
text(xl(1),yl(2)+diff(yl)*.1,[char('a'+I-1),' ',stnm_e{I}],'fontname',ftn,'fontsize',fts)
text(xt(end),yl(1)-diff(yl)*.1,loc_rmax{I,1},'fontname',ftn,'fontsize',fts,'Rotation',67.5,'HorizontalAlignment','right');

im0 = -5;
im1 = 0;
im2 = 10;
ymk_ps = yl(1)+diff(yl)*.08;
% hymk = plot([xl(1) xl(1)],[ymk_ps-diff(yl)*ymk_lt/2 ymk_ps+diff(yl)*ymk_lt/2],'k-', ...
%       [im1+.5 im1+.5],[ymk_ps-diff(yl)*ymk_lt/2 ymk_ps+diff(yl)*ymk_lt/2],'k-', ...
%       [LMI+.5 LMI+.5],[ymk_ps-diff(yl)*ymk_lt/2 ymk_ps+diff(yl)*ymk_lt/2],'k-','clipping','off');hold on
hxvc = plot([im0 im1],[ymk_ps ymk_ps],'k-',[im1+1 im2],[ymk_ps ymk_ps],'k-','clipping','off');hold on
text(mean([im0 im1]),ymk_ps+diff(yl)*.07,'Previous year','fontname',ftn,'fontsize',fts-.5,'HorizontalAlignment','center');
text(mean([im1+1 im2]),ymk_ps+diff(yl)*.07,'Current year','fontname',ftn,'fontsize',fts-.5,'HorizontalAlignment','center');
hdlt = diff(xl)*.03;
hdag = 2;
hdxm = diff(xl)*.015;
hpt(1) = patch([im0 im0+hdlt*cosd(hdag) im0+hdxm im0+hdlt*cosd(hdag)],[ymk_ps ymk_ps+hdlt*sind(hdag) ymk_ps ymk_ps-hdlt*sind(hdag)],'k','clipping','off');hold on
hpt(2) = patch([im1 im1-hdlt*cosd(hdag) im1-hdxm im1-hdlt*cosd(hdag)],[ymk_ps ymk_ps+hdlt*sind(hdag) ymk_ps ymk_ps-hdlt*sind(hdag)],'k','clipping','off');hold on
hpt(3) = patch([im1+1 im1+1+hdlt*cosd(hdag) im1+1+hdxm im1+1+hdlt*cosd(hdag)],[ymk_ps ymk_ps+hdlt*sind(hdag) ymk_ps ymk_ps-hdlt*sind(hdag)],'k','clipping','off');hold on
hpt(4) = patch([im2 im2-hdlt*cosd(hdag) im2-hdxm im2-hdlt*cosd(hdag)],[ymk_ps ymk_ps+hdlt*sind(hdag) ymk_ps ymk_ps-hdlt*sind(hdag)],'k','clipping','off');hold on

if I>=8
%     ad_ymv = yl(1)-diff(yl)*.24;
%     ad_lt = .06;
%     had = plot([xl(1) xl(1)],ad_ymv+diff(yl)*[-ad_lt ad_lt],'k-',...
%                [xsg1 xsg1],ad_ymv+diff(yl)*[-ad_lt ad_lt],'k-',...
%                [xsg2 xsg2],ad_ymv+diff(yl)*[-ad_lt ad_lt],'k-','clipping','off');hold on
%     ad_xmv = .0075;
%     vc1_x1 = xl(1)+diff(xl)*ad_xmv;
%     vc1_x2 = xsg1-diff(xl)*ad_xmv;
%     vc2_x1 = xsg1+diff(xl)*ad_xmv;
%     vc2_x2 = xsg2-diff(xl)*ad_xmv;
%     had2 = plot([vc1_x1 vc1_x2],[ad_ymv ad_ymv],'k-',...
%                 [vc2_x1 vc2_x2],[ad_ymv ad_ymv],'k-','clipping','off');hold on
%     av_lt = diff(xl)*.02;
%     av_ag = 4;
%     had3 = plot([vc1_x1+av_lt*cosd(av_ag) vc1_x1 vc1_x1+av_lt*cosd(av_ag)],[ad_ymv+av_lt*sind(av_ag) ad_ymv ad_ymv-av_lt*sind(av_ag)],'k-',...
%                 [vc1_x2-av_lt*cosd(av_ag) vc1_x2 vc1_x2-av_lt*cosd(av_ag)],[ad_ymv+av_lt*sind(av_ag) ad_ymv ad_ymv-av_lt*sind(av_ag)],'k-',...
%                 [vc2_x1+av_lt*cosd(av_ag) vc2_x1 vc2_x1+av_lt*cosd(av_ag)],[ad_ymv+av_lt*sind(av_ag) ad_ymv ad_ymv-av_lt*sind(av_ag)],'k-',...
%                 [vc2_x2-av_lt*cosd(av_ag) vc2_x2 vc2_x2-av_lt*cosd(av_ag)],[ad_ymv+av_lt*sind(av_ag) ad_ymv ad_ymv-av_lt*sind(av_ag)],'k-','clipping','off');hold on
%     text((xl(1)+xsg1)/2,ad_ymv-diff(yl)*.09,'Previous Year','fontname',ftn,'fontsize',fts,'HorizontalAlignment','center');
%     text((xsg1+xsg2)/2,ad_ymv-diff(yl)*.09,'Current Year','fontname',ftn,'fontsize',fts,'HorizontalAlignment','center');
    xlabel('Months','position',[mean(xl) yl(1)-diff(yl)*.25],'HorizontalAlignment','center');
end
if mod(I-1,NC)==0
    ylabel('Correlation Coefficient','fontname',ftn,'fontsize',fts);
else
    set(gca,'yticklabel',{});
end
if I==ns
    lgd = {'Temperature','Precipitation','Relative Humidity','\itp\rm = 0.05'};
    hl = legend(hb,lgd,'location','east','Orientation','Horizontal','box','off','fontname',ftn,'fontsize',fts,...
        'Position',[.393 .17 .571 .0236]);
end
end

set(gcf, 'PaperPositionMode', 'auto');
tool_isf_menu(4,3);
% print(gcf,'-painters','-depsc',get(gcf,'Name'));
