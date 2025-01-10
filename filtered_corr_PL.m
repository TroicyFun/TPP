clear

flt_lth = 1;
flt_meth = 2;

%% replace by importing your data
% TPP序列
[ chr,yr_chr,~ ] = read_chronology( 'D:\CG\For yours\QZGY\青藏高原_重建降水数据.xlsx',2 );
chr_nm = 'TPP';
chr_flt_c = auto_filter( chr,yr_chr,flt_meth,flt_lth );
chr_flt = auto_filter( chr_flt_c,yr_chr,flt_meth,flt_lth );
yr_union = yr_chr;
cp_union = chr_flt;
cp_nm{1,1} = chr_nm;

% 各种序列
fn = 'D:\CG\For yours\QZGY\Gen_data\树轮石笋冰芯湖泊等序列.xlsx';
max_sr_ps = 'EC';
[~, ~, ci] = xlsread(fn,1,['A1:',max_sr_ps,'1']);
ik = 0;
ix = 1;
in = 0;
for i = 1:length(ci)
    if i==ix
        ik = ik+1;
        ci1 = char(ci(i));
        IC = find(ci1=='\');
        pas(ik,:) = [str2num(ci1(1:IC-1)) str2num(ci1(IC+1:end)) ix];
        ix = ix+pas(ik,2)+2;
        for j=1:pas(ik,2)
            in = in+1;
            sr_nm(in,1) = ci(i+j);
        end
    end
end
pas(:,4) = [0;cumsum(pas(1:end-1,2))];
[~, ~, raw] = xlsread(fn,1,['A2:',max_sr_ps,num2str(max(pas(:,1))+1)]);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw);
raw(R) = {NaN};
sr = reshape([raw{:}],size(raw));

%% calculate and plot
np = size(pas,1);
nc = length(sr_nm);
NC = 5;
NR = 4;
width = 250;
left = 33+(0:NC-1)*(width+35);
ttw = left(end)+width+25;
height = 140;
bottom = 35+(NR-1:-1:0)*(height+49);
tth = bottom(1)+height+28;
nppf = NC*NR;
nfg = ceil(np./nppf);
ik = 1;
ifg = 0;
for I = 2:np
    yr_cp = flipud(sr(1:pas(I,1),pas(I,3):pas(I,3)));
    sr_cp = flipud(sr(1:pas(I,1),pas(I,3)+1:pas(I,3)+pas(I,2)));
    [LY,nsr] = size(sr_cp);

    [yr_cp_a,ia,ib] = unique(yr_cp);
    sr_cp_a = zeros(length(yr_cp_a),size(sr_cp,2));
    for i = 1:length(yr_cp_a)
        ibj = ib==i;
        sr_cp_a(i,:) = mean(sr_cp(ibj,:),1);
    end

    for j = 1:nsr
        cp_c = sr_cp_a(:,j);
        inn = ~isnan(cp_c);
        cp_c = cp_c(inn);
        yr_cp_c = yr_cp_a(inn);
        yr_cp_itp = (ceil(yr_cp_c(1)):floor(yr_cp_c(end)))';
        if length(yr_cp_itp)>flt_lth
            ik = ik+1;
            ic = ik-1;
            cp_nm_c = sr_nm{ik};
            cp_nm_c = cp_nm_c(~isspace(cp_nm_c));
            cp_itp = interp1(yr_cp_c,cp_c,yr_cp_itp);
            cp_flt_c = auto_filter( cp_itp,yr_cp_itp,flt_meth,flt_lth );
            cp_flt = auto_filter( cp_flt_c,yr_cp_itp,flt_meth,flt_lth );
            if ik==11 || ik==12
                IY0 = find(yr_cp_itp==1200,1,'first');
                yr_cp_itp = yr_cp_itp(IY0:end);
                cp_flt = cp_flt(IY0:end);
            end
            [yr_crx,ia,ib] = intersect(yr_chr,yr_cp_itp);
            r0 = corrcoef(chr_flt(ia),cp_flt(ib));
            R = r0(2);
            EDOF = length(yr_crx)/flt_lth;
            P = get_pval( R,EDOF );
            [ var1,yr_union,varnew1 ] = tool_yr_uniform( cp_union,yr_union,cp_flt,yr_cp_itp );
            cp_union = [var1,varnew1];
        
            ipar = find(cp_nm_c=='(',1,'first');
            if ~isempty(ipar)
                cp_nm_c = cp_nm_c(1:ipar-1);
            end
            ipar = find(cp_nm_c=='_');
            if ~isempty(ipar)
                cp_nm_c(ipar) = ' ';
            end

            if mod(ic,nppf)==1
                ifg = ifg+1;
                ipn = 1;
                imaname = ['combine_sr\',num2str(flt_lth),'年滤波\',num2str(flt_lth),'yr_FIG',num2str(ifg)];
                figure('Color','w','position',[1 1 ttw tth], 'Name', imaname);
                set(gcf, 'PaperPositionMode', 'auto');
                tool_isf_menu( 4,2 )
            else
                ipn = ipn+1;
            end
            subplot('Position',[left(mod(ipn-1,NC)+1)/ttw bottom(ceil(ipn/NC))/tth width/ttw height/tth]);
            if P<0.05
                p1 = '< 0.05';
            elseif P<0.1
                p1 = '< 0.1';
            else
                p1 = '> 0.1';
            end
        
            if R<0
                cp_flt = -cp_flt;
            end
            plot(yr_crx,zscore(chr_flt(ia)),yr_crx,zscore(cp_flt(ib)),'linewidth',.5);
            xl = [min(yr_crx) max(yr_crx)];
            xl = [xl(1)-diff(xl)*.025 xl(2)+diff(xl)*.025];
            yl = get(gca,'ylim');
            yl = [yl(1) yl(2)+diff(yl)*.2];
            axis([xl yl]);
            text(xl(2)-(xl(2)-xl(1))*.05,yl(2)-(yl(2)-yl(1))*.1,['L = ',num2str(flt_lth), ...
                ', EDOF = ',num2str(EDOF),', \itr\rm = ',num2str(roundn(R,-3)),', \itp\rm ',p1], ...
                'HorizontalAlignment','Right','fontsize',8);
            text(xl(1),yl(2)+(yl(2)-yl(1))*.1,[num2str(ic),' ',cp_nm_c],'fontsize',8);
            ipar = find(cp_nm_c==' ');
            if ~isempty(ipar)
                cp_nm_c(ipar) = '_';
            end
            if mod(ic,nppf)==0 || ic==nc-1
                tool_save_image( imaname,4,2 )
            end
        end
    end
    
end

%% output data
recname = ['combine_sr\',num2str(flt_lth),'年滤波\',num2str(flt_lth),'yr_TPP与序列对比.xlsx'];
HF = size(dir(recname));
if HF(1)
    delete(recname)
end
tbn = 1;
xlswrite(recname,{[num2str(length(yr_union)),'\',num2str(length(sr_nm))]},tbn,'A1');
xlswrite(recname,sr_nm',tbn,'B1');
xlswrite(recname,yr_union,tbn,'A2');
xlswrite(recname,cp_union,tbn,'B2');
