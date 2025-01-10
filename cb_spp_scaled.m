clear

%% replace by importing your data
% dimensions for the imported data "obs": year, month, variables, stations
yr = 1961:2010;
[ obs,yr_obs,spt,vrn ] = read_myobs_data( 'D:\CG\For yours\QZGY\solve\spp_extracted_station_data_comb.xlsx',1 );
[LY,LM,nv,ns] = size(obs);
[~,ia,~] = intersect(yr_obs,yr);

%% calculate scale and merge the data from different stations
% obs_c: the scaled data
% obs_sm: station-merged scaled data 
IS = 3;
isclzsc = 1;
obs_c = zeros(LY,LM,nv,ns)*nan;
for I = 1:ns
    for J = 1:nv
        if I==ns
            a1 = obs(7:end,1:7,J,I);
            a1s = sum(a1,2);
            a2 = obs(7:end,8:end,J,I);
            a2s = a2;
            for i = 1:size(a2,2)
                inn = ~isnan(a2(:,i));
                aX = a1s(inn);
                aY = a2(inn,i);
                b = polyfit(aX,aY,1);
                ye = b(1)*aX(~inn)+b(2);
                a2s(~inn,i) = ye;
            end
            qmls = [a1 a2s];
            obs(7:end,:,J,I) = qmls;
        end
        for m = 1:LM
            zs = mean(obs(ia,m,J,:),4,'omitnan');
            Xs = obs(ia,m,J,I);
            X = obs(:,m,J,I);
            if isclzsc==1
                zsm = mean(zs,'omitnan');
                zsv = var(zs,'omitnan');
                strhd = 'scaled';
            elseif isclzsc==2
                zsm = 0;
                zsv = 1;
                strhd = 'zscored';
            end
            Xm = mean(Xs,'omitnan');
            Xv = var(Xs,'omitnan');
            b = sqrt(Xv/zsv);
            a = b*zsm-Xm;
            X1 = (X+a)/b;
            obs_c(:,m,J,I) = X1;
        end
    end
end
obs_sm = mean(obs_c,4,'omitnan');

cut_yr = 1959;
ICY = find(yr_obs==cut_yr,1,'first');
obs_sm(1:ICY,:,:) = nan;

%% write out
filename = ['solve\',strhd,'_Õ¾µãÊý¾Ý_½Ø¶ÏÊ±¼ä',num2str(cut_yr-1900+1),'.xlsx'];
HF = size(dir(filename));
if HF(1)
    delete(filename)
end
wd = 1:12;
spt = '´ó²ñµ©+µÂÁî¹þ+²è¿¨+¸Õ²ì+¸ñ¶ûÄ¾+ÅµÄ¾ºé+¶¼À¼+Âê¶à+ÐËº£+ÇúÂéÀ³';
ix = 1;
for j = 1:nv
    obsm = mean(obs_sm(:,:,j),1,'omitnan');
    wt_cal = [{spt},{[]},vrn(j)];
    xlswrite(filename,wt_cal,1,['A',num2str(ix)]);
    xlswrite(filename,yr_obs',1,['A',num2str(ix+2)]);
    xlswrite(filename,{'AVERAGE'},1,['A',num2str(ix+2+LY)]);
    xlswrite(filename,{[num2str(LY),'\m']},1,['A',num2str(ix+1)]);
    xlswrite(filename,wd,1,['B',num2str(ix+1)]);
    xlswrite(filename,[obs_sm(:,:,j);obsm],1,['B',num2str(ix+2)]);
    ix = ix+LY+3;
end