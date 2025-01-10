function [ sr_flt ] = auto_filter( sr_data,yr_sr,fmval,fmpa )
%[ sr_flt ] = auto_filter( sr_data,yr_sr,fmval,fmpa )
%   Do filtering automatically.

for I = 1:size(sr_data,2)
    sr_c = sr_data(:,I);
    IC = find(~isnan(sr_c));
    sr_c = sr_c(IC);
    yr_sr_c = yr_sr(IC);
    sr_flt(:,I) = zeros(length(yr_sr),1)*nan;
    if ~isempty(sr_c)
        switch fmval
            case {1,2,3,4,5,6}  % smooth
                smmth = {'moving','lowess','loess','sgolay','rlowess','rloess'};
                if fmval~=4 || length(fmpa)==1
%                     if fmval==2 ||fmval==3
%                         fmpa = fmpa/length(yr_sr_c);
%                     end
                    yy = smooth(yr_sr_c,sr_c,fmpa,smmth{fmval});
                else
                    yy = smooth(yr_sr_c,sr_c,fmpa(1),smmth{fmval},fmpa(2));
                end
            case 7  % median
                yy = medfilt1(sr_c,fmpa,'omitnan');
            case 8  % Gaussain smooth
                yy = Gauss_filter( sr_c,fmpa );
            case 9  % highpass:difference
                yy = diff(sr_c);
                yy(end+1) = nan;
            case 10  % bandpass
                if length(fmpa)==1
                    yy = bp_filter( sr_c,1/fmpa );
                else
                    yy = bp_filter( sr_c,1/fmpa(1),fmpa(2) );
                end
            case 11  % FBR
                yy = FBR_filter( sr_c,fmpa,[0 0],3 );
            case 12  % Lanczos low-pass
                yy = lanczosfilter(sr_c,[],1/fmpa,[],'low');  
            case 13  % Lanczos high-pass
                yy = lanczosfilter(sr_c,[],1/fmpa,[],'high');  
        end
        sr_flt(IC,I) = yy;
    end
end

end

