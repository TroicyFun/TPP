function [ X_flt,B,A ] = calc_filter_anal( X,yr,fmeth,fmval,fmpa,ispar )
%[ X_flt,b,a,Ff ] = calc_filter_anal( X,yr,fmeth,fmval,fmpa,ispar )
%   Gridded data filtering analysis.
% [Input]:
% fmeth: 's','l','h','b' = 'smooth','lowpass','highpass','bandpass'
%  fmval  Method    fmpa(opt)
%  1) 'smooth'
%   1  'moving'     [span]
%   2  'lowess'     [span]
%   3  'loess'      [span]
%   4  'rlowess'    [span]
%   5  'rloess'     [span]
%   6  'sgolay'     [span (order)]
%  2) 'lowpass'
%   1  'butterworth'     [Fc(=1/Tc) (Fs, default=1) (order, default=5)]
%   2  'fir'             [Fc(=1/Tc) (Fs, default=1) (order, default=1.5*Tc)]
%   3  'lanczos'         [Fc(=1/Tc) (Fs, default=1) ((order-1)/2 <order must be odd>, default=100)]
%   4  'sgolay'          [span (order, default=1)]
%   5  'median'          [span]
%  3) 'highpass'
%   1  'butterworth'     [Fc(=1/Tc) (Fs, default=1) (order, default=5)]
%   2  'fir'             [Fc(=1/Tc) (Fs, default=1) (order, default=1.5*Tc)]
%   3  'lanczos'         [Fc(=1/Tc) (Fs, default=1) ((order-1)/2 <order must be odd>, default=100)]
%   4  'diff'            []
%  4) 'bandpass'
%   1  'butterworth'     [Fca(=1/Tca) Fcb(=1/Tcb) (Fs, default=1) (order, default=5)]
%   2  'fir'             [Fca(=1/Tca) Fcb(=1/Tcb) (Fs, default=1) (order, default=1.5*Tcb)]
%   3  'lanczos'         [Fca(=1/Tca) Fcb(=1/Tcb) (Fs, default=1) ((order-1)/2 <order must be odd>, default=100)]
% [Output]:
%  X_flt: Filtered series
%  B: window para1 (B & A only for butterworth, fir, lanczos)
%  A: window para2 (For lanczos: Ff, because A≡1)

if nargin<6
    ispar = 0;
end

[LY,ny,nx,nz,nt,nv] = size(X);
nal = ny*nx*nz*nt*nv;
Xr = reshape(X,LY,nal);
Xr_flt = nan(LY,nal);
B = [];
A = [];
if ~ispar
    for i = 1:nal
        Xr_c = Xr(:,i);
        IC = find(~isnan(Xr_c));
        Xr_c = Xr_c(IC);
        yr_c = yr(IC);
        Xr_flt_c = nan(LY,1);
        a = [];
        b = [];
        yy = [];
        if ~isempty(Xr_c)
            switch fmeth
                case {'s','smooth'}  % smooth
                    smmth = {'moving','lowess','loess','rlowess','rloess','sgolay'};
                    if ischar(fmval)
                        fmvalc = fmval;
                    else
                        fmvalc = smmth{fmval};
                    end
                    if fmval==6 && length(fmpa)>1
                        yy = smooth(yr_c,Xr_c,fmpa(1),fmvalc,fmpa(2));
                    else
                        yy = smooth(yr_c,Xr_c,fmpa,fmvalc);
                    end
                case {'l','lowpass'}  % lowpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth low
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                            [b, a] = butter(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'low');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir low
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                            [b, a] = fir1(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'low');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos low
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        end
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(2),fmpa_a(1),fmpa_a(3),'low');
                    elseif (isnumeric(fmval) && fmval==4) || strcmpi(fmval,'sgolay')  % sgolay low
                        if length(fmpa)<2
                            fmpa(2) = 1;
                        end
                        yy = sgolayfilt(Xr_c, fmpa(2), fmpa(1));
                    elseif (isnumeric(fmval) && fmval==5) || strcmpi(fmval,'median')  % median low
                        yy = medfilt1(Xr_c,fmpa,'omitnan');
                    end
                case {'h','highpass'}  % highpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth high
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                            [b, a] = butter(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'high');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir high
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                            [b, a] = fir1(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'high');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos high
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        end
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(2),fmpa_a(1),fmpa_a(3),'high');
                    elseif (isnumeric(fmval) && fmval==4) || strcmpi(fmval,'diff')  % diff high
                        yy = [nan;diff(Xr_c)];
                    end
                case {'b','bandpass'}  % bandpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth bandpass
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                            [b, a] = butter(fmpa_a(4), [fmpa_a(1) fmpa_a(2)]/(fmpa_a(3)/2), 'bandpass');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir bandpass
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                            [b, a] = fir1(fmpa_a(4), [fmpa_a(1) fmpa_a(2)]/(fmpa_a(3)/2), 'bandpass');
                        end
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos bandpass
                        if i==1
                            fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        end
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(3),fmpa_a(1),fmpa_a(4),'band',fmpa_a(2));
                    end
            end
            Xr_flt_c(IC) = yy;
        end
        if ~isempty(b)
            Br(:,i) = b;
        end
        if ~isempty(a)
            Ar(:,i) = a;
        end
        Xr_flt(:,i) = Xr_flt_c;
    end
else
    parfor i = 1:nal
        Xr_c = Xr(:,i);
        IC = find(~isnan(Xr_c));
        Xr_c = Xr_c(IC);
        yr_c = yr(IC);
        Xr_flt_c = nan(LY,1);
        a = [];
        b = [];
        yy = [];
        if ~isempty(Xr_c)
            switch fmeth
                case {'s','smooth'}  % smooth
                    smmth = {'moving','lowess','loess','rlowess','rloess','sgolay'};
                    if ischar(fmval)
                        fmvalc = fmval;
                    else
                        fmvalc = smmth{fmval};
                    end
                    if fmval==6 && length(fmpa)>1
                        yy = smooth(yr_c,Xr_c,fmpa(1),fmvalc,fmpa(2));
                    else
                        yy = smooth(yr_c,Xr_c,fmpa,fmvalc);
                    end
                case {'l','lowpass'}  % lowpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth low
                        fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                        [b, a] = butter(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'low');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir low
                        fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                        [b, a] = fir1(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'low');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos low
                        fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(2),fmpa_a(1),fmpa_a(3),'low');
                    elseif (isnumeric(fmval) && fmval==4) || strcmpi(fmval,'sgolay')  % sgolay low
                        fmpa_a = fmpa;
                        if length(fmpa_a)<2
                            fmpa_a(2) = 1;
                        end
                        yy = sgolayfilt(Xr_c, fmpa_a(2), fmpa_a(1));
                    elseif (isnumeric(fmval) && fmval==5) || strcmpi(fmval,'median')  % median low
                        yy = medfilt1(Xr_c,fmpa,'omitnan');
                    end
                case {'h','highpass'}  % highpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth high
                        fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                        [b, a] = butter(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'high');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir high
                        fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                        [b, a] = fir1(fmpa_a(3), fmpa_a(1)/(fmpa_a(2)/2), 'high');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos high
                        fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(2),fmpa_a(1),fmpa_a(3),'high');
                    elseif (isnumeric(fmval) && fmval==4) || strcmpi(fmval,'diff')  % diff high
                        yy = [nan;diff(Xr_c)];
                    end
                case {'b','bandpass'}  % bandpass
                    if (isnumeric(fmval) && fmval==1) || strcmpi(fmval,'butterworth')  % butterworth bandpass
                        fmpa_a = check_fmpa(fmpa,fmeth,'butterworth');
                        [b, a] = butter(fmpa_a(4), [fmpa_a(1) fmpa_a(2)]/(fmpa_a(3)/2), 'bandpass');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==2) || strcmpi(fmval,'fir')  % fir bandpass
                        fmpa_a = check_fmpa(fmpa,fmeth,'fir');
                        [b, a] = fir1(fmpa_a(4), [fmpa_a(1) fmpa_a(2)]/(fmpa_a(3)/2), 'bandpass');
                        yy = filtfilt(b, a, Xr_c); % 使用零相位延迟滤波
                    elseif (isnumeric(fmval) && fmval==3) || strcmpi(fmval,'lanczos')  % lanczos bandpass
                        fmpa_a = check_fmpa(fmpa,fmeth,'lanczos');
                        [yy,~,b,~,a] = lanczos_filter(Xr_c,1/fmpa_a(3),fmpa_a(1),fmpa_a(4),'band',fmpa_a(2));
                    end
            end
            Xr_flt_c(IC) = yy;
        end
        if ~isempty(b)
            Br(:,i) = b;
        end
        if ~isempty(a)
            Ar(:,i) = a;
        end
        Xr_flt(:,i) = Xr_flt_c;
    end
end
X_flt = reshape(Xr_flt,LY,ny,nx,nz,nt,nv);
if exist('Br','var')
    LB = size(Br,1);
    B = reshape(Br,LB,ny,nx,nz,nt,nv);
end
if exist('Ar','var')
    LA = size(Ar,1);
    A = reshape(Ar,LA,ny,nx,nz,nt,nv);
end
end

function fmpa_out = check_fmpa(fmpa_in,fmeth,fmvalc)
fmpa_out = fmpa_in;
lf = length(fmpa_out);
if strcmpi(fmeth(1),'b')
    mxlf = 4;
else
    mxlf = 3;
end
if lf==mxlf-2
    fmpa_out(lf+1) = 1;
    lf = length(fmpa_out);
end
if lf==mxlf-1
if strcmpi(fmvalc,'butterworth')
    fmpa_out(lf+1) = 5;
elseif strcmpi(fmvalc,'fir')
    fmpa_out(lf+1) = floor(1.5/fmpa_in(lf-1));
elseif strcmpi(fmvalc,'lanczos')
    fmpa_out(lf+1) = 100;
end
lf = length(fmpa_out);
end
if lf==mxlf
    return
else
    error('Number of input fmpa is wrong!')
end

end