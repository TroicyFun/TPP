function [ var1,yrsp1,varnew1 ] = tool_yr_uniform( var0,yrsp0,varnew,yrnew )
%[ var1,yrsp1,varnew1 ] = tool_yr_uniform( var0,yrsp0,varnew,yrnew )
%   Uniform the old and new variables into a same year span.

A0 = size(var0);
An = size(varnew);

if isequal(yrsp0,yrnew)
    yrsp1 = yrsp0;
    var1 = var0;
    varnew1 = varnew;
else
    yrmin = min(min(yrsp0),min(yrnew));
    yrmax = max(max(yrsp0),max(yrnew));
    yrsp1 = union(yrsp0,yrnew);
    if length(yrsp1)>length(yrsp0)
        A0(1) = length(yrsp1);
        var1 = zeros(A0)*nan;
        for k = 1:length(yrsp0)
            IY0(k) = find(yrsp1==yrsp0(k));
        end
        var1(IY0,:,:,:,:,:) = var0;
    else
        var1 = var0;
        yrsp1 = yrsp0;
    end
    for k = 1:length(yrnew)
        IYN(k) = find(yrsp1==yrnew(k));
    end
    An(1) = length(yrsp1);
    varnew1 = zeros(An)*nan;
    varnew1(IYN,:,:,:,:,:) = varnew;
end

end

