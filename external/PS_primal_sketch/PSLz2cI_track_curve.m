function    [passed,string_str,ener_str,scale_str] = PSLz2cI_track_curve(index_start,passed,indexes,ener,scl,lookup,succ,pred,iters);
% [passed,string_str,ener_str,scale_str] = PSLz2cI_track_curve(index_start,passed,indexes,ener,scl,lookup,succ,pred,iters)
%
% Get continuous curves by hopping from point to point as suggested by the successor/predecessor
% information. 
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

ener_temp   = zeros(1,50);
string_temp = zeros(1,50);
scale_temp  = zeros(1,50);

string_temp(1) = indexes(index_start);
ener_temp(1)   = ener(index_start);
scale_temp(1)  = scl(index_start);

for iter =1:iters,
    cnt = 1;
    cont=1;
    index = index_start;
    %% iter 1: start searching in one direction of starting point
    %% iter 2: other direction.
    while cont
        %% If we have already started tracking the curve,
        %% `next_location' should  be the next non-occupied point
        %% It has nothing to do with `iter' any longer.
        next_location = succ(index);
        if (next_location<0|(passed(next_location)==1))
            next_location = pred(index);
        end
        t=  indexes(index);
        passed(t) = 1;

        if (next_location>0)&(passed(next_location)~=1)
            cnt = cnt + 1;
            %% find next point in line
            index = lookup(next_location);
            string_temp(cnt) = next_location;
            ener_temp(cnt)  = ener(index);
            scale_temp(cnt) = scl(index);
        else
            cont = 0;
        end
    end

    if iter==1,
        string_str = [string_temp(1:cnt)];
        scale_str  = [scale_temp(1:cnt)];
        ener_str   = [ener_temp(1:cnt)];
    else
        string_str = [fliplr(string_temp(2:cnt)),string_str];
        scale_str  = [fliplr(scale_temp(2:cnt)),scale_str];
        ener_str   = [fliplr(ener_temp(2:cnt)),ener_str];
    end
end
