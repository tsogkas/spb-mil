function res = keep_points(inp,keep,mode);
res = inp; if ~exist('mode','var'), mode =0; end
fieldnm = fieldnames(inp);
nfields =  length(fieldnm);
for k =1:nfields,
    eval(sprintf('is_cell = iscell(inp.%s);',fieldnm{k}));
    if k==nfields,  is_descr = strcmp(fieldnm{k},'descriptor'); else is_descr = 0; end
    if ~is_descr,
        if is_cell,
            try
                eval(sprintf('res.%s = inp.%s(keep);',fieldnm{k},fieldnm{k}));
            catch
                eval(sprintf('res.%s = [];',fieldnm{k},fieldnm{k}));
            end
        else
            try,
                if (~isempty(findstr('patch',fieldnm{k}))|findstr('_d_',fieldnm{k}))&(mode==0),
                    eval(sprintf('res.%s = inp.%s(keep,:,:);',fieldnm{k},fieldnm{k}));
                elseif mode==1,
                    eval(sprintf('res.%s = inp.%s(keep,:);',fieldnm{k},fieldnm{k}));
                else
                    eval(sprintf('res.%s = inp.%s(:,keep);',fieldnm{k},fieldnm{k}));
                end
            catch,
                try,
                    eval(sprintf('res.%s = inp.%s(keep);',fieldnm{k},fieldnm{k}));
                catch
                    eval(sprintf('res.%s = inp.%s;',fieldnm{k},fieldnm{k}));
                end
            end
        end
    else
        ndescr = length(inp.descriptor);
        for k=1:ndescr,
            res.descriptor{k} =  keep_points(inp.descriptor{k},keep,1);
            res.descriptor{k}.slide = keep_points(inp.descriptor{k}.slide,keep,1);
        end
    end
end

