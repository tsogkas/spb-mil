function line_features = PSLz3b_embed_lines(line_parsing,thresh_line_ener);
% line_features = PSLz3b_embed_lines(line_parsing,thresh_line_ener)
%
% For each of the straight line tokens, get its pose descriptor and the
% average strength of the differential operator along its corresponding curve.
%
% Iasonas Kokkinos <jkokkin@stat.ucla.edu>
% 10/10/2007

lines = line_parsing.lines;
merit_geom = line_parsing.merit_geom;
line_count = 0;
lines_kept = [];

for cnt=1:length(lines), 
    pts = lines{cnt};
    idxs = pts(3,:);

    
	ener_wt = max(pts(4,:));
    if (ener_wt>thresh_line_ener),
        line_count = line_count+1;

        lines_kept=  [lines_kept,cnt];
        
        ener(line_count) =  ener_wt;        
        c_m(line_count)  = mean(pts(1,:));
        c_n(line_count)  = mean(pts(2,:));
        scales(line_count) = sqrt(sum((pts(1:2,1) - pts(1:2,end)).^2))/2;            
        orientations(line_count)=  get_orientation(pts,c_m(line_count),c_n(line_count));

        scales_clipped =  clip(pts(5,:),.01,inf);
        ratios(line_count) = exp(mean(log(scales_clipped)))/max(scales(line_count),1);        
    end
end

fields_wt = fieldnames(line_parsing);
for k=1:length(fields_wt),
    eval(sprintf('%s = %s(lines_kept);',fields_wt{k},fields_wt{k}));
end
if isempty(lines_kept),
    orientations =[]; ratios = []; scales =[]; c_m =[]; c_n =[]; ener = [];
end

fields_wt = [fields_wt(:)',{'orientations','ratios','scales','c_m','c_n','ener'}];
compress_structure;    
line_features = structure;

function    theta = get_orientation(pts,mean_c_x,mean_c_y);
%% form least squares approximation 
%% to the points belonging to the line
%% and use the parameters to determine the angle's orientation
c_x = pts(1,:);
c_y = pts(2,:);
c_x_2= pow_2(c_x);
c_y_2 = pow_2(c_y);

mean_c_x_2 = mean(c_x_2);
mean_c_y_2 = mean(c_y_2);
mean_c_xy = mean(c_x.*c_y);

pow_2_mcx = mean_c_x*mean_c_x;
pow_2_mcy = mean_c_y*mean_c_y;

design(1,1) = mean_c_x_2 - pow_2_mcx;
design(1,2) = mean_c_xy - mean_c_x*mean_c_y;
design(2,1) = design(1,2);
design(2,2) = mean_c_y_2 - pow_2_mcy;

[eigv,t] = eig(design);
theta = atan2(-eigv(2,1),eigv(1,1));


       