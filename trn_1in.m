function trn_1in(mf_n, epoch_n)

if nargin < 2, epoch_n = 400; end;
if nargin < 1, mf_n = 9; end;

% ====== Generate training data
THRESH=0.0;
s=load('data.mat'); u=s.u; f=s.f; data=[u,f]; clear s u f; 
    data=downsample(data,100);
    data_n=size(data,1);
    data(:,2)=-THRESH+data(:,2); idx=find(data(:,2)<0);data(idx,2)=0;data(:,2)=data(:,2)/max(data(:,2));
trn_data = data;
chk_data = [];

% ====== Training options
ss = 0.05; %ss = 0.01;
error_goal = 0;
mf_type ='gaussmf';% 'trimf';% 'gbellmf';
in_fismat = genfis1(trn_data, mf_n, mf_type);

[trn_out_fismat, trn_error, step_size, chk_out_fismat, chk_error] = ...
	anfis(trn_data, in_fismat, ...
	[epoch_n error_goal ss nan nan], [1,1,1,1], chk_data);

trn_error(find(trn_error == -1)) = NaN*find(trn_error == -1);
if ~isempty(chk_data),
	chk_error(find(chk_error == -1)) = NaN*find(chk_error == -1);
end
step_size(find(step_size == -1)) = NaN*find(step_size == -1);

in = trn_data(:, 1);
out = evalfis(in, trn_out_fismat);
    trn_RMSE = norm(out - trn_data(:,2))/sqrt(size(in, 1))

if ~isempty(chk_data),
	in = chk_data(:, 1);
	out = evalfis(in, chk_out_fismat);
	chk_RMSE = norm(out - chk_data(:,2))/sqrt(size(in, 1))
end

figure;

subplot(221);
    in_range = getfis(in_fismat, 'inrange');
    mf_type = getfis(in_fismat, 'inmftypes');
    mf_para = getfis(in_fismat, 'inmfparams');

    x = linspace(in_range(1, 1), in_range(1, 2),  size(data,1));
    mf = evalmmf(x, mf_para, mf_type);
    plot(x', mf'); axis([in_range 0 1.2]); title(['Initial MFs,n=',num2str(mf_n)]);

subplot(222);
    in_range = getfis(trn_out_fismat, 'inrange');
    mf_type = getfis(trn_out_fismat, 'inmftypes');
    mf_para = getfis(trn_out_fismat, 'inmfparams');
    x = linspace(in_range(1, 1), in_range(1, 2), size(data,1));
    mf = evalmmf(x, mf_para, mf_type);
    plot(x', mf'); axis([in_range 0 1.2]);title('Final MFs');

subplot(223);
if ~isempty(chk_data),
	data_n = size(trn_data, 1);
	new_data_n = data_n + 10*(data_n-1);
	dense_x = linspace(min(trn_data(:,1)), max(trn_data(:,1)), new_data_n)';
	dense_y = evalfis(dense_x, trn_out_fismat);
	plot(trn_data(:, 1), trn_data(:, 2), '+', ...
		chk_data(:, 1), chk_data(:, 2), 'o', dense_x, dense_y);
else
	data_n = size(trn_data, 1);
	new_data_n = data_n + 10*(data_n-1);
	dense_x = linspace(min(trn_data(:,1)), max(trn_data(:,1)), new_data_n)';
	dense_y = evalfis(dense_x, trn_out_fismat);
	plot(trn_data(:, 1), trn_data(:, 2), '+', dense_x, dense_y);
end
title('f(u) and ANFIS Outputs');

subplot(224);
    max_mf = ones(mf_n, 1)*max(mf); index = find(mf ~= max_mf);
    tmp = getfis(trn_out_fismat,'outmfparams'); tmp = tmp(:, 1:2);	% Bug due to getfis.m in FLT v2. Fixed in my getfis.m
    rule_output = tmp*[x; ones(size(x))]; 
    rule_output(index) = NaN*index;
    plot(x', rule_output'); title('Each Rule''s Outputs');

figure;

subplot(221);
    mf_sum = sum(mf); normalized_mf = mf./mf_sum(ones(mf_n, 1), :);
    plot(x', normalized_mf'); axis([in_range 0 1.2]); title('Normalized Initial MFs');

subplot(222); 
    mf_sum = sum(mf); normalized_mf = mf./mf_sum(ones(mf_n, 1), :);
    plot(x', normalized_mf'); axis([in_range 0 1.2]); title('Normalized Final MFs');

subplot(223); tmp = [trn_error chk_error]; plot(tmp); title('Error Curves');
subplot(224); plot(step_size); title('Step Sizes');

