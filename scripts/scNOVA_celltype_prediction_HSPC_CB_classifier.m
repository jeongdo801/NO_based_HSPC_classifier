%%-------------------------------------------------------------------------------------
%% Prediction based on PLS (cord blood) 
%%-------------------------------------------------------------------------------------
%%This is to infer cell type of HSPC using PLS-DA (Genebody matrix will be used)
%%Original script : /Users/jeong/Documents/Strand_Seq/Computational_code_HSPC_script_PLS_Pred
%%Last update : 2021-11-19
%%Dependency : Matlab

%%model: /Users/jeong/Downloads/HSPC_matlab/matlab_scMNase_final_model_CB.mat	

%%What's new in this analysis flow
%% Model : pls(X_GB_RPM_log_new_auto(:, Hit> 0), Y_new_auto, 100); %%175 X 18851 matrix (ref matrix)
%% Strand_seq data : will be RPM normalized and 18851 genes will be extracted and autoscaled using mean and stdev of ref matrix
%Hit = PLS_cor_GB_RPM_log_new.vip>1.2451;
Hit = (Hit5> 0 & RNA>0);


data_GB_Strand = data_GB_Strand_CNnorm_sc_stringent_ploidyassignR;
[m,n] = size(data_GB_Strand);
data_GB_Strand_RPM = zeros(m,n);

for i = 1:n
	data_GB_Strand_RPM(:,i) = data_GB_Strand(:,i)*1000000/sum(data_GB_Strand(:,i));
end

hist(data_GB_Strand_RPM(:,1),200);
hist(data_GB_Strand_RPM(:,2),200);

data_GB_Strand_RPM_log = log2(data_GB_Strand_RPM+1);
      
data_GB_Strand_RPM_log_sub = data_GB_Strand_RPM_log(data_GB_RPM_log_CV_new>0 & chr<=22,:);    
X_GB_Strand_RPM_log = data_GB_Strand_RPM_log_sub';


[mstrand, nstrand] = size(X_GB_Strand_RPM_log);
Ref_mean = mean(X_GB_RPM_log_new);
Ref_std = std(X_GB_RPM_log_new);
X_GB_Strand_RPM_log_auto = (X_GB_Strand_RPM_log-Ref_mean(ones(mstrand,1),:))./Ref_std(ones(mstrand,1),:);

[ypred_train, that1_train] = Pred_PLS(X_GB_RPM_log_new_auto(:, Hit> 0), Y_new_auto, X_GB_RPM_log_new_auto(:, Hit>0), 10);
[ypred_test, that1_test] = Pred_PLS(X_GB_RPM_log_new_auto(:, Hit> 0), Y_new_auto, X_GB_Strand_RPM_log_auto(:, Hit>0), 10);


ypred_max_train = zeros(length(ypred_train),1);
for i = 1:length(ypred_train)
	[M,I] = max(ypred_train(i,:));
	ypred_max_train(i,1)=I;
end

ypred_max_test = zeros(length(ypred_test),1);
for i = 1:length(ypred_test)
	[M,I] = max(ypred_test(i,:));
	ypred_max_test(i,1)=I;
end	




ypred_max_test_rank12 = zeros(length(ypred_test),2);
for i = 1:length(ypred_test)
	[M,I] = sort(ypred_test(i,:));
	ypred_max_test_rank12(i,1)=I(8);
	ypred_max_test_rank12(i,2)=I(7);
end	

ypred_max_test_rank12 = zeros(length(ypred_test),2);
for i = 1:length(ypred_test)
	[M,I] = sort(ypred_test(i,:));
	ypred_max_test_rank12(i,1)=M(8);
	ypred_max_test_rank12(i,2)=M(7);
end	
