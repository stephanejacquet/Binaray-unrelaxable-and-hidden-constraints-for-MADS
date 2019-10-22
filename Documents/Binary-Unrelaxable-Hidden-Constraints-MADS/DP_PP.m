function DP_PP()

clear all;

format('long')

%%***********
%INPUT DATA
%%***********

cBlue = [28 96 243] ./ 255;
cRed = [239 11 3] ./ 255;
cGreen = [96 239 10] ./ 255;
cOrange = [241 126 32] ./ 255;
cViolet = [150 11 150] ./ 255;
cRose = [255 95 145] ./ 255;
cNoir = [0 0 0] ./ 255;

cLightBlue = [72 242 253] ./ 255;
cLightRed = [231 100 148] ./ 255;
cLightGreen = [56 220 190] ./ 255;
cLightOrange = [245 206 149] ./ 255;
cLightViolet = [187 133 189] ./ 255;


directory = 'Problems90/';
%directory = 'Problemes90_3_ordi_perso/';
n = 90; % number of problems
m = 3; % number of solvers
%solvers = {'NOMAD-Epane','NOMAD-Direction-Success','NOMAD-Gaussian','NOMAD-Uniform','NOMAD-Epane-exploration','NOMAD-Epane-denominateur'};
%solvers = {'NOMAD-Epane','NOMAD-Direction-Success','NOMAD-Gaussian','NOMAD-Epane-exploration','NOMAD-Epane-exploration-normalised'};
%solvers = {'NOMAD-Direction-Success','NOMAD-Epane-exploration'};
%solvers = {'NOMAD-Radius','NOMAD-Epane','NOMAD-Default','NOMAD-Gaussian','NOMAD-EpanRadius'};
%solvers = {'NOMAD-Epane','NOMAD-Uniform','NOMAD-Gaussian'};
%solvers = {'NOMAD-Direction-Success','NOMAD-Epane-exploration'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.9'};
solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.4','NOMAD-3.8.0.1.13'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.7','NOMAD-3.8.0.1.8','NOMAD-3.8.0.1.9','NOMAD-3.8.0.1.13'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.4','NOMAD-3.8.0.1.5','NOMAD-3.8.0.1.6','NOMAD-3.8.0.1.7','NOMAD-3.8.0.1.8','NOMAD-3.8.0.1.9','NOMAD-3.8.0.1.13'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.7','NOMAD-3.8.0.1.8','NOMAD-3.8.0.1.9'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.4','NOMAD-3.8.0.1.5','NOMAD-3.8.0.1.6'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.4','NOMAD-3.8.0.1.7'};
%solvers = {'NOMAD-3.8.0.1.5','NOMAD-3.8.0.1.8','NOMAD-3.8.0.1.13'};
%solvers = {'NOMAD-3.8.0.1.0','NOMAD-3.8.0.1.6','NOMAD-3.8.0.1.9'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.02','NOMAD-3.8.07','NOMAD-3.8.04','NOMAD-3.8.05','NOMAD-3.8.06'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.08','NOMAD-3.8.09','NOMAD-3.8.10'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.08','NOMAD-3.8.09','NOMAD-3.8.10'};

%solvers = {'NOMAD-3.8.0','NOMAD-3.8.10'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.03','NOMAD-3.8.04','NOMAD-3.8.05'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.02','NOMAD-3.8.03','NOMAD-3.8.04'};
%solvers = {'NOMAD-3.8.04','NOMAD-3.8.08'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.02','NOMAD-3.8.03','NOMAD-3.8.04','NOMAD-3.8.05','NOMAD-3.8.06','NOMAD-3.8.08'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.08','NOMAD-3.8.09','NOMAD-3.8.10'};
%solvers = {'NOMAD-3.8.0','NOMAD-3.8.08'};
%solvers = {'NOMAD-3.8.0m','NOMAD-3.8.11','NOMAD-3.8.12'};
%solvers = {'NOMAD-Epane','NOMAD-Gaussian'};
%solvers = {'NOMAD-Direction-Success','NOMAD-Epane'};
%colors = {cBlue,cRed,cGreen,cOrange,cViolet};
colors = {cBlue,cOrange,cViolet,cRed,cGreen,cRose,cNoir,cLightBlue};
markers = ['o','s', '^', 'v','d','>','<','^'];
lines = {'--','--','--','--','--','--','--','--'};
tau = [1e-1]; % list of convergence test tolerances
%tau = [1e-2,1e-3]; % list of convergence test tolerances
%tau = [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6];
l_tau = size(tau,2); % number of convergence test tolerances
logplot = 0;
%espa = [8,10,12];
espa = randperm(m,m)+ (10-m+3) + 8;

%%*******************
%READING INPUT DATA
%%*******************

probSize = 4*ones(1,n);
for i=1:n
	for j=1:m        
		pwd = strcat(directory,num2str(i),'/',solvers(j),'/res.txt');
		M{i,j} = dlmread(pwd{1});
		maxEvalProb((m*(i-1))+j) = M{i,j}(end,1);
	end
	probSize = [probSize,length(M{i,1}(1,:))-2];
end
probSizeSty = 8*ones(1,n/3); 
probSizeMDO = 7*ones(1,n/3); 
probSizeLock = 6*ones(1,n/3);
probSize = [probSizeSty probSizeMDO probSizeLock];
maxEval = max(maxEvalProb);

F = NaN(maxEval,m,n);

for i=1:n
	for j=1:m
		for k=1:length(M{i,j}(:,1))
			F(M{i,j}(k,1),j,i) = M{i,j}(k,2);
		end
		for l=2:maxEval
			if isnan(F(l,j,i))
				F(l,j,i) = F(l-1,j,i);
			end
			if F(l,j,i) > F(l-1,j,i);
				F(l,j,i) = F(l-1,j,i);
            end
        end
        idx = 1;
        while isnan(F(idx,j,i))
           idx = idx+1; 
        end
        fx0(i,j) = F(idx,j,i);
		%fx0(i,j) = F(1,j,i);
	end
end

%%************************
%COMPUTING t_PS matrices
%%************************

for l=1:l_tau
    
	t_PS{l} = zeros(n,m);
    
	for i=1:n
		fL = min(F(end,1:end,i));
		for j=1:m
			crit = fx0(i,j) - ((1-tau(l))*(fx0(i,j) - fL));
			k = find(F(:,j,i) <= crit,1);
			if (isempty(k))
				t_PS{l}(i,j) = NaN;
				t_PS_new{l}(i,j) = NaN;
			else
				t_PS{l}(i,j) = k;
				t_PS_new{l}(i,j) = k/(probSize(i)+1);
			end
		end
	end
	r_PS{l} = t_PS{l}./repmat(min(t_PS{l},[],2),1,m);
end

for l=1:l_tau
	t_PS_max(l) = max(max(t_PS_new{l}));
	t_PS_new{l}(isnan(t_PS_new{l})) = 2*t_PS_max(l);
	t_PS_new{l} = sort(t_PS_new{l});
	r_PS_max(l) = max(max(r_PS{l}));
	r_PS{l}(isnan(r_PS{l})) = 2*r_PS_max(l);
	r_PS{l} = sort(r_PS{l});
end

%%*********************
%Plotting PPs
%%*********************

for l=1:l_tau
	hl{l} = zeros(m,1);
	figure;
	hold on;
	grid on;
    
	for j=1:m
		[xs{l,j},ys{l,j}] = stairs(r_PS{l}(:,j),(1:n)/n);
		indx = find(xs{l,j}<=r_PS_max(l),1,'last');
		if isempty(indx)
			xs{l,j} = 1;
			ys{l,j} = 0;
        else ceil(r_PS_max(l))
			xs{l,j} = xs{l,j}(1:indx);
			ys{l,j} = ys{l,j}(1:indx);
		end
		if (xs{l,j}(1)==1)
			vv = find(xs{l,j}==1,1,'last');
			xs{l,j} = xs{l,j}(vv:end);
			ys{l,j} = ys{l,j}(vv:end);
		else
			xs{l,j} = [1;xs{l,j}(1);xs{l,j}];
			ys{l,j} = [0;0;ys{l,j}];
		end
		xs{l,j} = [xs{l,j};ceil(r_PS_max(l))];
		ys{l,j} = [ys{l,j};ys{l,j}(end)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:m
        nj = length(xs{l,j});
        tmp = linspace(xs{l,j}(1), xs{l,j}(nj), espa(j));
        for r=1:espa(j)
            idxtmp1 = find(xs{l,j} == tmp(r),1);
            if isempty(idxtmp1)
                idxX{l,j}(r) =  find(xs{l,j} >= tmp(r),1)-1;
            else
                idxX{l,j}(r) = idxtmp1;
            end  
        end
        option = [markers(j)];
        L{l,j} = ys{l,j}(idxX{l,j});
        if (logplot)
            plot(xs{l,j}(espa(j) + round(linspace(1,nj-espa(j)),espa(j))),log(L{l,j}),option,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7)
        else
            plot(tmp,L{l,j},option,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7)
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:m
		option = [lines{j}];
		if (logplot)
			hl{l}(j) = semilogx(xs{l,j},ys{l,j},option,'color',colors{j},'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7);
		else
			hl{l}(j) = plot(xs{l,j},ys{l,j},option,'color',colors{j},'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7);
        end
        
	end
	title ( strcat('precision \tau = 1e- ',num2str(-log10(tau(l)))) );
    rect1=[0.15, 0.75, 0.15, 0.15];
    %leg ('Epane','Default','Gaussian','Uniform','NOMAD-exploration','NOMAD-denominateur');
    %legend ('Epane','Default','Gaussian','NOMAD-exploration','NOMAD-normalised');
	%leg = legend ('Gaussian','Epan.(knn)','Default','Gaussian(knn)','Epan.');
    %leg = legend ('Gaussian','Epan.(knn)','Gaussian(knn)','Epan.');
    %leg = legend ('Distance à la cache','Direction de dernier succès');
    %leg = legend ('Epanechnikov','Direction of last success');
    %leg = legend ('Distance to cache','Direction of last success');
    %leg = legend ('Defaut','Ma version');
    %leg = legend ('Defaut','Max','Median no updates','Median with updates');
    %leg = legend ('Defaut','Max','Median no updates','Median with updates','h_max update','Max with update','med. h_max');
    %leg = legend ('Defaut','Mediane ','Violation','Maximum');
    %leg = legend ('Ma version5','Ma version8');
    %leg = legend ('Defaut','Mediane all updates');
    %leg = legend ('Defaut','Ma version7','Ma version8','Ma version9');
    %leg = legend ('Defaut','Ma version4','Ma version5','Ma version6');
    %leg = legend ('Defaut','1)','2)','3)');
    %leg = legend ('Defaut','2)+sans mise à l échelle','2)+avec mise à l échelle');
    %leg = legend ('Defaut','Ma version4','Ma version7');
    %leg = legend ('Defaut','Ma version5','Ma version8');
    %leg = legend ('Defaut','Ma version6','Ma version9');
    %leg = legend ({'Ma version6','Ma version9'}, 'Position',rect1);
    %leg = legend ('Epanechnikov','Uniforme','Gaussian');
	%set(legend,'interpreter','none');
	xlabel ( 'ratio' );
	ylabel ( 'ratio of solved problems' );  
    %xlabel ( 'performance ratio \alpha' );
    %ylabel ( 'ratio of solved problems' ); 
	if (logplot)
		axis([1 ceil(r_PS_max(l)) 0 1.01]);
		twop = floor(log2(1.1*r_PS_max(l)));
		set(gca,'XTick',2.^[0:twop])
    else
        if (ceil(r_PS_max(l))==1)
          axis([1 2 0 1.01]); 
        else
            
		axis([1 ceil(r_PS_max(l)) 0 1.01]);
        end
	end
	hold off;
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',20)
	%print(strcat('WeldedPP1e-',num2str(-log10(tau(l)))), '-dpng');
     saveas(gcf,sprintf('PP%d.pdf',l));
end

%%*********************
%Plotting DPs
%%*********************

for l=1:l_tau
	hll{l} = zeros(m,1);
	figure;
	hold on;
	grid on;
    for j=1:m
        [xs2{l,j},ys2{l,j}] = stairs(t_PS_new{l}(:,j),(1:n)/n);
        indx = find(xs2{l,j}<=t_PS_max(l),1,'last');
        if isempty(indx)
            xs2{l,j} = 0;
            ys2{l,j} = 0;
        else
            xs2{l,j} = xs2{l,j}(1:indx);
            ys2{l,j} = ys2{l,j}(1:indx);
        end
        if (xs2{l,j}(1)~=0)
            xs2{l,j} = [0;xs2{l,j}(1);xs2{l,j}];
            ys2{l,j} = [0;0;ys2{l,j}];
        end
        xs2{l,j} = [xs2{l,j};ceil(t_PS_max(l))];
        ys2{l,j} = [ys2{l,j};ys2{l,j}(end)];
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:m
        nj = length(xs2{l,j});
        tmp2 = linspace(xs2{l,j}(1), xs2{l,j}(nj), espa(j));
        for r=1:espa(j)
            idxtmp1 = find(xs2{l,j} == tmp2(r),1);
            if isempty(idxtmp1)
                idxX2{l,j}(r) =  find(xs2{l,j} >= tmp2(r),1)-1;
            else
                idxX2{l,j}(r) = idxtmp1;
            end
        end
        option = [markers(j)];
        L2{l,j} = ys2{l,j}(idxX2{l,j});
        plot(tmp2,L2{l,j},option,'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7)
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:m
        option = [lines{j}];
        hll{l}(j) = plot(xs2{l,j},ys2{l,j},option,'color',colors{j},'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',colors{j},'MarkerSize',7);
    end
	title ( strcat('precision \tau = 1e- ',num2str(-log10(tau(l)))) );
	axis ( [1 ceil(t_PS_max(l))  0  1.01] );
    %leg ('Epane','Default','Gaussian','Uniform','NOMAD-exploration','NOMAD-denominateur');
    %legend ('Epane','Default','Gaussian','NOMAD-exploration','NOMAD-normalised');
	%leg = legend ('Gaussian','Epan.(knn)','Default','Gaussian(knn)','Epan.');
    %leg = legend ('Epanechnikov','NOMAD-Defaut');
    %leg = legend ('Epanechnikov','Uniforme','Gaussien');
    %leg = legend ('Defaut','Ma version');
    %leg = legend ('Defaut','Max','Median no updates','Median with updates','h_max update','Max with update','med. h_max');
    %leg = legend ('Défaut','Médiane ','1-ère Violation','Maximum');
    %leg = legend ('Default','Median ','1st Violation','Maximum');
    %leg = legend ('Defaut','Mediane all updates');
    %leg = legend ('$h_1$','$h_7$');
    leg = legend ('Default','$\tilde{h}_1$','$\tilde{h}_7$');
    %leg = legend ('Default','$\tilde{h}_1$','$\tilde{h}_2$','$\tilde{h}_3$');
    %leg = legend ('Default','$\tilde{h}_4$','$\tilde{h}_5$','$\tilde{h}_6$','$\tilde{h}_7$');
    %leg = legend ('Default','$h_4$','$h_5$','$h_6$','$h_7$');
    %leg = legend ('Ma version6','Ma version9');
    %leg = legend ('Defaut','Ma version7','Ma version8','Ma version9');
    %leg = legend ('Defaut','Ma version4','Ma version5','Ma version6');
    %leg = legend ('Defaut','Ma version4','Ma version7');
    %leg = legend ('Defaut','2)+sans mise à l échelle','2)+avec mise à l échelle');
    %leg = legend ('Defaut','Ma version5','Ma version8');
    %leg = legend ('Defaut','Ma version6','Ma version9');
    %leg = legend ('Direction de dernier succès','Epanechnikov');
    %leg = legend ('Epanechnikov','Direction of last success'); 
    %leg = legend ('Gaussian','Epan.(knn)','Gaussian(knn)','Epan.');
    %leg = legend ('Default','Gaussian','GaussRadius');
	%set(legend,'interpreter','latex');
    %leg = legend ('Default','$\tilde{h}$',' $\tilde{h} \mbox{ and } h$');
    set(legend,'interpreter','latex');
	xlabel ( 'groups of (n+1) evaluations' );
	ylabel ( 'ratio of solved problems' );  
	%xlabel ( 'groupes de (n+1) évaluations' );
	%ylabel ( 'ratio de problèmes résolus' );    
	hold off;
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',20)
	%print(strcat('WeldedDP1e-',num2str(-log10(tau(l)))), '-dpng');
    saveas(gcf,sprintf('DP%d.pdf',l));
end

