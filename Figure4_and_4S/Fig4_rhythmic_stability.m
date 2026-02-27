
function Fig4_rhythmic_stability(ran_num)

rng(ran_num);

global A K B phi % 4 parameters



result1=[]; % save the parameter
result2=[]; % save the error


protein=readmatrix('protein_time_course.csv');
protein_time = protein(1,:);
protein_value = protein(2,:);

for wb=1:100000
    
    A=10*rand();
    K=10*rand();
    B=rand();
    phi = 24* rand();
    
    days=7;
    plevel = [];
    C1=[0];
    for j=1:days
        light = 1;
        tspan = 24*(j-1):1:24*(j-1)+24;
        options = odeset('NonNegative',1);
        [T1,C1] = ode15s(@(t,C) ODE_rhythmic_degra(t,C),tspan,C1(end,:),options);
        if j==days
            plevel = [plevel; C1];
        end
    end
    protein_simul=plevel(protein_time+1);
    
    pmax=max(protein_simul);
    
    protein_simul = protein_simul./pmax;
    error = norm(protein_value-protein_simul')
    
    result1=[result1; [A K B phi]];
    result2=[result2; [error,protein_simul']];

    csvwrite(['para_',num2str(ran_num),'.csv'],result1)
    csvwrite(['error_',num2str(ran_num),'.csv'],result2)
    end
    
    
end