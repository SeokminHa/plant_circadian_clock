%% ODE function assuming hill type form of regulation

function dC = ODE_rhythmic_degra(t,C)

global A K B phi

%  ODE function

dC= A * input(t)- K * (1+B* cos(2*pi/24*(t-phi)))*C;


    function [outputArg1] = input(inputArg1)
        input_d = readmatrix('mRNA_time_course.csv');
        outputArg1=max(makima(input_d(1,:),input_d(2,:),mod(inputArg1,24)),0);
    end

end

