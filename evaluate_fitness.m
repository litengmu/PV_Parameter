%% 光伏模块参数评估函数合集

function  obj = evaluate_fitness(x,func_flag)

switch  func_flag
    case 1
        % SDM
        obj = PV_model_single(x);
    case 2
        % DDM
        obj = PV_model_double(x);
    case 3
        %  Photowatt-PWP201 PV module
        obj = PV_module1(x);
    case 4
        %  STM6-40:36 module
        obj = PV_module2(x);
    case 5
        %  STP6-120:36 PV module
        obj = PV_module3(x);
    case 6
        %  SM55 T25 PV module
        obj = SM55_T25(x);
    case 7
        %  SM55 T40 PV module
        obj = SM55_T40(x);
    case 8
        %  SM55 T60 PV module
        obj = SM55_T60(x);
    case 9
        %  SM55 T25 G200 PV module
        obj = SM55_G200(x);
    case 10
        %  SM55 T25 G400 PV module
        obj = SM55_G400(x);
    case 11
        %  SM55 T25 G600 PV module
        obj = SM55_G600(x);
    case 12
        %  SM55 T25 G800 PV module
        obj = SM55_G800(x);
    case 13
        %  SM55 T25 G1000 PV module
        obj = SM55_G1000(x);
    case 14
        %  SM55 T25 PV module
        obj = KC200GT_T25(x);
    case 15
        %  SM55 T40 PV module
        obj = KC200GT_T50(x);
    case 16
        %  SM55 T60 PV module
        obj = KC200GT_T75(x);
    case 17
        %  SM55 T25 G200 PV module
        obj = KC200GT_G200(x);
    case 18
        %  SM55 T25 G400 PV module
        obj = KC200GT_G400(x);
    case 19
        %  SM55 T25 G600 PV module
        obj = KC200GT_G600(x);
    case 20
        %  SM55 T25 G800 PV module
        obj = KC200GT_G800(x);
    case 21
        %  SM55 T25 G1000 PV module
        obj = KC200GT_G1000(x);
    case 22
        %  SM55 T25 G1000 PV module
        obj = PV_model_triple(x);
    case 23
        %  SM55 T25 G1000 PV module
        obj = SM55_STC(x);
    case 24
        %  SM55 T25 G1000 PV module
        obj = KC200GT_STC(x);
end

end

%% ****************************************************************
% -------------------------  SDM
function  result = calculate_objective_single(x,V_L,I_L)

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 33.0;		%  33 centi-degree RTC France solar cell;  25 centi-degree PVM752 GaAs solar cell
V_t = k * T / q;
result = I_ph - I_SD * ( exp( (V_L + I_L*R_s) / (V_t*n) ) - 1.0 ) - ( (V_L + I_L*R_s)/R_sh ) - I_L;
end

function obj = PV_model_single(x)

a = load('RTC France solar cell.txt');  %PVM752 GaAs solar cell(SDM).txt
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_single(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end
% -------------------------

%% ****************************************************************
% -----------------------  DDM
function  result = calculate_objective_double(x,V_L,I_L)

I_ph	= x(1);
I_SD1	= x(2);
R_s		= x(3);
R_sh	= x(4);
n1		= x(5);
I_SD2	= x(6);
n2		= x(7);

q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 33.0;		%  33 centi-degree RTC France solar cell;  25 centi-degree PVM752 GaAs solar cell

result = I_ph - I_SD1 * ( exp( (q*(V_L + I_L*R_s)) / (n1*k*T) ) -1.0 ) - I_SD2 * ( exp( (q*(V_L + I_L*R_s)) / (n2*k*T) ) -1.0 ) - ( (V_L + I_L*R_s)/R_sh ) - I_L;

end

function obj = PV_model_double(x)

a = load('RTC France solar cell.txt'); %PVM752 GaAs solar cell(DDM).txt
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_double(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end


% -----------------------  DDM
function  result = calculate_objective_triple(x,V_L,I_L)

I_ph	= x(1);
I_SD1	= x(2);
R_s		= x(3);
R_sh	= x(4);
n1		= x(5);
I_SD2	= x(6);
n2		= x(7);
I_SD3	= x(8);
n3		= x(9);

q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 33.0;		%  33 centi-degree RTC France solar cell;  25 centi-degree PVM752 GaAs solar cell

result = I_ph - I_SD1 * ( exp( (q*(V_L + I_L*R_s)) / (n1*k*T) ) -1.0 ) - I_SD2 * ( exp( (q*(V_L + I_L*R_s)) / (n2*k*T) ) -1.0 )-I_SD3 * ( exp( (q*(V_L + I_L*R_s)) / (n3*k*T) ) -1.0 ) - ( (V_L + I_L*R_s)/R_sh ) - I_L;

end

function obj = PV_model_triple(x)

a = load('RTC France solar cell.txt'); %PVM752 GaAs solar cell(DDM).txt
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_triple(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end



%% ****************************************************************
% -------------------------  PV module: Photowatt-PWP201 PV module
function  result = calculate_objective_module1(x,V_L,I_L)

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 45.0;		%   45 centi-degree
V_t = k * T / q;
NS=1;
NP=1;
result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;

end

function obj = PV_module1(x)

a = load('Photowatt-PWP201 PV module.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_module1(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end


%% ****************************************************************
% -------------------------  PV module:STM6-40:36 module
function  result = calculate_objective_module2(x,V_L,I_L)%%%

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 51.0;		%   51 centi-degree
V_t = k * T / q;
NS=36;
NP=1;
result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;

end

function obj = PV_module2(x)

a = load('STM6-40:36 module.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_module2(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end
% -------------------------  PV module

%% ****************************************************************
% -------------------------  PV module: STP6-120:36 PV module
function  result = calculate_objective_module3(x,V_L,I_L)%%%

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + 55.0;		%   55 centi-degree
V_t = k * T / q;
NS=36;
NP=1;
result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;
end

function obj = PV_module3(x)

a = load('STP6-120:36 PV module.txt');
actual_V_data =  a(:,1);
actual_I_data =  a(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = calculate_objective_module3(x,actual_V_data(j), actual_I_data(j));
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end


%% Below SM55
function  result = SM55(x,V_L,I_L,t)%%%

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + t;		%   55 centi-degree
V_t = k * T / q;
NS=36;
NP=1;
%result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;
result = NP*I_ph  - NP*I_SD * ( exp( (V_L*NP + (I_L*R_s*NS)) / (V_t*n*NS*NP) ) - 1.0 ) -  (V_L*NP + (I_L*R_s*NS))/(R_sh*NS)  - I_L;

end

function obj = SM55_T25(x)
t=25;
load("SM55.mat",'SM55_T25');
actual_V_data =  SM55_T25(:,1);
actual_I_data =  SM55_T25(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end

function obj = SM55_T40(x)
t=40;
load("SM55.mat",'SM55_T40');
actual_V_data =  SM55_T40(:,1);
actual_I_data =  SM55_T40(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end

function obj = SM55_T60(x)
t=60;
load("SM55.mat",'SM55_T60');
actual_V_data =  SM55_T60(:,1);
actual_I_data =  SM55_T60(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = SM55_G200(x)
t=25;
load("SM55.mat",'SM55_G200');
actual_V_data =  SM55_G200(:,1);
actual_I_data =  SM55_G200(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = SM55_G400(x)
t=25;
load("SM55.mat",'SM55_G400');
actual_V_data =  SM55_G400(:,1);
actual_I_data =  SM55_G400(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = SM55_G600(x)
t=25;
load("SM55.mat",'SM55_G600');
actual_V_data =  SM55_G600(:,1);
actual_I_data =  SM55_G600(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = SM55_G800(x)
t=25;
load("SM55.mat",'SM55_G800');
actual_V_data =  SM55_G800(:,1);
actual_I_data =  SM55_G800(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = SM55_G1000(x)
t=25;
load("SM55.mat",'SM55_G1000');
actual_V_data =  SM55_G1000(:,1);
actual_I_data =  SM55_G1000(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end


%% Below KC200GT
function  result = KC200GT(x,V_L,I_L,t)%%%

I_ph = x(1);
I_SD = x(2);
R_s	 = x(3);
R_sh = x(4);
n	 = x(5);
q = 1.60217646e-19;
k = 1.3806503e-23;
T = 273.15 + t;		%   55 centi-degree
V_t = k * T / q;
NS=54;
NP=1;
%result = NP*I_ph - NP*I_SD * ( exp( (V_L/NS + (I_L*R_s/NP)) / (V_t*n) ) - 1.0 ) - ( NP*(V_L/NS + (I_L*R_s/NP))/R_sh ) - I_L;
% result = I_ph - I_SD * ( exp( (V_L + (I_L*R_s)) / (V_t*n*NS) ) - 1.0 ) -  (V_L + (I_L*R_s))/R_sh  - I_L;
result = NP*I_ph  - NP*I_SD * ( exp( (V_L*NP + (I_L*R_s*NS)) / (V_t*n*NS*NP) ) - 1.0 ) -  (V_L*NP + (I_L*R_s*NS))/(R_sh*NS)  - I_L;

end

function obj = KC200GT_T25(x)
t=25;
load("KC200GT.mat",'KC200GT_T25');
actual_V_data =  KC200GT_T25(:,1);
actual_I_data =  KC200GT_T25(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end

function obj = KC200GT_T50(x)
t=50;
load("KC200GT.mat",'KC200GT_T50');
actual_V_data =  KC200GT_T50(:,1);
actual_I_data =  KC200GT_T50(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);

end

function obj = KC200GT_T75(x)
t=75;
load("KC200GT.mat",'KC200GT_T75');
actual_V_data =  KC200GT_T75(:,1);
actual_I_data =  KC200GT_T75(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_G200(x)
t=25;
load("KC200GT.mat",'KC200GT_G200');
actual_V_data =  KC200GT_G200(:,1);
actual_I_data =  KC200GT_G200(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_G400(x)
t=25;
load("KC200GT.mat",'KC200GT_G400');
actual_V_data =  KC200GT_G400(:,1);
actual_I_data =  KC200GT_G400(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_G600(x)
t=25;
load("KC200GT.mat",'KC200GT_G600');
actual_V_data =  KC200GT_G600(:,1);
actual_I_data =  KC200GT_G600(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_G800(x)
t=25;
load("KC200GT.mat",'KC200GT_G800');
actual_V_data =  KC200GT_G800(:,1);
actual_I_data =  KC200GT_G800(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_G1000(x)
t=25;
load("KC200GT.mat",'KC200GT_G1000');
actual_V_data =  KC200GT_G1000(:,1);
actual_I_data =  KC200GT_G1000(:,2);
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = sqrt(fitness/data_len);
end

function obj = KC200GT_STC(x)
t=25;
actual_V_data =  [0;26.3;32.9];
actual_I_data =  [8.21;7.61;0];
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = KC200GT(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = fitness ;

end

function obj = SM55_STC(x)
t=25;
actual_V_data =  [0;17.4;21.7];
actual_I_data =  [3.45;3.15;0];
data_len = length(actual_V_data);
for j=1:data_len
    error_value(j) = SM55(x,actual_V_data(j), actual_I_data(j),t);
end
fitness = sum(error_value.^2);
obj = fitness ;

end