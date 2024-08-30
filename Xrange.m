
%  SDM RTC France solar cell
if func_flag==1
    Xmin = [ 0.0  0.0      0.0  0.0    1.0];
    Xmax = [ 1.0  1.0e-06  0.5  100.0  2.0];
%     Xmin = [0.760775662  0.323154e-06  0.03637551  53.72563852    1.481225178];
%     Xmax = [0.760775662  0.323154e-06  0.03637551  53.72563852    1.481225178];

    D = 5;
    known_optimal   = 0.0;
end
 
%  DDM RTC France solar cell
if func_flag==2
    Xmin = [ 0.0  0.0      0.0   0.0   1.0  0.0     1.0];
    Xmax = [ 1.0  1.0e-06  0.5  100.0  2.0  1.0e-06  2];
%     Xmin = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
%     Xmax = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
    D = 7;
    known_optimal   = 0.0;
end

%  Photowatt-PWP201 PV module
if func_flag==3
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2.0  50.0e-06  2.0  2000.0  50.0];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
  
    D = 5;
    known_optimal   = 0.0;
end


% STM6-40:36 module
if func_flag==4
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 8  50.0e-06  0.36  1000  60.0];
    D = 5;
    known_optimal   = 0.0;
end
 
%  STP6-120:36 PV module
if func_flag==5
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 8  50.0e-06  0.36  1500  50.0];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
  
    D = 5;
    known_optimal   = 0.0;
end 

%  SM55 T25 PV module
if func_flag==6
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

%  SM55 T50 PV module
if func_flag==7
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(40-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06   2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

%  SM55 T75 PV module
if func_flag==8
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(60-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06   2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% SM55 T25 G200 PV module
if func_flag==9
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(200/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% SM55 T25 G400 PV module
if func_flag==10
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(400/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% SM55 T25 G600 PV module
if func_flag==11
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(600/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% SM55 T25 G800 PV module
if func_flag==12
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(800/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% SM55 T25 G1000 PV module
if func_flag==13
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 


%  KC200GT T25 PV module
if func_flag==14
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

%  KC200GT T50 PV module
if func_flag==15
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(50-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06   2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

%  KC200GT T75 PV module
if func_flag==16
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(75-25);
    I_sc_2=(1000/1000)*(I_sc_std+temp_c*(75-25));
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06   2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% KC200GT T25 G200 PV module
if func_flag==17
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(200/1000)+temp_c*(25-25);
    I_sc_2=(200/1000)*(I_sc_std+temp_c*(25-25));
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% KC200GT T25 G400 PV module
if func_flag==18
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(400/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% KC200GT T25 G600 PV module
if func_flag==19
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(600/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% KC200GT T25 G800 PV module
if func_flag==20
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(800/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

% KC200GT T25 G1000 PV module
if func_flag==21
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 5000  4];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 


if func_flag==22

    Xmin = [ 0.0  0.0      0.0   0.0   1.0  0.0     1.0    0.0     1.0     ];
    Xmax = [ 1.0  1.0e-06  0.5  100.0  2.0  1.0e-06  2   1.0e-06  2       ];
%     Xmin = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
%     Xmax = [0.760781  0.225974e-06  0.03674  55.485441 1.451017 0.749346e-06 2];
    D = 9;
    known_optimal   = 0.0;

end 

if func_flag==23
    I_sc_std=3.45;
    temp_c=1.2e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 1000  1.5];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 

if func_flag==24
    I_sc_std=8.21;
    temp_c=3.18e-3;
    I_sc=I_sc_std*(1000/1000)+temp_c*(25-25);
    Xmin = [ 0.0  0.0       0.0  0    1.0];
    Xmax = [ 2*I_sc  100e-06  2 1000  1.5];
%     Xmin = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];
%     Xmax = [1.030514  3.482263e-06  1.201271  981.982240  48.642835];    
    D = 5;
    known_optimal   = 0.0;
end 
