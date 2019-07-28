clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77GHz
% Max Range = 200m
% Range Resolution = 1 m
% Max Velocity = 100 m/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%speed of light = 3e8
%% User Defined Range and Velocity of target
% *%TODO* :
% define the target's initial position and velocity. Note : Velocity
% remains contant
target_start_position = 60; %target range =100m
target_velocity = -50; %target doppler velocity = -20m/s


%% FMCW Waveform Generation
% *%TODO* :
%Design the FMCW waveform by giving the specs of each of its parameters.
% Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
% chirp using the requirements above.

%Operating carrier frequency of Radar 
fc= 77e9;             %carrier freq


max_range=200;
range_resolution=1; % range resolution can determine bandwith
max_velocity=100;  % max velocity can determine Ts_x-axis
c=3e8;

wavelength = c/fc; %wavelength

bandwidth= c/(2*range_resolution);  %d_res = c/2*Bsweep

%Range= c*Ts*fb / (2* Bsweep)
%for max range, we assume fb = Bsweep/5.5
%so Ts= 5.5*2* maxRange / c
Tchirp= 5.5*2*max_range/c;

slope=bandwidth/Tchirp;

%% for 128 chrips, each chrips send 1024 times.

%The number of chirps in one sequence. Its ideal to have 2^ value for the ease of running the FFT
%for Doppler Estimation. 
Nd=128;                   % #of doppler cells OR #of sent periods % number of chirps

%The number of samples on each chirp. 
Nr=1024;                  %for length of time OR # of range cells

% Timestamp for running the displacement scenario for every sample on each
% chirp
% y = linspace(x1,x2,n) generates n points. The spacing between the points is (x2-x1)/(n-1).
t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples

%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));



%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

target_position=target_start_position;

for i=1:length(t)         
    
    % *%TODO* :
    %For each time stamp update the Range of the Target for constant velocity. 
    target_position= target_start_position + target_velocity*t(i);
    
    time_delay=(2*target_position)/c; %compute time_delay.
    
    % *%TODO* :
    %For each time sample we need update the transmitted and
    %received signal. 
    Tx(i) = cos(2*pi*(fc*t(i)+slope*t(i)*t(i)/2));
    Rx(i) = cos(2*pi*(fc*(t(i)-time_delay)+slope*(t(i)-time_delay)*(t(i)-time_delay)/2));
    
    % *%TODO* :
    %Now by mixing the Transmit and Receive generate the beat signal
    %This is done by element wise matrix multiplication of Transmit and
    %Receiver Signal
    Mix(i) = Tx(i)*Rx(i);
    
end

%% RANGE MEASUREMENT

 % *%TODO* :
%reshape the vector into Nr*Nd array. Nr and Nd here would also define the size of
%Range and Doppler FFT respectively.
Mix = reshape(Mix,[Nr,Nd]);

 % *%TODO* :
%run the FFT on the beat signal along the range bins dimension (Nr) and
%normalize.
Y=fft(Mix,Nr,1);

 % *%TODO* :
% Take the absolute value of FFT output
P2=abs(Y/Nr);

 % *%TODO* :
% Output of FFT is double sided signal, but we are interested in only one side of the spectrum.
% Hence we throw out half of the samples.
P1=P2(1:Nr/2);

P1=fftshift(P1); %important

%plotting the range
figure ('Name','Range from First FFT')
subplot(2,1,1)

% *%TODO* :
% plot FFT output 
range = linspace(-200,200,Nr/2)*((Nr/2)/400);

plot(range, P1);
axis ([0 200 0 1]);


%% RANGE DOPPLER RESPONSE
% The 2D FFT implementation is already provided here. This will run a 2DFFT
% on the mixed signal (beat signal) output and generate a range doppler
% map.You will implement CFAR on the generated RDM


% Range Doppler Map Generation.

% The output of the 2D FFT is an image that has reponse in the range and
% doppler FFT bins. So, it is important to convert the axis from bin sizes
% to range and doppler based on their Max values.

Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure,surf(doppler_axis,range_axis,RDM);

%% CFAR implementation

%Slide Window through the complete Range Doppler Map

% *%TODO* :
%Select the number of Training Cells in both the dimensions.
Tr=10;
Td=8;

% *%TODO* :
%Select the number of Guard Cells in both dimensions around the Cell under 
%test (CUT) for accurate estimation
Gr=4;
Gd=4;

% *%TODO* :
% offset the threshold by SNR value in dB
offset=6;

% *%TODO* :
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


% *%TODO* :
%design a loop such that it slides the CUT across range doppler map by
%giving margins at the edges for Training and Guard Cells.
%For every iteration sum the signal level within all the training
%cells. To sum convert the value from logarithmic to linear using db2pow
%function. Average the summed values for all of the training
%cells used. After averaging convert it back to logarithimic using pow2db.
%Further add the offset to it to determine the threshold. Next, compare the
%signal under CUT with this threshold. If the CUT level > threshold assign
%it a value of 1, else equate it to 0.

   % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
   % CFAR



% *%TODO* :
% The process above will generate a thresholded block, which is smaller 
%than the Range Doppler Map as the CUT cannot be located at the edges of
%matrix. Hence,few cells will not be thresholded. To keep the map size same
% set those values to 0. 
RESULT=zeros(Nr/2,Nd);

for i= Tr+Gr+1: Nr/2-(Gr+Tr)
    for j= Td+Gd+1:Nd-(Gd+Td)
        noise_level=0;
        train_cell=0;
        %sliding the grid
        for p= i-(Tr+Gr):i+Tr+Gr
            for q= j-(Td+Gd):j+Td+Gd
                if(abs(i-p)>Gr||abs(j-q)>Gd)
                    noise_level=noise_level+db2pow(RDM(p,q));
                    train_cell=train_cell+1;
                end
            end
        end
        
        noise_average=noise_level/train_cell;
        noise_threshold=pow2db(noise_average);
        noise_threshold=noise_threshold+offset;
        
        CUT=RDM(i,j);
        if(CUT<noise_threshold)
            RESULT(i,j)=0;
        else
            RESULT(i,j)=1;
        end
    end
end



% *%TODO* :
%display the CFAR output using the Surf function like we did for Range
%Doppler Response output.
figure,surf(doppler_axis,range_axis,RESULT);
colorbar;


 
 