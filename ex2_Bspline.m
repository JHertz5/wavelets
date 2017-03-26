close all;
clear all;
clc;

%Define program parameters, including length of signal in terms of samples,
%the sampling kernel, the number of iterations of the filter in signal, the
%time axis vector.
M = 256;
T = 2;
alpha = 1/128;

%As want to reproduce up to degree 3 need scaling function with 3+1
%vanishing moments, hence use dB4
sampling_kernel = 'dB4';
N = 5;
t = 0:(M-1);
t = t./T;


%Find the Dual basis of the B-Spline 
syms z;
temp = 2-z-z^-1;
P_sum = factor(1+temp+(10/16)*temp^2+(20/64)*temp^3);
H = expand(((z+1)^4)*(P_sum));
%Shift up by 3 as sum2poly only works on positive indices
dual_phi_T = sym2poly((z^3)*H);

%Calculate the B spline scaling function
[beta0, psi_T, xval] = wavefun('haar', 1);
phi_T = conv(beta0, conv(beta0,conv(beta0,conv(beta0, beta0))));
phi_T = phi_T(4:end-2);

figure(2)
subplot(1,2,1)
plot(phi_T), title('B3-Spline Scaling Function'), xlabel('Tap Number'), ylabel('Amplitude of Tap');
axis tight;
subplot(1,2,2)
plot(dual_phi_T), title('Dual basis of B-Spline Scaling Function'), xlabel('Tap Number'), ylabel('Amplitude of Tap');
axis tight;

%Create Sampling and Reconstruction matrices
S_sample = sampling_matrix2(dual_phi_T, T, M);
S_recon = sampling_matrix2(phi_T, T, M);

%Create c using the inner product between B spline and the polynomials
c = calc_c(S_recon, N, T, M);

%Initialise reproduced/reconstructed polynomial
rep_poly = zeros(size(c,1), M);

%Loop calculates the reconstructed time monomial and for orders 0 -(N-1). 
for j = 1:size(c,1)
    t1 = sprintf('Reconstructed monomial, order = %.1f',(j-1));
    t2 = sprintf('Error, order = %.1f',(j-1));
    for k = 1:size(S_sample,1)
        temp = alpha.*c(j, k).*S_sample(k,:);
        rep_poly(j,:) = rep_poly(j,:) + temp;
        
        %Plot the individual sampling kernels multiplied by their
        %respective coefficients
        %if(k<=60)
            figure(1)
            subplot(2, 2, j);
            plot(temp);
            axis tight;
            hold on;
        %end
       
        
    end
    %Calculate the error between the reconstructed monomial and the
    %original monomial. 
    error = abs((t.^(j-1)) - rep_poly(j,:));
    
    %Plot results
    figure(1)
    subplot(2, 2, j);
    plot(rep_poly(j,:), 'g'), title(t1), xlabel('Time'), ylabel('Reconstructed t^m');
    hold on;
    plot(t.^(j-1), 'r');
    axis tight
   
end
