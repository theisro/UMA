%This function generate the resync filter with the same alghorithm of 
%fix_unsync_ir starting from an estimation of rx_speed in ppm

function irR = generate_resync_filter(rx_speed,N,Fs)

k = rx_speed/7.3288e+04;

f = linspace(0,Fs/2,N/2+1);

%Crea il filtro rifasatore utilizzando la sola componente di fase relativa
%al termine logaritmico riconducibile allo skew dei clock
pR = (k.*(f.*(log(f)-1)));
pR(1) = 0;
pR(N/2+1) = 0;
HR = exp(-j*pR);
HR = [HR flipdim(conj(HR(2:end-1)),2)];


%Trasformalo nel dominio del tempo
irR = ifft(HR);
irR = [ irR(N/2:end) irR(1:N/2-1) ];


end
