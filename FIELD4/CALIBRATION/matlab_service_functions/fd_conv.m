function c = fd_conv(a,b)
 
la = length(a);
lb = length(b);
lc = la + lb - 1;

pa = [a(:) ; zeros(lc-la,1)];
pb = [b(:) ; zeros(lc-lb,1)];

fa = fft(pa);
fb = fft(pb);
fc = fa.*fb;

c = ifft(fc);
end

