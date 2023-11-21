function outData = oa_conv(inData,h,N)

    overlapBlock = zeros(N,1);

    %Prepare frequency domain FIR coeff
    H = fft([h(1:N);zeros(N,1)]);
     
    %Process the data block by block 
    D = length(inData);
    d = 1;
    while(d+N-1 <= D)
        inBlock = inData(d:d+N-1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Block process
        
        fd_inBlock = fft([inBlock;zeros(N,1)]);
        
        fd_convBlock = fd_inBlock .* H;
        
        convBlock = ifft(fd_convBlock);
        
        outBlock = overlapBlock + convBlock(1:N);
        
        overlapBlock = [ convBlock(N+1:2*N-1) ; 0 ];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        outData(d:d+N-1) = outBlock;
        d = d+N;
    end
    
end
