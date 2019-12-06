function [x, fs] = test_signal(k)
% x = test_signa(k)
% k = 1 - low-oscillatory signal    (length 512)
% k = 2 - high-oscillatory signal   (length 512)
% k = 3 - speech waveform ('I'm')   (length 2048)
% k = 4 - low-oscillatory signal    (length 250)
% k = 11 - low-oscillatory complex signal    (length 512) (like k = 1)
% k = 12 - high-oscillatory complex signal   (length 512) (like k = 2);

switch k
    
    case 1
        
        % Low-resonance signal - length 512
        
        N2 =  40; v2 = sin(0.10*pi*(1:N2)) .* blackman(N2)';
        N4 =  20; v4 = sin(0.20*pi*(1:N4)) .* blackman(N4)';
        N6 =  80; v6 = sin(0.05*pi*(1:N6)) .* blackman(N6)';
        
        N = 2^9; % N = 512;
        x = zeros(1,N);
        x(50+(1:N2)) = v2;
        x(300+(1:N4)) = x(300+(1:N4)) + v4;
        x(150+(1:N6)) = x(150+(1:N6)) + v6;

        fs = 1;
        
    case 2

        % High-resonance signal - length 512

        N1 = 160; v1 = sin(0.10*pi*(1:N1)) .* blackman(N1)';
        N3 =  80; v3 = sin(0.20*pi*(1:N3)) .* blackman(N3)';
        N5 = 320; v5 = sin(0.05*pi*(1:N5)) .* blackman(N5)';
        
        N = 2^9; % N = 512
        x = zeros(1,N);
        x(50+(1:N1)) = v1;
        x(200+(1:N3)) = x(200+(1:N3)) + v3;
        x(180+(1:N5)) = x(180+(1:N5)) + v5;
        
        fs = 1;

    case 3
        
        % Load speech signal - length 2048
        x = load('speech1.txt');
        fs = 16000;                
        x = x(:)';
                
    case 4
        
        % Low-resonance signal - length 256
        
        N2 =  40; v2 = sin(0.10*pi*(1:N2)) .* blackman(N2)';
        N4 =  20; v4 = sin(0.20*pi*(1:N4)) .* blackman(N4)';
        N6 =  80; v6 = sin(0.05*pi*(1:N6)) .* blackman(N6)';
        
        N = 256;
        x = zeros(1,N);
        x(50+(1:N2)) = v2;
        x(100+(1:N6)) = x(100+(1:N6)) + v6;
        x(200+(1:N4)) = x(200+(1:N4)) + v4;

        fs = 1;
        
        
    case 5

        % High-resonance signal - length 512
        
        N1 = 160; v1 = sin(0.10*pi*(1:N1)) .* blackman(N1)';
        N3 =  80; v3 = sin(0.20*pi*(1:N3)) .* blackman(N3)';
        N5 = 320; v5 = sin(0.05*pi*(1:N5)) .* blackman(N5)';
        
        N = 2^9; % N = 512
        x = zeros(1,N);
        
        x(10+(1:N1)) = v1;
        x(150+(1:N3)) = x(150+(1:N3)) + v3;
        x(180+(1:N5)) = x(180+(1:N5)) + v5;
        
        fs = 1;
        
        
    case 11
        
        % Low-resonance signal - length 512
        
        I = sqrt(-1);
        N2 =  40; v2 = exp(I*0.10*pi*(1:N2)) .* blackman(N2)';
        N4 =  20; v4 = exp(I*0.20*pi*(1:N4)) .* blackman(N4)';
        N6 =  80; v6 = exp(I*0.05*pi*(1:N6)) .* blackman(N6)';
        
        N = 2^9; % N = 512;
        x = zeros(1,N);
        x(50+(1:N2)) = v2;
        x(300+(1:N4)) = x(300+(1:N4)) + v4;
        x(150+(1:N6)) = x(150+(1:N6)) + v6;

        fs = 1;
        
    case 12

        % High-resonance signal - length 512

        I = sqrt(-1);
        N1 = 160; v1 = exp(I*0.10*pi*(1:N1)) .* blackman(N1)';
        N3 =  80; v3 = exp(I*0.20*pi*(1:N3)) .* blackman(N3)';
        N5 = 320; v5 = exp(I*0.05*pi*(1:N5)) .* blackman(N5)';
        
        N = 2^9; % N = 512
        x = zeros(1,N);
        x(50+(1:N1)) = v1;
        x(200+(1:N3)) = x(200+(1:N3)) + v3;
        x(180+(1:N5)) = x(180+(1:N5)) + v5;
        
        fs = 1;
        
end


        
        