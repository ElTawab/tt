clear;
clc;
close all;

%data bits 
bits_data = randi([0,1],16,1000);      %random data to be repetition coded
data_rep = zeros(16,1000*3);           %data to be transmitted

%% coding
%reprtition code
for k = 1:16
   counter = 1;
   for i=1:length(bits_data)
       if bits_data(k,i)== 0
           data_rep(k,counter:counter+2)= [0 0 0];
       elseif bits_data(k,i)== 1
           data_rep(k,counter:counter+2)= [1 1 1];
       end
       counter = counter + 3;
   end
end
%no code
data_no_code = bits_data;

%% Interleaver
data_intrlvd_no_code = matintrlv(data_no_code,4,4);     %no code
data_intrlvd_rep = matintrlv(data_rep,4,4);             %repetition

%% Mapper
%% no code
%*********************************QPSKs********************************%
QPSK_bits_no_code=zeros(16,length(data_intrlvd_no_code)/2);
binary_4 = [1 1; 0 1; 0 0; 1 0 ];
for j=1:16
    for i=1:2:length(data_intrlvd_no_code)-1
        for k = 1 : 4
            if(data_intrlvd_no_code(j,i:i+1) == binary_4(k,:) )
                 QPSK_bits_no_code(j,floor(i/2)+1) =  sqrt(2)*exp(1i*((pi/4)+(pi/2)*(k-1))); 
            end
        end
    end
end

%*********************************16QAM********************************%
QAM16_bits_no_code = zeros(16,length(data_intrlvd_no_code)/4);
binary_16_3  = [1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];
binary_16_1  = [1 1 1 0; 1 1 1 1; 1 1 0 1; 1 1 0 0];
binary_16_n1 = [0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0];
binary_16_n3 = [0 0 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0];

for j=1:16
    for i = 1:4:length(data_intrlvd_no_code)-3
        d = 3;
          for k = 1 : 4
            if(data_intrlvd_no_code(j,i:i+3) == binary_16_3(k,:) )
                 QAM16_bits_no_code(j,floor(i/4)+1) = 3+1i*d; 
            end
            if(data_intrlvd_no_code(j,i:i+3) == binary_16_1(k,:) )
                 QAM16_bits_no_code(j,floor(i/4)+1) = 1+1i*d; 
            end
            if(data_intrlvd_no_code(j,i:i+3) == binary_16_n1(k,:) )
                 QAM16_bits_no_code(j,floor(i/4)+1) = -1+1i*d; 
            end
            if(data_intrlvd_no_code(j,i:i+3) == binary_16_n3(k,:) )
                 QAM16_bits_no_code(j,floor(i/4)+1) = -3+1i*d; 
            end
            d = d-2;   
          end
    end  
end

%% repetition
%*********************************QPSKs********************************%
QPSK_bits_rep=zeros(16,length(data_intrlvd_rep)/2);
binary_4 = [1 1; 0 1; 0 0; 1 0 ];
for j=1:16
    for i=1:2:length(data_intrlvd_rep)-1
        for k = 1 : 4
            if(data_intrlvd_rep(j,i:i+1) == binary_4(k,:) )
                 QPSK_bits_rep(j,floor(i/2)+1) =  sqrt(2)*exp(1i*((pi/4)+(pi/2)*(k-1))); 
            end
        end
    end
end

%*********************************16QAM********************************%
QAM16_bits_rep = zeros(16,length(data_intrlvd_rep)/4);
binary_16_3  = [1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];
binary_16_1  = [1 1 1 0; 1 1 1 1; 1 1 0 1; 1 1 0 0];
binary_16_n1 = [0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0];
binary_16_n3 = [0 0 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0];

for j=1:16
    for i = 1:4:length(data_intrlvd_rep)-3
        d = 3;
          for k = 1 : 4
            if(data_intrlvd_rep(j,i:i+3) == binary_16_3(k,:) )
                 QAM16_bits_rep(j,floor(i/4)+1) = 3+1i*d; 
            end
            if(data_intrlvd_rep(j,i:i+3) == binary_16_1(k,:) )
                 QAM16_bits_rep(j,floor(i/4)+1) = 1+1i*d; 
            end
            if(data_intrlvd_rep(j,i:i+3) == binary_16_n1(k,:) )
                 QAM16_bits_rep(j,floor(i/4)+1) = -1+1i*d; 
            end
            if(data_intrlvd_rep(j,i:i+3) == binary_16_n3(k,:) )
                 QAM16_bits_rep(j,floor(i/4)+1) = -3+1i*d; 
            end
            d = d-2;   
          end
    end  
end

%% channel
%% no code
%*********************************QPSK*********************************%
QPSK_received_no_code = zeros(16,length(bits_data));
Noisy_QPSK_no_code = zeros(16,length(QPSK_bits_no_code));
BER_QPSK_no_code = 0;
EB_QBSK_no_code = 1;

for d = 1:20
    %% channel
    No = EB_QBSK_no_code / (10^(d/10));
    v1 = randn(16,length(QPSK_bits_no_code));
    v2 = randn(16,length(QPSK_bits_no_code));
    v  = randn(16,length(QPSK_bits_no_code)).*sqrt(No/2)+1i*randn(16,length(QPSK_bits_no_code)).*sqrt(No/2);
    R  = sqrt(v1.^2 + v2.^2)/sqrt(2);
    Noisy_QPSK_no_code = R.*QPSK_bits_no_code + v;
    
    %% receiver
    Noisy_QPSK_no_code = Noisy_QPSK_no_code./R;
    for j=1:16
        k=1;    
        for g = 1:length(Noisy_QPSK_no_code)
            if(real(Noisy_QPSK_no_code(j,g)) >= 0)
                if(imag(Noisy_QPSK_no_code(j,g)) >= 0)
                    QPSK_received_no_code(j,k:k+1) = binary_4(1,:);
                elseif(imag(Noisy_QPSK_no_code(j,g)) <= 0)
                    QPSK_received_no_code(j,k:k+1) = binary_4(4,:);
                end
            elseif(real(Noisy_QPSK_no_code(j,g)) <= 0)  
                if(imag(Noisy_QPSK_no_code(j,g)) >= 0)
                    QPSK_received_no_code(j,k:k+1) = binary_4(2,:);
                elseif(imag(Noisy_QPSK_no_code(j,g)) <= 0)
                    QPSK_received_no_code(j,k:k+1) = binary_4(3,:);
                end
            end
        k = k+2;
        end
    end
        QPSK_received_no_code = matintrlv(QPSK_received_no_code,4,4);
        [NoisyQPSK,error] = symerr( QPSK_received_no_code,bits_data);
        BER_QPSK_no_code(d) = error;
end

%*********************************16QAM*********************************%
QAM16_received_no_code = zeros(16,length(bits_data));
Noisy_QAM16_no_code = zeros(16,length(QAM16_bits_no_code));
BER_QAM16_no_code = 0;
EB_QAM16_no_code = 2.5;

for d = 1:20
    %% channel
    No = EB_QAM16_no_code / (10^(d/10));
    v1 = randn(16,length(QAM16_bits_no_code));
    v2 = randn(16,length(QAM16_bits_no_code));
    v  = randn(16,length(QAM16_bits_no_code)).*sqrt(No/2)+1i*randn(16,length(QAM16_bits_no_code)).*sqrt(No/2);
    R  = sqrt(v1.^2 + v2.^2)/sqrt(2);
    Noisy_QAM16_no_code = R.*QAM16_bits_no_code + v;
    
    %% receiver
    Noisy_QAM16_no_code = Noisy_QAM16_no_code./R;
    for j=1:16
        k = 1;
        for g = 1:length(Noisy_QAM16_no_code)
            if(real(Noisy_QAM16_no_code(j,g)) > 2)
                a = 2;
                for r = 2:4
                    if(imag(Noisy_QAM16_no_code(j,g)) < a && imag(Noisy_QAM16_no_code(j,g)) > a-2)
                        QAM16_received_no_code(j,k:k+3) = binary_16_3(r,:); break;
                    elseif(imag(Noisy_QAM16_no_code(j,g)) > a )
                        QAM16_received_no_code(j,k:k+3) = binary_16_3(1,:); break;
                    end
                    a = a-2;
                end
            elseif(real(Noisy_QAM16_no_code(j,g)) > 0 && real(Noisy_QAM16_no_code(j,g)) < 2 )
                 a=2;
                for r = 2:4
                   if(imag(Noisy_QAM16_no_code(j,g)) < a && imag(Noisy_QAM16_no_code(j,g)) > a-2)
                        QAM16_received_no_code(j,k:k+3) = binary_16_1(r,:);break;
                    elseif(imag(Noisy_QAM16_no_code(j,g)) > a )
                        QAM16_received_no_code(j,k:k+3) = binary_16_1(1,:);break;
                    end
                    a=a-2;
                end
             elseif(real(Noisy_QAM16_no_code(j,g)) > -2 && real(Noisy_QAM16_no_code(j,g)) < 0 )
                 a=2;
                for r = 2:4
                    if(imag(Noisy_QAM16_no_code(j,g)) < a && imag(Noisy_QAM16_no_code(j,g)) > a-2)
                        QAM16_received_no_code(j,k:k+3) = binary_16_n1(r,:);break;
                    elseif(imag(Noisy_QAM16_no_code(j,g)) > a )
                        QAM16_received_no_code(j,k:k+3) = binary_16_n1(1,:);break;
                    end
                    a=a-2;
                end
                elseif(real(Noisy_QAM16_no_code(j,g)) < -2)
                a=2;
                for r = 2:4
                    if(imag(Noisy_QAM16_no_code(j,g)) < a && imag(Noisy_QAM16_no_code(j,g)) > a-2)
                        QAM16_received_no_code(j,k:k+3) = binary_16_n3(r,:);break;
                    elseif(imag(Noisy_QAM16_no_code(j,g)) > a )
                        QAM16_received_no_code(j,k:k+3) = binary_16_n3(1,:);break;
                    end
                    a=a-2;
                end   
            end
            k = k+4; 
        end
    end
        QAM16_received_no_code = matintrlv(QAM16_received_no_code,4,4);
        [NoisyQPSK,error] = symerr( QAM16_received_no_code,bits_data);
        BER_16QAM_no_code(d) = error;
end


%% repetition
%*********************************QPSK*********************************%
QPSK_received_rep = zeros(16,length(bits_data));
Noisy_QPSK_rep = zeros(16,length(QPSK_bits_rep));
BER_QPSK__rep = 0;
EB_QBSK__rep = 1;
QPSK_received_de_rep=zeros(16,length(bits_data));

for d = 1:20
    %% channel
    No = EB_QBSK__rep / (10^(d/10));
    v1 = randn(16,length(QPSK_bits_rep));
    v2 = randn(16,length(QPSK_bits_rep));
    v  = randn(16,length(QPSK_bits_rep)).*sqrt(No/2)+1i*randn(16,length(QPSK_bits_rep)).*sqrt(No/2);
    R  = sqrt(v1.^2 + v2.^2)/sqrt(2);
    Noisy_QPSK_rep = R.*QPSK_bits_rep + v;
    
    %% receiver
    Noisy_QPSK_rep = Noisy_QPSK_rep./R;
    for j=1:16
        k=1;    
        for g = 1:length(Noisy_QPSK_rep)
            if(real(Noisy_QPSK_rep(j,g)) >= 0)
                if(imag(Noisy_QPSK_rep(j,g)) >= 0)
                    QPSK_received_rep(j,k:k+1) = binary_4(1,:);
                elseif(imag(Noisy_QPSK_rep(j,g)) <= 0)
                    QPSK_received_rep(j,k:k+1) = binary_4(4,:);
                end
            elseif(real(Noisy_QPSK_rep(j,g)) <= 0)  
                if(imag(Noisy_QPSK_rep(j,g)) >= 0)
                    QPSK_received_rep(j,k:k+1) = binary_4(2,:);
                elseif(imag(Noisy_QPSK_rep(j,g)) <= 0)
                    QPSK_received_rep(j,k:k+1) = binary_4(3,:);
                end
            end
        k = k+2;
        end
 
       counter = 1;
       for i=1:length(bits_data)
           if sum(QPSK_received_rep(j,counter:counter+2)) > 1
                QPSK_received_de_rep(j,i) = 1;
           else
               QPSK_received_de_rep(j,i) = 0;
           end
       counter = counter + 3;
       end
    end
     QPSK_received_de_rep = matintrlv(QPSK_received_de_rep,4,4);
        [NoisyQPSK,error] = symerr( QPSK_received_de_rep,bits_data);
        BER_QPSK_rep(d) = error;
end

%*********************************16QAM*********************************%
QAM16_received_rep = zeros(16,length(bits_data));
Noisy_QAM16_rep = zeros(16,length(QAM16_bits_rep));
BER_QAM16_rep = 0;
EB_QAM16_rep = 2.5;
QAM16_received_de_rep=zeros(16,length(bits_data));

for d = 1:20
    %% channel
    No = EB_QBSK__rep / (10^(d/10));
    v1 = randn(16,length(QAM16_bits_rep));
    v2 = randn(16,length(QAM16_bits_rep));
    v  = randn(16,length(QAM16_bits_rep)).*sqrt(No/2)+1i*randn(16,length(QAM16_bits_rep)).*sqrt(No/2);
    R  = sqrt(v1.^2 + v2.^2)/sqrt(2);
    Noisy_QAM16_rep = R.*QAM16_bits_rep + v;
    
    %% receiver
    Noisy_QAM16_rep = Noisy_QAM16_rep./R;
    for j=1:16
        k = 1;
        for g = 1:length(Noisy_QAM16_rep)
            if(real(Noisy_QAM16_rep(j,g)) > 2)
                a = 2;
                for r = 2:4
                    if(imag(Noisy_QAM16_rep(j,g)) < a && imag(Noisy_QAM16_rep(j,g)) > a-2)
                        QAM16_received_rep(j,k:k+3) = binary_16_3(r,:); break;
                    elseif(imag(Noisy_QAM16_rep(j,g)) > a )
                        QAM16_received_rep(j,k:k+3) = binary_16_3(1,:); break;
                    end
                    a = a-2;
                end
            elseif(real(Noisy_QAM16_rep(j,g)) > 0 && real(Noisy_QAM16_rep(j,g)) < 2 )
                 a=2;
                for r = 2:4
                   if(imag(Noisy_QAM16_rep(j,g)) < a && imag(Noisy_QAM16_rep(j,g)) > a-2)
                        QAM16_received_rep(j,k:k+3) = binary_16_1(r,:);break;
                    elseif(imag(Noisy_QAM16_rep(j,g)) > a )
                        QAM16_received_rep(j,k:k+3) = binary_16_1(1,:);break;
                    end
                    a=a-2;
                end
             elseif(real(Noisy_QAM16_rep(j,g)) > -2 && real(Noisy_QAM16_rep(j,g)) < 0 )
                 a=2;
                for r = 2:4
                    if(imag(Noisy_QAM16_rep(j,g)) < a && imag(Noisy_QAM16_rep(j,g)) > a-2)
                        QAM16_received_rep(j,k:k+3) = binary_16_n1(r,:);break;
                    elseif(imag(Noisy_QAM16_rep(j,g)) > a )
                        QAM16_received_rep(j,k:k+3) = binary_16_n1(1,:);break;
                    end
                    a=a-2;
                end
                elseif(real(Noisy_QAM16_rep(j,g)) < -2)
                a=2;
                for r = 2:4
                    if(imag(Noisy_QAM16_rep(j,g)) < a && imag(Noisy_QAM16_rep(j,g)) > a-2)
                        QAM16_received_rep(j,k:k+3) = binary_16_n3(r,:);break;
                    elseif(imag(Noisy_QAM16_rep(j,g)) > a )
                        QAM16_received_rep(j,k:k+3) = binary_16_n3(1,:);break;
                    end
                    a=a-2;
                end   
            end
            k = k+4; 
        end
        
           counter = 1;
           for i=1:length(bits_data)
               if sum(QAM16_received_rep(j,counter:counter+2)) > 1
                    QAM16_received_de_rep(j,i) = 1;
               else
                   QAM16_received_de_rep(j,i) = 0;
               end
           counter = counter + 3;
           end      
    end
            QAM16_received_de_rep = matintrlv(QAM16_received_de_rep,4,4);
            [NoisyQPSK,error] = symerr( QAM16_received_de_rep,bits_data);
            BER_16QAM_rep(d) = error; 
end

x_axis=-4:15;
figure;
semilogy(x_axis,BER_QPSK_no_code,'r');
hold on
semilogy(x_axis,BER_QPSK_rep,'b');
hold on
legend('QPSK no code BER','QPSK repetition code BER');
xlabel('Eb/No');
ylabel('BER');
title('QPSK');

figure;
semilogy(x_axis,BER_16QAM_no_code,'r');
hold on
semilogy(x_axis,BER_16QAM_rep,'b');
hold on
legend('16QAM no code BER','16QAM repetition code BER');
xlabel('Eb/No');
ylabel('BER');
title('16QAM');

figure;
semilogy(x_axis,BER_QPSK_no_code,'r');
hold on
semilogy(x_axis,BER_16QAM_no_code,'g');
hold on

semilogy(x_axis,BER_QPSK_rep,'b');
hold on
semilogy(x_axis,BER_16QAM_rep,'k');
hold on
legend('QPSK no code BER','16QAM no code BER','QPSK repetition code BER','16QAM repetition code BER');
xlabel('Eb/No');
ylabel('BER');
title('QPSK & 16QAM');
