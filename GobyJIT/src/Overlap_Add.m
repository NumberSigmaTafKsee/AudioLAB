%%
%Programa que calcula o Overlap entre dois vetores
%%
%Descrição
%Foram utilizados dois vetores bases de valores sequenciais para visualisar
%melhor a resposta do programa, mas o programa foi feito de forma a
%calcular o método do Overlap Add para quaisquer dois vetores

clc 
clear all

%Vetor Entrada, representando a entrada
Entrada = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

%Vetor H representando a resposta
H = [ 0.5 0.5]; 

%Tamanho da resposta da convolução
N = size(Entrada,2) +size(H,2) - 1;

%Vetor da resposta da convolucao
Resposta = zeros(1,N);

%Fator de repartição do Overlapp
DivisaoOverlap= 2;
%%
%Preenchimento com zeros da entrada e da resposta

    Entradazeros = zeros(N,1);
   for contador = 1:size(Entrada,2)
        Entradazeros(contador) = Entrada(contador);
   end
   Entrada=Entradazeros;

   Hzero = zeros(N,1);
   for contador = 1:size(H,2)
        Hzero(contador) = H(contador);
   end
   H=Hzero;

   
%%Overlap Add   

%Convolucao Circular
    for indice = 1:N
        for  cont = 1:N
                Resposta(1,indice) = Entrada(cont,1)*H(mod(N+indice-cont,N)+1)+Resposta(1,indice);
        end
    end

for indice = 1:DivisaoOverlap -1
    Resposta(1,floor(N*indice/DivisaoOverlap)) = 2*Resposta(1,floor(N/DivisaoOverlap));
end

%Representação gráfica para entrada e saída
figure(1)
stem(Entrada)

figure(2)
stem(H)

figure(3)
stem(Resposta)