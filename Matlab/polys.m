%a program to compare original data curve with n degree polynomial curve
%preparing data
cp_data=load ('data'); %data file values of temp vs cp
temp=cp_data(:,1); %calling temperature data points
cp=cp_data(:,2); %calling specific heat data points

%curve fit
%linear polynomial - cp=a*T + b
coeff1=polyfit(temp,cp,1);
predicted_cp1=polyval(coeff1,temp);

%second degree polynomial - cp=a*T^2 + b*T + c
coeff2=polyfit(temp,cp,2);
predicted_cp2=polyval(coeff2,temp);

%third degree polynomial - cp=a*T^3 + b*T^2 + c*T + d
coeff3=polyfit(temp,cp,3);
predicted_cp3=polyval(coeff3,temp);

% plotting basic data curve
figure(1)
plot(temp,cp,'linewidth',2)
xlabel("Temperature (K)")
ylabel("Specific Heat (Kj/Kmol-K)")
title("Basic Curve")
hold off

%compare actual with predicted curves
%plotting linear curve with basic curve
figure(2)
plot(temp,cp,'linewidth',2);
hold on
plot(temp,predicted_cp1,'linewidth',2,'color','m')
legend('Original data','Linear data','location','northwest')
xlabel("Temperature (K)")
ylabel("Specific Heat (Kj/Kmol-K)")
title("Basic curve vs Linear data curve")
hold off

%plotting second degree polynomial curve with basic curve
figure(3)
plot(temp,cp,'linewidth',2);
hold on
plot(temp,predicted_cp2,'linewidth',2,'color','g')
legend('Original data','Second degree data','location','northwest')
xlabel("Temperature (K)")
ylabel("Specific Heat (Kj/Kmol-K)")
title("Basic curve vs Second degree data curve")
hold off

%plotting third degree polynomial curve with basic curve
figure(4)
plot(temp,cp,'linewidth',2);
hold on
plot(temp,predicted_cp3,'linewidth',2,'color','r')
legend('Original data','Third degree data','location','northwest')
xlabel("Temperature (K)")
ylabel("Specific Heat (Kj/Kmol-K)")
title("Basic curve vs Third degree data curve")
