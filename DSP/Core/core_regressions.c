
#define LEN(a) (sizeof(a) / sizeof(*a))

double a;
double b;
double c;

double sum( double arr[], int n) {
   if (n < 0) {
     return 0;
   } 
   else{
     return arr[n] + sum(arr, n-1);
    }
}

double mean(double arr[], int n){
	double avg;
	avg = sum( arr, n-1 )/ n;
	return avg;
}

double standard_deviation(double arr[], double mean, int n){
	double deviation_score[n];
	double deviation_score_squared[n];
	for(int i = 0; i<n; i++){
		deviation_score[i] = arr[i]-mean;
		deviation_score_squared[i] = deviation_score[i]*deviation_score[i];
	}
	
	double stdev;
	stdev = sqrt(sum(deviation_score_squared, n-1)/(n-1));
	return stdev;
}

int linear_regression(double x[], double y[], int n){
	double mean_x = mean(x, n);
	double mean_y = mean(y, n);
	double deviation_score_x[n];
	double deviation_score_y[n];
	double deviation_score_squared_x[n];
	double deviation_score_squared_y[n];
	double product[n];
	
	for(int i = 0; i<n; i++){
		deviation_score_x[i] = x[i]-mean_x;
		deviation_score_y[i] = y[i]-mean_y;
		deviation_score_squared_x[i] = deviation_score_x[i]*deviation_score_x[i];
		deviation_score_squared_y[i] = deviation_score_y[i]*deviation_score_y[i];
		product[i] =  deviation_score_x[i]*deviation_score_y[i];
	}
	
	double r;
	r = sum(product,n-1)/sqrt(sum(deviation_score_squared_x,n-1)*sum(deviation_score_squared_y,n-1));
	
	b = r*standard_deviation(y,mean_y, n)/standard_deviation(x,mean_x, n);
	a = mean_y - b*mean_x;
	
	printf("Values of regression coefficients are: a = %0.2lf and b = %0.2lf",a,b);
    printf("\nEquation of best fit is: y = %0.2lf + %0.2lfx",a,b);
    
    return 0;
}

int quadratic_regression(double x[], double y[], int n){
	double x2[n];
	double x3[n];
	double x4[n];
	double xy[n];
	double x2y[n];
	
	for(int i = 0; i<n; i++){
		x2[i] = x[i]*x[i];
		x3[i] = x2[i]*x[i];
		x4[i] = x2[i]*x2[i];
		xy[i] = x[i]*y[i];
		x2y[i] = x2[i]*y[i];
	}
	
	double x_x = sum(x2,n-1) - (sum(x,n-1)*sum(x,n-1))/n;
	double x_y = sum(xy,n-1) - (sum(x,n-1)*sum(y,n-1))/n;
	double x_x2 = sum(x3,n-1) - (sum(x,n-1)*sum(x2,n-1))/n;
	double x2_y = sum(x2y,n-1) - (sum(x2,n-1)*sum(y,n-1))/n;
	double x2_x2 = sum(x4,n-1) - (sum(x2,n-1)*sum(x2,n-1))/n;
	
	a = (x2_y*x_x - x_y*x_x2)/(x_x*x2_x2 - x_x2*x_x2);
	b = (x_y*x2_x2 - x2_y*x_x2)/(x_x*x2_x2 - x_x2*x_x2);
	c = mean(y,n) - b*mean(x,n) - a*mean(x2,n);
	
	printf("Values of regression coefficients are: a = %0.2lf, b = %0.2lf and c = %0.2lf",a,b,c);
    printf("\nEquation of best fit is: y = %0.2lfx^2 + %0.2lfx + %0.2lf",a,b,c);
    
    return 0;

}

int exponential_regression(double x[], double y[], int n){
	double z[n];
	double x_z[n];
	double x_squared[n];
	
	for(int i = 0; i<n; i++){
		z[i] = log(y[i]);
		x_z[i] = x[i]*z[i];
		x_squared[i] = x[i]*x[i];
		
	}
	
	double a_1 = (n*sum(x_z,n-1) - sum(x,n-1)*sum(z,n-1))/(n*sum(x_squared,n-1) - sum(x,n-1)*sum(x,n-1));
	double a_0 = (mean(z,n)) - a_1*(mean(x,n));
	
	a = exp(a_0);
	b = a_1;
	
	printf("Values of regression coefficients are: a = %0.5lf and b = %0.5lf",a,b);
    printf("\nEquation of best fit is: y = %0.5lfe^(%0.5lfx)",a,b);
    
    return 0;
}