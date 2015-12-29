#include<stdio.h>
#include<stdlib.h>

void tavg_sigmat(char *, char *);

int main()
{
	char *readfile = malloc(80*sizeof(char));
	char *writefile = malloc(100*sizeof(char));
	printf("Name the file whose switch time values you wish to\nsquare and cube.\n");
	scanf("%s", readfile);
	printf("Enter the name of the writefile\n");
	scanf("%s", writefile);
	
	
	
	tavg_sigmat(readfile, writefile);
	printf(".....................Done!\n");

}

void tavg_sigmat(char *readfile, char *writefile)
{
	FILE *fp1 = fopen(readfile, "r");
	
	FILE *fp2 = fopen(writefile, "w");
	
	double tsavg=0, tsqavg=0, tcubavg=0; int no;
	double mtsavg, mtsqavg, mtcubavg;
	
	while(!feof(fp1))
	{
		fscanf(fp1, "%d%lf%lf%lf", &no, &mtsavg, &mtsqavg, &mtcubavg);
		tsavg+=mtsavg;
		tsqavg+=mtsqavg;
		tcubavg+=mtcubavg;
	}

	tsavg = tsavg/no;
	tsqavg = tsqavg/no;
	tcubavg = tcubavg/no;

	double sigma_t = tsqavg - tsavg*tsavg;

	double skewness_t = tcubavg - 3*sigma_t*tsavg + 2*tsavg*tsavg*tsavg;

	fprintf(fp2, "%.4f\t%.4f\t%.4f", tsavg, sigma_t, skewness_t);
	
	fclose(fp1);
	fclose(fp2);


}
