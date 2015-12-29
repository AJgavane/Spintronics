#include<stdio.h>
#include<stdlib.h>

void computeTsqTcub(char *, char *);

int main()
{
	char *readfile = malloc(80*sizeof(char));
	char *writefile = malloc(100*sizeof(char));
	printf("Name the file whose switch time values you wish to\nsquare and cube.\n");
	scanf("%s", readfile);
	printf("Enter the name of the writefile\n");
	scanf("%s", writefile);
	
	
	
	computeTsqTcub(readfile, writefile);
	printf(".....................Done!\n");

}

void computeTsqTcub(char *readfile, char *writefile)
{
	
	FILE *fp1 = fopen(readfile, "r");
	if(fp1==NULL)
	{
		printf("No such file\n");
		return;
	}
	FILE *fp2 = fopen(writefile, "w");
	int ind; double ti;

	while(!feof(fp1))
	{
		fscanf(fp1, "%d%lf", &ind, &ti);
		fprintf(fp2, "%d\t%.4f\t%.4f\t%.4f\n", ind, ti, ti*ti, ti*ti*ti);
	}
	
	fclose(fp1); fclose(fp2);
	return;

}
