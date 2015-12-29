#include<stdio.h>

int main()
{
	int n, i;
	printf("Enter n\n");
	scanf("%d", &n);
	int a[n];
	srand(time(NULL));
	for(i=0; i<n; i++)
	{
		a[i] = rand()%100;
		printf("%d\t", a[i]);
	}

}
