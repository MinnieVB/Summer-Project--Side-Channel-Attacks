#include <stdio.h>
#include <math.h>

unsigned long mean(unsigned long array[], unsigned long limit)
{
	unsigned long sum = 0, average;
	int i;
	for (i=0; i<limit; i++)
	{
		sum = sum + array[i];
	}
	average = sum/limit;
	return average;
}


unsigned long sd(unsigned long array[], unsigned long limit, unsigned long mean)
{
	int i;
	unsigned long sum = 0, sd;
	for (i=0; i<limit; i++)
	{
		sum = sum + (array[i]-mean)*(array[i]-mean);
	}
	sd = (unsigned long)sqrt((double)sum/limit);
	return sd;
}


unsigned long min1(unsigned long array[], unsigned long limit)
{
	int i;
	unsigned long min = array[0];
	
	for (i=1; i<limit; i++)
	{
		if (array[i] < min)
		{
			min = array[i];
		}
	}
	
	return min;
}


unsigned long min2(unsigned long array[], unsigned long limit)
{
	int i, min1Index = 0;
	unsigned long min1 = array[0], min2;
	
	//Finding first minimum
	for (i=1; i<limit; i++)
	{
		if (array[i] < min1)
		{
			min1 = array[i];
			min1Index = i;
		}
	}
	array[min1Index] = 0;
	
	//Finding second minimum
	if (min1Index == 0)
	{
		min2 = array[1];
	}
	else
	{
		min2 = array[0];
	}
	
	for (i=1; i<limit; i++)
	{
		if (i != min1Index && array[i] < min2)
		{
			min2 = array[i];
		}
	}

	return min2;
}

int main()
{
	unsigned long start, add[1000], addMean, addSD, addMin1, addMin2;
	int numAdd = 0;
	
	add[0] = 7;
	numAdd++;
	add[1] = 2;
	numAdd++;
	add[2] = 4;
	numAdd++;
	add[3] = 2;
	numAdd++;
	add[4] = 3;
	numAdd++;
	add[5] = 8;
	numAdd++;
	
	addMean = mean(add, numAdd);
	addSD = sd(add, numAdd, addMean);
	addMin1 = min1(add, numAdd);
	addMin2 = min2(add, numAdd);
	
	printf("& %lu ", addMean);
	printf("& %lu ", addSD);
	printf("& %lu ", addMin1);
	printf("& %lu ", addMin2);
	printf("\\\\ \n");
	return 0;

}
