void ftocc(int *nin, double *x)
{
int n = nin[0];
int i;
for (i=0; i<n; i++)
x[i] = (x[i] - 32) * 5 / 9;
}