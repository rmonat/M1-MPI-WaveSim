#include <stdio.h>
#include <stdlib.h>

int main()
{
    size_t n = 256, m = 3*n;
    char *file = "young.in";
    double speed = 0.35;
    char type = 1;


    FILE *fp;
    fp = fopen(file, "wb");

    // header of the file
    fwrite(&type, sizeof(char), 1, fp);
    fwrite(&n, sizeof(size_t), 1, fp);
    fwrite(&m, sizeof(size_t), 1, fp);
    fwrite(&speed, sizeof(double), 1, fp);
    
    char cell_type;
    double cell_value;

    // now the real part
    for(size_t i = 0; i < n; i++)
    {
	for(size_t j = 0; j < m; j++)
	{
//	    if(j == 0) // walls at the left
//	    {
//		cell_type = 1;
//		cell_value = 0;
//	    }
//	    else if((i == 0 || i == (n-1) && j <= m/3)) // walls on top and bottom and beginning
//	    {
//		cell_type = 1;
//		cell_value = 0;
//	    }
//	    else
		if(j == m/3 && !(i >= n/3-2 && i <= n/3+2) && !(i >= 2*n/3-2 && i <= 2*n/3+2)) // walls except at two points
	    {
		cell_type = 1;
		cell_value = 0;
	    }
	    else if(j == (m/2-1)) // sensors at the end
	    {
		cell_type = 2;
		cell_value = 0;
	    }
	    else if((j >= m/6 && j <= m/6+8) && (i >= (n/2-1-4) && i <= (n/2+4))) // source
	    {
		cell_type = 0;
		cell_value = 3;
	    }
	    else // Basic case
	    {
		cell_type = 0;
		cell_value = 0;
	    }
	    fwrite(&cell_type, sizeof(char), 1, fp);
	    fwrite(&cell_value, sizeof(double), 1, fp);
	}
    }
    fclose(fp);
}
