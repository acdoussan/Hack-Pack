// Ross 04/30/2017

public class MatrixChain
{
	public static void main(String [] args)
	{
		//Matrix A is 2x4, Matrix B is 4x2, Matrix C is 2x3 and so on.
		int matr[] = new int[] {2, 4, 2, 3, 1, 4};
		int size = matr.length;

		int res = MCM(matr, size);
		System.out.println(res);
	}

	static int MCM(int[] matr, int size)
	{
		int a;
		int memo[][] = new int[size][size];

		for(int i=0; i<size; i++)
		{
			memo[i][i] = 0;
		}

		for(int sub=2; sub <= (size-1); sub++)
		{
			for(int j=1; j <= ((size-1) - sub + 1); j++)
			{
				memo[j][sub+j-1] = Integer.MAX_VALUE;
				for(int k=j; k <= (sub+j-2); k++)
				{
					a = memo[j][k] + memo[k+1][j+sub-1] + (matr[j-1] * matr[k] * matr[sub+j-1]);
					if (a < memo[j][j+sub-1])
					{
						memo[j][j+sub-1] = a;
					}
				}
			}
		}

		return memo[1][size-1];
	}
}
