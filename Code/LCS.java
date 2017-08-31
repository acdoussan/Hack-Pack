// Ross 04/30/2017

public class LCS
{
	static int LCS(String str1, String str2, int len1, int len2)
	{
		int[][] memo = new int[len1+1][len2+1];

		for(int i=1; i<=len1; i++)
		{
			for(int j=1; j<=len2; j++)
			{
				if((i == 0) || (j == 0))
				{
					memo[i][j] = 0;
				}
				else if(str1.charAt(i-1) == str2.charAt(j-1))
				{
					memo[i][j] = 1 + memo[i-1][j-1];
				}
				else
				{
					int sub = Math.max(memo[i][j-1], memo[i-1][j]);
					memo[i][j] = sub;
				}
			}
		}
		return memo[len1][len2];
	}
}
