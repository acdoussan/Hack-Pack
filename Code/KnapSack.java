// Adam Doussan AD844156 04/14/2017

public class KnapSack
{
	// unlim = true means can take unlimited times, otherwise each thing can be
	// taken only once. Add multiple times if multiple limited number of same
	// object. kp[money] contains the answer for money.
	public static long[] solve(int [] cost, int [] value, int money, boolean unlim)
	{
		long [] kp = new long [money+1];
		int n = cost.length;

		if(unlim)
		{
			for(int j = 0; j < n; j++)
			{
				for(int k = cost[j]; k < money+1; k++)
				{
					if(kp[k-cost[j]] + value[j] > kp[k])
						kp[k] = kp[k - cost[j]] + value[j];
				}
			}
		}
		
		else
		{
			for(int j = 0; j < n; j++)
			{
				for(int k = money+1; k >= cost[j] ; k--)
				{
					if(kp[k-cost[j]] + value[j] > kp[k])
						kp[k] = kp[k - cost[j]] + value[j];
				}
			}
		}

		return kp;
	}
}
