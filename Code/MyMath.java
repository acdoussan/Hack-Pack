// Adam Doussan AD844156 04/13/2017

public class MyMath
{
	public static long lcm(long a, long b)
	{
		return a / gcd(a,b) * b;
	}

	public static long gcd(long a, long b)
	{
		return (b == 0) ? a : gcd(b, a%b);
	}

	public static int MCSS(int [] days)
	{
		int maxSum = days[0];
		int runningSum = days[0];

		for(int i = 0; i < days.length; i++)
		{
			if(runningSum < 0)
				runningSum = 0;

			runningSum += days[i];

			if(runningSum > maxSum)
				maxSum = runningSum;
		}

		return maxSum;
	}

	public static boolean[] primeSieve(int n)
	{
		boolean[] isPrime = new boolean[n+1];
		Arrays.fill(isPrime, true);
		isPrime[0] = false;
		isPrime[1] = false;

		for(int i = 2; i*i <= n; i++)
			if(isPrime[i])
				for(int j = i + i; j <= n; j += i)
					isPrime[j] = false;

		return isPrime;
	}

	public static long[] extendedEuclid(long a, long b)
	{
		if(b == 0)
			return new long []{1,0,a};
		else
		{
			long q = a / b;
			long r = a % b;
			long [] rec = extendedEuclid(b,r);
			return new long [] {rec[1], rec[0] - q*rec[1], rec[2]};
		}
	}

	public static BigInteger choose(long x, long y)
	{
		if(y<0 || y>x)
			return BigInteger.ZERO;
		if(y == 0 || y == x)
			return BigInteger.ONE;

		BigInteger ans = BigInteger.ONE;

		for (long i = x - y + 1; i <= x; i++)
		    ans = ans.multiply(BigInteger.valueOf(i));
		for (long j = 1; j <= y; j++)
		    ans = ans.divide(BigInteger.valueOf(j));
		return ans;
	}

	public static long [][] chooseDP(int maxDepth)
	{
		long[][] choose = new long[maxDepth][maxDepth];

		for(int i = 0; i < maxDepth; i++)
		{
			choose[i][0] = choose[i][i] = 1;

			for(int j = 1; j < i; j++) 
				choose[i][j] = choose[i-1][j] + choose[i-1][j-1];
		}

		return choose;
	}

	public static BigInteger sumOfDivisors(ArrayList<pair> primeFac)
	{
		BigInteger sum = new BigInteger("1");

		for(pair a : primeFac)
		{
			BigInteger prime = new BigInteger("" + a.prime);

			BigInteger temp = prime.pow((int)a.exp + 1).subtract(BigInteger.ONE);
			temp = temp.divide(prime.subtract(BigInteger.ONE));

			sum = sum.multiply(temp);
		}

		return sum;
	}

	public static ArrayList<pair> primeFactorize(int n)
	{
		ArrayList<pair> res = new ArrayList<>();
		int div = 2;

		while(div*div <= n)
		{
			int exp = 0;
			while(n % div == 0)
			{
				n /= div;
				exp++;
			}

			if(exp > 0) res.add(new pair(div, exp));
			div++;
		}

		if(n > 1) res.add(new pair(n,1));
		return res;
	}
}

class pair
{
	public int prime, exp;

	public pair(int p, int e)
	{
		prime = p; exp = e;
	}
}
