// Adam Doussan AD844156 04/13/2017

public class Floyd
{
	public static int [][] path;
	final public static int oo = (int)1e9;

	// gets answer without destroying passed in array
	public static int [][] run(int [][] adj)
	{
		int n = adj.length;
		int [][] m = copy(adj);
		path = new int[adj.length][adj.length];

		for(int i = 0; i < adj.length; i++)
			for(int j = 0; j < adj.length; j++)
				if(adj[i][j] != oo)
					path[i][j] = j;
				else
					path[i][j] = oo;

		for(int k = 0; k < n; k++)
			for(int i = 0; i < n; i++)
				for(int j = 0; j < n; j++)
					if(m[i][k] != oo && m[k][j] != oo)
						if(m[i][k] + m[k][j] < m[i][j])
						{
							m[i][j] = m[i][k] + m[k][j];
							path[i][j] = path[i][k];
						}

		return m;
	}

	public static int [][] copy(int [][] a)
	{
		int [][] res = new int [a.length][a[0].length];
			for(int i = 0; i < a.length; i++)
				for(int j = 0; j < a[0].length; j++)
					res[i][j] = a[i][j];
		return res;
	}

	public static ArrayList<Integer> getPath(int i, int j)
	{
		if(path[i][j] == oo)
			return null;

		ArrayList<Integer> myPath = new ArrayList<>();

		myPath.add(i);

		while(i != j)
		{
			i = path[i][j];
			myPath.add(i);
		}

		return myPath;
	}
}
