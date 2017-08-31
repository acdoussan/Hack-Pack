// Adam Doussan AD844156 05/01/2017
import java.util.*;
public class myFordFulkerson
{
	public int n, s, t;
	public int [][] g;

	public myFordFulkerson(int n)
	{
		this.n = n+2; this.s = n; this.t = n+1;
		g = new int [n+2][n+2];
	}

	public void addEdge(int v1, int v2, int cap)
	{
		g[v1][v2] = cap;
	}

	public int flow()
	{
		int ans = 0;
		int [] prev = new int [n];

		while(bfs(prev))
		{
			int myNext = prev[t];
			int myFlow = g[myNext][t];

			while(myNext != s)
			{
				myFlow = Math.min(myFlow, g[prev[myNext]][myNext]);
				myNext = prev[myNext];
			}

			ans += myFlow;
			myNext = t;

			while(myNext != s)
			{
				int temp = prev[myNext];
				g[myNext][temp] += myFlow;
				g[temp][myNext] -= myFlow;
				myNext = temp;
			}
		}

		return ans;
	}

	public boolean bfs(int[] prev)
	{
		Queue<Integer> q = new ArrayDeque<>();
		boolean [] visited = new boolean [n];

		q.add(s);
		visited[s] = true;

		while(!q.isEmpty())
		{
			int myNext = q.remove();

			for(int i = 0; i < n; i++)
			{
				if(!visited[i] && g[myNext][i] > 0)
				{
					visited[i] = true;
					q.add(i);
					prev[i] = myNext;
				}
			}
		}

		return visited[t];
	}
}
