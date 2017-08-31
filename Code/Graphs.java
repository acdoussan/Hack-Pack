// Adam Doussan AD844156 04/11/2017

public class Graphs
{
	public static void dfs(ArrayList[] graph, boolean [] visited, int v)
	{
		visited[v] = true;
		for(Integer next : ((ArrayList<Integer>[])graph)[v])
			if(!visited[next])
				dfs(graph, visited, next);
	}

	public static int[] bfs(ArrayList[] graph, int v)
	{
		int n = graph.length;
		int [] distance = new int [n];

		Arrays.fill(distance, -1);
		distance[n] = 0;

		Queue<Integer> q = new ArrayDeque<Integer>();
		q.add(v);

		while(!q.isEmpty())
		{
			int cur = q.remove();
			for(Integer next : ((ArrayList<Integer>[])graph)[cur])
			{
				if(distance[next] == -1)
				{
					distance[next] = distance[cur] + 1;
					q.add(next);
				}
			}
		}

		return distance;
	}
}
