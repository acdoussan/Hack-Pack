// Adam Doussan AD844156 04/13/2017

public class Dijkstras
{
	final public static int oo = (int) 1e9;
	private static int n;

	// pass in a value of true to stop if you want early termination for
	// destination d, false if you want the entire dist array to be filled in
	public static int[] run(ArrayList<Edge>[] g, int s, int d, boolean stop)
	{
		n = g.length;

		boolean [] visited = new boolean[n];
		int [] dist = new int [n];
		Arrays.fill(dist, oo);

		PriorityQueue<Edge> pq = new PriorityQueue<Edge>();
		pq.add(new Edge(s, 0));
		dist[s] = 0;

		while(!pq.isEmpty())
		{
			Edge at = pq.remove();
			if(visited[at.e]) continue;
			visited[at.e] = true;

			if(stop && at.e == d) return dist;

			for(Edge adj : g[at.e])
				if(!visited[adj.e] && adj.w + at.w < dist[adj.e])
					pq.add(new Edge(adj.e, dist[adj.e] = adj.w + at.w));
		}

		return dist;
	}
}

class Edge implements Comparable<Edge>
{
	int e, w;

	public Edge(int e, int w)
	{
		this.e = e; this.w = w;
	}

	public int compareTo(Edge o)
	{
		return this.w - o.w;
	}
}
