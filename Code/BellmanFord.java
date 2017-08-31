// Adam Doussan AD844156 04/13/2017

public class BellmanFord
{
	final public static int oo = (int) 1e9;

	// returns null if there is a negative cycle, edges are directional
	public static int[] run(Edge[] eList, int n, int s)
	{
		int [] dist = new int [n];
		Arrays.fill(dist, oo);
		dist[s] = 0;

		for(int i = 0; i < n - 1; i++)
			for(Edge e : eList)
				if(dist[e.v1] + e.w < dist[e.v2])
					dist[e.v2] = dist[e.v1] + e.w;

		for(Edge e : eList)
			if(dist[e.v1] + e.w < dist[e.v2])
				return null;

		return dist;
	}
}

class Edge
{
	public int v1, v2, w;

	public Edge(int v1, int v2, int w)
	{
		this.v1 = v1; this.v2 = v2; this.w = w;
	}
}
